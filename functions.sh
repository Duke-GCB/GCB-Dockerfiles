#!/bin/bash

function current_branch_name() {
  git rev-parse --abbrev-ref HEAD
}

function get_num_parents() {
  sha=$1
  parents=$(git cat-file -p $sha | grep parent | wc -l)
  echo $parents
}

# Given a range, produce the list of file paths changed
function changed_paths_in_range() {
  compare_range=$1
  # diff-filter=d excludes deleted files
  git diff --name-only --diff-filter=d $compare_range
}

# Check if DEPLOY_BRANCH is set and current branch can be tested
function check_deploy_branch() {
  if [[ "$DEPLOY_BRANCH" == "" ]]; then
    echo "Error: DEPLOY_BRANCH is empty"
    echo "Please ensure DEPLOY_BRANCH is set to the name of the git branch that should be considered for deploying, typically 'master'";
    exit 1;
  else
    echo "Using Deploy branch as $DEPLOY_BRANCH..."
  fi
  current_branch=$(current_branch_name)
  if [[ "$current_branch" == "$DEPLOY_BRANCH" ]]; then
    # Current branch is the deploy branch, we must have a merge commit to test/build
    num_parents=$(get_num_parents "HEAD")
    if [[ "$num_parents" -lt "2" ]]; then
      echo "Error: Current branch is the deploy branch ($DEPLOY_BRANCH), but HEAD is not a merge commit. Build failed." && exit 1
    fi
  fi
}

# If the current branch is the deploy branch, return a range representing
# the two parents of the HEAD's merge commit. If not, return a range comparing
# the current HEAD with the deploy_branch
function get_compare_range() {
  current_branch=$(current_branch_name)
  if [[ "$current_branch" == "$DEPLOY_BRANCH" ]]; then
      # On the deploy branch (e.g. master)
      # check_deploy_branch should have verified it is a merge commit
      range_start="HEAD^1" # alias for first parent
      range_end="HEAD^2" # alias for second parent
  else
    # Not on the deploy branch (e.g. master)
    # When not on the deploy branch, always compare with the deploy branch
    # Circle resets master to the tested commit, so we have to use origin/master
    range_start="origin/$DEPLOY_BRANCH"
    range_end="HEAD"
  fi
  echo "$range_start $range_end"
}

# Given a docker repo owner, image name, and version, produce a local docker build command
function build_docker_cmd() {
  owner=$1
  tool=$2
  version=$3
  echo "docker build -f $tool/$version/Dockerfile" -t "$owner/$tool:$version" "$tool/$version"
}

# Given a docker repo owner, image name, and version, produce a docker push command
function push_docker_cmd() {
  owner=$1
  tool=$2
  version=$3
  echo "docker push $owner/$tool:$version";
}

# Given a docker repo owner, image name, and version, produce a docker pull command
function pull_docker_cmd() {
  owner=$1
  tool=$2
  version=$3
  echo "docker pull $owner/$tool:$version";
}

# Given a docker repo owner, image name, source and dest tags produce a docker tag command
function tag_docker_cmd() {
  owner=$1
  tool=$2
  src=$3
  tag=$4
  echo "docker tag $owner/$tool:$version $owner/$tool:$tag"
}

# Given a docker repo owner, image name, and version, produce a command that returns the image id if it exists locally
function docker_image_id_cmd() {
  owner=$1
  tool=$2
  tag=$3
  echo "docker images $owner/$tool:$tag -q"
}

# Given a docker repo owner, image name, and version, check if it exists locally and pull if necessary
function ensure_local_image() {
  owner=$1
  tool=$2
  version=$3
  local_image_id=$($(docker_image_id_cmd $owner $tool $version))
  if [[ "$local_image_id" == "" ]]; then
    echo "Image $owner/$tool:$version does not exist locally for tagging, pulling..."
    $(pull_docker_cmd $owner $tool $version)
  fi
}

# Given
# 1. a Docker repo owner (e.g. "dukegcb") and
# 2. a list of relative paths to Dockerfiles (e.g. "fastqc/0.11.4/Dockerfile bwa/0.7.12/Dockerfile",
# issue a docker build command and tag any versions with a latest symlink
function build_images() {
  echo "Building changed Dockerfiles..."; echo
  owner="$1"
  changed_paths="$2"
  # Check for Dockerfile changes first
  for changed_path in $changed_paths; do
    IFS='/' read -r -a f <<< "$changed_path"
    tool="${f[0]}"
    version="${f[1]}"
    filename="${f[2]}"
    if [[ "$filename" == "Dockerfile" && "$version" != "latest" ]]; then
      attempted_build="1"
      echo "Building $owner/$tool:$version..."
      $(build_docker_cmd $owner $tool $version)
      # Check if there's a symlink $tool/latest pointing to THIS version
      if [[ "$tool/latest/Dockerfile" -ef "$tool/$version/Dockerfile" ]]; then
        echo "Tagging $owner/$tool:$version as $owner/$tool:latest"
        $(tag_docker_cmd $owner $tool $version "latest")
      fi
    fi
  done;

  # After building all Dockerfiles, check for any changes to latest
  echo "Updating latest tags..."; echo
  for changed_path in $changed_paths; do
    IFS='/' read -r -a f <<< "$changed_path"
    tool="${f[0]}"
    version="${f[1]}"
    if [[ -L "$changed_path" && "$filename" == "" && "$version" == "latest" ]]; then
      attempted_build="1"
      # The changed file is a symlink called latest, e.g. "fastqc/latest"
      # Determine the version it's pointing to
      dest_version=$(readlink $changed_path)
      # In order to tag to version, it must exist locally. If it wasn't built in previous loop,
      # need to pull it
      ensure_local_image $owner $tool $dest_version
      echo "Tagging $owner/$tool:$dest_version as $owner/$tool:latest"
      $(tag_docker_cmd $owner $tool $dest_version "latest")
    fi
  done;

  if [[ "$attempted_build" == "" ]]; then
    echo "No changes to Dockerfiles or latest symlinks detected, nothing to build";
  fi
}

# Given
# 1. a Docker repo owner (e.g. "dukegcb") and
# 2. a list of relative paths to Dockerfiles (e.g. "fastqc/0.11.4/Dockerfile bwa/0.7.12/Dockerfile",
# issue a docker push command for the images built by build_images
function push_images() {
  owner="$1"
  changed_paths="$2"
  for changed_path in $changed_paths; do
    IFS='/' read -r -a f <<< "$changed_path"
    tool="${f[0]}"
    version="${f[1]}"
    filename="${f[2]}"
    if [[ "$filename" == "Dockerfile" && "$version" != "latest" ]]; then
      attempted_push="1"
      echo "Pushing $owner/$tool:$version..."
      $(push_docker_cmd $owner $tool $version)
      # Check if there's a symlink $tool/latest pointing to THIS version
      if [[ "$tool/latest/Dockerfile" -ef "$tool/$version/Dockerfile" ]]; then
        echo "Pushing $owner/$tool:latest..."
        $(push_docker_cmd $owner $tool "latest")
      fi
    fi
  done;

  # After pushing all Dockerfiles, check for any changes to latest and push those
  echo "Pushing latest tags..."; echo
  for changed_path in $changed_paths; do
    IFS='/' read -r -a f <<< "$changed_path"
    tool="${f[0]}"
    version="${f[1]}"
    if [[ -L "$changed_path" && "$filename" == "" && "$version" == "latest" ]]; then
      attempted_push="1"
      # The changed file is a symlink called latest, e.g. "fastqc/latest"
      # Determine the version it's pointing to
      echo "Pushing $owner/$tool:latest..."
      $(push_docker_cmd $owner $tool "latest")
    fi
  done;

  if [[ "$attempted_push" == "" ]]; then
    echo "No changes to Dockerfiles or latest symlinks detected, nothing to push";
  fi
}

function print_changed() {
  range="$1"
  paths="$2"
  echo "Changed files in ($range)"
  echo
  for changed_path in $paths; do
    echo "  $changed_path"
  done
  echo
}

function check_org() {
  if [[ "$DOCKERHUB_ORG" == "" ]]; then
    echo "Error: DOCKERHUB_ORG is empty"
    echo "Please ensure DOCKERHUB_ORG is set to the name of the Docker Hub organization";
    exit 1;
  else
    echo "Using Docker Hub org as $DOCKERHUB_ORG..."
  fi
}
