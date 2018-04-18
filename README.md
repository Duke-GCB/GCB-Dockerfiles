GCB-Dockerfiles
==========

This repository houses `Dockerfile`s and tests for producing versioned container images of computational genomics tools.

Images produced are are pushed to [Docker Hub](https://hub.docker.com/) under the [dukegcb organization](https://hub.docker.com/u/dukegcb/) and suitable for use in [CWL](https://commonwl.org) workflows or under [Singularity](http://singularity.lbl.gov).

The directory structure of this repository uses a convention to determine the names and tags of built docker images, which is explained below.

## Directory Structure

As this repository houses many `Dockerfile`s, the directory structure allows those to be organized by tool name and version, allowing automatic builds and versioning conventions.

`Dockerfile`s are organized into subdirectories based on the tool and version. For example, the Dockerfile for FastQC version 0.11.4 is located at `fastqc/0.11.4/Dockerfile`.

### Versioning and Tags

To preserve versions and support the Docker `latest` tag, we adopt two conventions when building:

1. The name of the directory immediately containing the `Dockerfile` is used as the Docker image tag when building
2. Creating a `latest` symlink to point to one of the version directories will cause that version to also be tagged as latest

For example, the following directory structure:

```
$ ls -ld fastqc/*
drwxr-xr-x  4 dcl9  staff  128 Feb 14 13:57 fastqc/0.11.4
lrwxr-xr-x  1 dcl9  staff    6 Feb 14 10:52 fastqc/latest -> 0.11.4
```

Will build a Docker image `dukegcb/fastqc:0.11.4` from the `Dockerfile` in `fastqc/0.11.4`, and also tag that image as `dukegcb/fastqc:latest`.

### Naming Restrictions

When naming directories, it's important to use names that are also valid for Docker images and tags. Sometimes, tool names and versions may contain characters (e.g. `+`) that are not valid for Docker image names or tags.

From [docs.docker.com](https://docs.docker.com/engine/reference/commandline/tag/#extended-description), regarding image names:

> Name components may contain **lowercase** letters, digits and separators. A separator is defined as a period, one or two underscores, or one or more dashes. A name component may not start or end with a separator.

Tags:

> A tag name must be valid ASCII and may contain lowercase and uppercase letters, digits, underscores, periods and dashes. A tag name may not start with a period or a dash and may contain a maximum of 128 characters.

## Continuous Integration

This repository includes scripts to automate the building, testing, and pushing of Docker images to the public Docker hub registry.

The scripts that run these processes are [build-ci.sh](build-ci.sh), [test-ci.sh](test-ci.sh), and [deploy-ci.sh](deploy-ci.sh). Each depends on the conventions denoted above, and use functions from [functions.sh](functions.sh) to determine what to process at each stage.

These processes happen on [CircleCI](https://circleci.com/gh/Duke-GCB/GCB-Dockerfiles), as configured in  [.circle/config.yml](.circle/config.yml), and expect the following 4 environment variables:

- `DEPLOY_BRANCH` - When changes are made on this git branch (typically `master`), built images will be pushed to Docker Hub
- `DOCKERHUB_ORG` - The Docker organization name to use when naming and pushing built images (e.g. `dukegcb`)
- `DOCKERHUB_USERNAME` - The Docker Hub username to use when authenticating with Docker Hub. See below.
- `DOCKERHUB_PASSWORD` - The Docker Hub password to use when authenticating with Docker Hub. See below.

### Docker Hub Credentials

For the CI service to push to Docker Hub, it must login with `docker login`. It uses the username/password set in the above variables.

To avoid using personal account credentials, we use a Docker Hub account created just for this purpose. This account is a member of the **cibuild** team, which has write access to the docker repos.

### Testing

Beside each `Dockerfile` there should be a `unittest.yml` file describing how to test the built image.  These tests are intended to confirm that the expected tool has been installed in the image and executes when run.

The [unittest.yml](fastqc/0.11.4/unittest.yml) for FastQC simply asserts that running the `fastqc -h` command produces usage text:

```
commands:
  - cmd: "fastqc -h"
    expect_text: ".*fastqc seqfile1 seqfile2.*"
```

### Build

Dockerfiles are built with the command noted in [functions.sh](functions.sh):

```
docker build -f $tool/$version/Dockerfile" -t "$owner/$tool:$version" "$tool/$version"
```

So this build command can easily be run locally to confirm the image builds correctly.

## Contributing

### Adding a new tool

1. Create a repo under the dukegcb Docker Hub organization with the tool name
2. Under collaborators, add the **cibuild** team with **write** access
3. Create a branch for the new tool in your local copy of the repo
4. Follow steps 2-10 under **Adding a new tool version**

### Adding a new tool version

1. Create a branch for new version in your local copy of the repo.
2. Create a directory for the version in the tool's directory.
3. Create and edit the `Dockerfile` in that directory to install that version of the tool.
4. Confirm that the `Dockerfile` builds as expected.
5. Write a `unittest.yml` file to test that the tool runs as expected.
6. If you wish to update the `latest` tag to this version, update the `latest` symlink to point to the new version `rm latest && ln -s new-version latest`.
7. Push your branch to GitHub
8. CircleCI will build and test changed `Dockerfile`s from your branch, but images will not be pushed to Docker Hub.
9. Submit a pull request to merge those changes into `master`
10. After review and merge to `master`, the newly built image(s) will be pushed to Docker Hub.
