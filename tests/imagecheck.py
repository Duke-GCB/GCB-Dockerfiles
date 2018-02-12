#!/usr/bin/env python3
# Recursively looks for files in current directory.
# For each unittest.yml found
#   1) builds a docker image in the directory containing unittest.yml
#   2) runs tests specified in unittest.yml inside the docker image using docker
#   3) runs tests specified in unittest.yml inside the docker image using singularity

import os
import sys
import subprocess
import yaml
import re

UNITTEST_FILENAME = "unittest.yml"
TEST_WORKDIR = "/tmp"
TEST_USER = "nobody"


def run_bash_cmd(command, ignore_non_zero_exit_status=False):
    """
    Run bash command and return output
    :param command: str bash command to run
    :param ignore_non_zero_exit_status: bool: when true returns output no matter the exit status
    :return: str: output from command
    """
    try:
        return str(subprocess.check_output(["bash", "-c", command]))
    except subprocess.CalledProcessError as exc:
        if ignore_non_zero_exit_status:
            return str(exc.output)
        else:
            print("ERROR", exc.returncode, exc.output, exc.stderr)
            raise


def get_test_list(filename):
    """
    Read YAML test configuration from a file
    :param filename: str: path to YAML file with test settings for a single docker file
    :return: [(command, expected_text)]: array of command and expected output pairs
    """
    with open(filename) as infile:
        testinfo_list = []
        unittest_config = yaml.load(infile)
        for testinfo in unittest_config['commands']:
            cmd = testinfo['cmd']
            expect_text = testinfo['expect_text']
            testinfo_list.append((cmd, expect_text))
        return testinfo_list


def build_docker_image(imagename, build_path):
    print("Building docker image for path {}".format(build_path))
    run_bash_cmd("docker build -t {} {}".format(imagename, build_path))


def run_docker_get_output(imagename, cmd, workdir=None, user=None):
    options = ""
    if workdir:
        options += "--workdir {} ".format(workdir)
    print("Testing image {} in docker with cmd {} {}".format(imagename, cmd, options))
    if user:
        options += "--user {} ".format(user)
    options += "-it "
    docker_cmd = "docker run {}{} {}".format(options, imagename, cmd)
    return run_bash_cmd(docker_cmd, ignore_non_zero_exit_status=True)


def run_tests(imagename, filename):
    had_error = False
    for cmd, expect_text in get_test_list(filename):
        expect_pattern = re.compile(expect_text, re.DOTALL)
        docker_output = run_docker_get_output(imagename, cmd)
        if not re.match(expect_pattern, docker_output):
            print_test_error(cmd, expect_text, docker_output)
            had_error = True
        docker_output_with_options = run_docker_get_output(imagename, cmd, workdir=TEST_WORKDIR, user=TEST_USER)
        if not re.match(expect_pattern, docker_output_with_options):
            print_test_error(cmd + " (with workdir and user options)", expect_text, docker_output_with_options)
            had_error = True
    return had_error


def print_test_error(cmd, expect_text, cmd_output):
    print("Error running {}".format(cmd))
    print("Expected Regex: {}".format(expect_text))
    print("Actual: {}".format(cmd_output))


def find_unittest_info():
    """
    Recursively walk the current directory looking for tests.
    :return: [(test_dir, imagename, test_filename)]: list of testing data tuples
    """
    test_info = []
    for root, dirs, files in os.walk("."):
        for name in files:
            if name == UNITTEST_FILENAME:
                test_directory = root
                imagename = "localhost:5000/test_{}".format(test_directory.replace("./", "").replace("/", ":").replace("+","_"))
                test_filename = os.path.join(root, name)
                test_info.append((test_directory, imagename, test_filename))
    return test_info


def get_test_file_path(file_path):
    filename = os.path.basename(file_path)
    if filename in ["Dockerfile", UNITTEST_FILENAME]:
        parent_directory = os.path.dirname(file_path)
        return "{}/{}".format(parent_directory, UNITTEST_FILENAME)
    return None


def get_unittest_file_paths(path_list):
    unittest_paths = set()
    for file_path in path_list:
        unittest_path = get_test_file_path(file_path)
        if unittest_path:
            unittest_paths.add(unittest_path)
    return unittest_paths


def find_and_run_tests(owner, changed_paths):
    had_errors = False
    tested_images = 0
    images_with_errors = 0
    for unittest_path in get_unittest_file_paths(changed_paths):
        parts = unittest_path.split(sep="/")
        if len(parts):
            tool, tag, _ = parts
            imagename = "{}/{}:{}".format(owner, tool, tag).replace("+","_")
            build_docker_image(imagename, "{}/{}".format(tool, tag))
            had_error = run_tests(imagename, unittest_path)
            tested_images += 1
            if had_error:
                images_with_errors += 1
                had_errors = True
        else:
            print("Skipping {}".format(unittest_path))
    print("Tested {} images. Images with errors: {}".format(tested_images, images_with_errors))
    return had_errors


def main():
    if len(sys.argv) < 3:
        print("Usage python3 tests/imagecheck.py <docker_owner> <unittest_or_dockerfile_path>...")
        sys.exit(1)
    else:
        owner = sys.argv[1]
        changed_paths = sys.argv[2:]
        had_errors = find_and_run_tests(owner, changed_paths)
        if had_errors:
            sys.exit(2)


if __name__ == "__main__":
    main()
