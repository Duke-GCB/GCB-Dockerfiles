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


def start_local_docker_registry():
    run_bash_cmd("docker run -d -p 5000:5000 --restart=always --name registry registry:2")


def stop_local_docker_registry():
    run_bash_cmd("docker rm --force registry")


def build_and_push_docker_image(imagename, build_path):
    print("Building docker image for path {}".format(build_path))
    run_bash_cmd("docker build -t {} {}".format(imagename, build_path))
    print("Pushing docker image {}".format(imagename))
    run_bash_cmd("docker push {}".format(imagename))


def run_docker_get_output(imagename, cmd):
    print("Testing image {} in docker with cmd {}".format(imagename, cmd))
    docker_cmd = "docker run -it {} {}".format(imagename, cmd)
    return run_bash_cmd(docker_cmd, ignore_non_zero_exit_status=True)


def run_singularity_get_output(imagename, cmd):
    print("Testing image {} in singularity with cmd {}".format(imagename, cmd))
    singularity_imagename = "docker://{}".format(imagename)
    singularity_cmd = "SINGULARITY_NOHTTPS=yes singularity exec {} {} 2>&1".format(singularity_imagename, cmd)
    singularity_output = run_bash_cmd(singularity_cmd, ignore_non_zero_exit_status=True)
    return singularity_output


def run_tests(imagename, filename):
    had_error = False
    for cmd, expect_text in get_test_list(filename):
        expect_pattern = re.compile(expect_text, re.MULTILINE)
        docker_output = run_docker_get_output(imagename, cmd)
        if not re.match(expect_pattern, docker_output):
            print_test_error(cmd, expect_text, docker_output)
            had_error = True
        sing_output = run_singularity_get_output(imagename, cmd)
        if not re.match(expect_pattern, sing_output):
            print_test_error(cmd, expect_text, sing_output)
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


def main():
    start_local_docker_registry()
    try:
        had_error = False
        for test_directory, imagename, test_filename in find_unittest_info():
            build_and_push_docker_image(imagename, test_directory)
            had_error = run_tests(imagename, test_filename)
        if had_error:
            sys.exit(1)
    finally:
        stop_local_docker_registry()


if __name__ == "__main__":
    main()
