import os
import sys
import subprocess
import yaml
import re


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
    docker_cmd = "docker build -t {} {}".format(imagename, build_path)
    return run_bash_cmd(docker_cmd)


def run_docker_get_output(imagename, cmd):
    docker_cmd = "docker run -it {} {}".format(imagename, cmd)
    return run_bash_cmd(docker_cmd, ignore_non_zero_exit_status=True)


def run_singularity_get_output(imagename, cmd):
    try:
        run_bash_cmd("docker create --name for_export {}".format(imagename))
        try:
            run_bash_cmd("sudo singularity image.create imported.img")
            run_bash_cmd("docker export for_export | sudo singularity image.import imported.img")
            singularity_cmd = "singularity exec imported.img {} 2>&1".format(cmd)
            singularity_output = run_bash_cmd(singularity_cmd, ignore_non_zero_exit_status=True)
        finally:
            os.remove("imported.img")
    finally:
        run_bash_cmd("docker rm for_export")
    return singularity_output


def run_tests(imagename, filename):
    had_error = False
    for cmd, expect_text in get_test_list(filename):
        expect_pattern = re.compile(expect_text, re.MULTILINE)
        print("running {} inside {}".format(cmd, imagename))
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


def find_unittest_info(dir_name):
    test_info = []
    for root, dirs, files in os.walk("."):
        for name in files:
            if name == "unittest.yml":
                test_directory = root
                image_name = "test_{}".format(test_directory.replace("./", "").replace("/", ":"))
                test_filename = os.path.join(root, name)
                test_info.append((test_directory, image_name, test_filename))
    return test_info


if __name__ == "__main__":
    base_directory = "."
    for test_directory, image_name, test_filename in find_unittest_info(base_directory):
        print(test_directory, image_name, test_filename)

    build_docker_image(image_name, test_directory)
    had_error = run_tests(image_name, test_filename)
    if had_error:
        sys.exit(1)
