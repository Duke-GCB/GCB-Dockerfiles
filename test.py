import os
import sys
import subprocess
import yaml
import re


def get_test_list(filename):
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
  return str(subprocess.check_output(["bash", "-c", docker_cmd]))


def run_docker_get_output(imagename, cmd):
  docker_cmd = "docker run -it {} {}".format(imagename, cmd)
  try:
     return str(subprocess.check_output(["bash", "-c", docker_cmd]))
  except subprocess.CalledProcessError as exc:
     return str(exc.output)


def run_singularity_get_output(imagename, cmd):
  run_bash_cmd("docker create --name for_export {}".format(imagename))
  run_bash_cmd("sudo singularity image.create imported.img")
  run_bash_cmd("docker export for_export | sudo singularity image.import imported.img")
  run_bash_cmd("docker rm for_export")
  docker_cmd = "singularity exec imported.img {} 2>&1".format(cmd)
  try:
     print("RUN", docker_cmd)
     ret = str(subprocess.check_output(["bash", "-c", docker_cmd]))
     os.remove("imported.img")
     return ret
  except subprocess.CalledProcessError as exc:
     ret = str(exc.output)
     os.remove("imported.img")
     return ret


def run_bash_cmd(command):
  return str(subprocess.check_output(["bash", "-c", command]))


def run_tests(imagename, filename):
  had_error = False
  for cmd, expect_text in get_test_list(filename):
     expect_pattern = re.compile(expect_text, re.MULTILINE)
     print("running {} inside {}".format(cmd, imagename))
     docker_output = run_docker_get_output(imagename, cmd)
     if not re.match(expect_pattern, docker_output):
        print("Error running {}".format(cmd))
        print("Expected Regex: {}".format(expect_text))
        print("Actual: {}".format(docker_output))
        had_error = True
     sing_output = run_singularity_get_output(imagename, cmd)
     if not re.match(expect_pattern, sing_output):
        print("Error running {}".format(cmd))
        print("Expected Regex: {}".format(expect_text))
        print("Actual: {}".format(sing_output))
        had_error = True

  return had_error


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
