import unittest
from unittest.mock import patch, mock_open, call
import subprocess
import imagecheck


class TestRunFunctions(unittest.TestCase):
    @patch('imagecheck.subprocess')
    def test_run_bash_cmd(self, mock_subprocess):
        mock_subprocess.check_output.return_value = "<file listing>"
        resp = imagecheck.run_bash_cmd("ls")
        self.assertEqual(resp, "<file listing>")
        mock_subprocess.check_output.assert_called_with(['bash', '-c', 'ls'])

    @patch('imagecheck.subprocess')
    def test_run_bash_cmd_raise_error(self, mock_subprocess):
        mock_subprocess.CalledProcessError = subprocess.CalledProcessError
        mock_subprocess.check_output.side_effect = subprocess.CalledProcessError(returncode=1, cmd="foo",
                                                                                 output="someoutput")

        # by default or requested raise on non-zero exit status
        with self.assertRaises(subprocess.CalledProcessError):
            imagecheck.run_bash_cmd("ls")
        with self.assertRaises(subprocess.CalledProcessError):
            imagecheck.run_bash_cmd("ls", ignore_non_zero_exit_status=False)

        # ignore error if explicitly requested
        output = imagecheck.run_bash_cmd("ls", ignore_non_zero_exit_status=True)
        self.assertEqual(output, "someoutput")

    def test_get_test_list(self):
        unittest_yaml = """
        commands:
          - cmd: "bwa"
            expect_text: ".*Usage:.*bwa.*"
          - cmd: "samtools"
            expect_text: ".*Usage:.*samtools.*"
        """
        with patch("builtins.open", mock_open(read_data=unittest_yaml)) as mock_file:
            testinfo_list = imagecheck.get_test_list("somefile.txt")
            self.assertEqual(len(testinfo_list), 2)
            self.assertEqual(testinfo_list[0], ("bwa", ".*Usage:.*bwa.*"))
            self.assertEqual(testinfo_list[1], ("samtools", ".*Usage:.*samtools.*"))

    @patch('imagecheck.run_bash_cmd')
    def test_run_docker_get_output(self, mock_run_bash_cmd):
        mock_run_bash_cmd.return_value = "/tmp/somedir"
        result = imagecheck.run_docker_get_output("someimage", "pwd")
        self.assertEqual(result, "/tmp/somedir")
        mock_run_bash_cmd.assert_has_calls([
            call('docker run -it --rm someimage pwd', ignore_non_zero_exit_status=True),
        ])

    @patch('imagecheck.run_bash_cmd')
    def test_run_docker_get_output_with_options(self, mock_run_bash_cmd):
        mock_run_bash_cmd.return_value = "/tmp/somedir"
        result = imagecheck.run_docker_get_output("someimage", "pwd", "/work", "ubuntu")
        self.assertEqual(result, "/tmp/somedir")
        mock_run_bash_cmd.assert_has_calls([
            call('docker run --workdir /work --user ubuntu -it --rm someimage pwd', ignore_non_zero_exit_status=True),
        ])

    @patch('imagecheck.get_test_list')
    @patch('imagecheck.run_docker_get_output')
    @patch('imagecheck.print_test_error')
    def test_run_tests_no_error(self, mock_print_test_error,
                                mock_run_docker_get_output, mock_get_test_list):
        mock_get_test_list.return_value = [
            ("bwa", ".*Usage: bwa.*")
        ]
        mock_run_docker_get_output.return_value = "Usage: bwa"
        had_error = imagecheck.run_tests("someimage", "/tmp/fakefile.yml")
        mock_print_test_error.assert_not_called()
        self.assertEqual(had_error, False)

    @patch('imagecheck.get_test_list')
    @patch('imagecheck.run_docker_get_output')
    @patch('imagecheck.print_test_error')
    def test_run_tests_docker_error(self, mock_print_test_error, mock_run_docker_get_output, mock_get_test_list):
        mock_get_test_list.return_value = [
            ("bwa", ".*Usage: bwa.*")
        ]
        mock_run_docker_get_output.return_value = "Error: no such file"
        had_error = imagecheck.run_tests("someimage", "/tmp/fakefile.yml")
        mock_print_test_error.assert_has_calls([
            call('bwa', '.*Usage: bwa.*', 'Error: no such file')
        ])
        self.assertEqual(had_error, True)

@patch('imagecheck.sys')
@patch('imagecheck.find_and_run_tests')
class TestMain(unittest.TestCase):

    def test_exits_1_without_enough_args(self, mock_run, mock_sys):
        mock_sys.argv = ['imagecheck.py']
        imagecheck.main()
        mock_sys.exit.assert_has_calls([call(1),])
        self.assertFalse(mock_run.called)

    def test_exits_2_if_failure(self, mock_run, mock_sys):
        mock_sys.argv = ['imagecheck.py', 'dockerorg']
        mock_run.return_value = True
        imagecheck.main()
        mock_sys.exit.assert_has_calls([call(2),])
        self.assertTrue(mock_run.called)

    def test_calls_find_and_run_without_changed_files(self, mock_run, mock_sys):
        mock_sys.argv = ['imagecheck.py', 'dockerorg']
        mock_run.return_value = False
        imagecheck.main()
        self.assertFalse(mock_sys.exit.called)
        self.assertTrue(mock_run.called)
        mock_run.assert_has_calls([
            call('dockerorg', [])
        ])

    def test_calls_find_and_run_with_changed_files(self, mock_run, mock_sys):
        mock_sys.argv = ['imagecheck.py', 'dockerorg', 'path1', 'path2']
        mock_run.return_value = False
        imagecheck.main()
        self.assertFalse(mock_sys.exit.called)
        self.assertTrue(mock_run.called)
        mock_run.assert_has_calls([
            call('dockerorg', ['path1','path2'])
        ])
