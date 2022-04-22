import glob
import sys
from pathlib import Path
from unittest.mock import patch

from tbpore.constants import repo_root
from tbpore.external_tools import ExternalTool


class TestExternalTools:
    @patch.object(
        ExternalTool,
        ExternalTool._build_command.__name__,
        return_value=["mocked", "command", "arg"],
    )
    @patch.object(Path, Path.mkdir.__name__)
    def test___constructor(self, mkdir_mock, build_command_mock):
        logdir = Path("logs")

        external_tool = ExternalTool("tool", "input", "output", "params", logdir)

        assert external_tool.command == ["mocked", "command", "arg"]
        assert external_tool.command_as_str == "mocked command arg"
        assert (
            external_tool.out_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.out"
        )
        assert (
            external_tool.err_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.err"
        )

        build_command_mock.assert_called_once_with("tool", "input", "output", "params")
        mkdir_mock.assert_called_once_with(parents=True, exist_ok=True)

    def test___build_command___simple_command(self):
        expected_escaped_command = ["tool", "param1", "param2", "-o", "out", "-i", "in"]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i in", "-o out", "param1 param2"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___single_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i 'escaped in'", "-o 'escaped out'", "params with 'escaped arg'"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___double_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", '-i "escaped in"', '-o "escaped out"', 'params with "escaped arg"'
        )
        assert expected_escaped_command == actual_escaped_command

    def test___run(self):
        logsdir = repo_root / "tests/helpers/logs"
        logsdir.mkdir(parents=True, exist_ok=True)
        for file in logsdir.iterdir():
            file.unlink()

        external_tool = ExternalTool(
            sys.executable,
            "input",
            "output",
            str(repo_root / "tests/helpers/run_test.py"),
            logsdir,
        )

        external_tool.run()

        out_file = glob.glob(f"{logsdir}/*.out")[0]
        with open(out_file) as out_file_fh:
            lines = out_file_fh.readlines()
            assert lines == ["out\n"]

        err_file = glob.glob(f"{logsdir}/*.err")[0]
        with open(err_file) as err_file_fh:
            lines = err_file_fh.readlines()
            assert lines == ["err\n"]
