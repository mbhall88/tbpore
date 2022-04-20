import subprocess
from typing import List
import shlex
import os
from loguru import logger
from pathlib import Path
import hashlib


class ExternalTool:
    def __init__(self, tool: str, input: str, output: str, params: str, logdir: Path = Path("logs")):
        self.command: List[str] = self._build_command(tool, input, output, params)
        os.makedirs(logdir, exist_ok=True)
        command_hash = hashlib.sha256(" ".join(self.command).encode("utf-8")).hexdigest()
        logfile_prefix: Path = logdir / f"{tool}_{command_hash}"
        self.out_log = f"{logfile_prefix}.out"
        self.err_log = f"{logfile_prefix}.err"

    @property
    def command_as_str(self) -> str:
        return " ".join(self.command)

    @staticmethod
    def _build_command(tool: str, input: str, output: str, params: str) -> List[str]:
        command = " ".join([tool, params, output, input])
        escaped_command = shlex.split(command)
        return escaped_command

    def run(self) -> None:
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            logger.info(f"Started running {self.command_as_str} ...")
            subprocess.check_call(self.command, stdout=stdout_fh, stderr=stderr_fh)
            logger.info(f"Done running {self.command_as_str}")
