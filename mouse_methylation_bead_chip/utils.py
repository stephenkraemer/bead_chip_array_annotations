import subprocess

def run_shell(command, check=True):
    return subprocess.run(
        command,
        shell=True,
        capture_output=True,
        encoding="utf-8",
        check=check,
    )

def run_shell_print(command):
    print(subprocess.run(
        command,
        shell=True,
        capture_output=True,
        encoding="utf-8",
        check=True,
    ).stdout)
