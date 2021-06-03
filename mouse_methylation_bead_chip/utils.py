import subprocess
from pkg_resources import resource_filename

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

def get_alignment_snakefile_fp():
    return resource_filename("mouse_methylation_bead_chip", "bead_chip_probe_alignments.smk")
