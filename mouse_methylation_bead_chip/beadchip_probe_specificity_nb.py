# * Setup

import pandas as pd
import jupyter_utils
import subprocess
import os

import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
import mouse_methylation_bead_chip.beadchip_probe_annotation_lib as anno_lib
import mouse_methylation_bead_chip.beadchip_probe_specificity_lib as lib
import mouse_methylation_bead_chip.utils as utils

# * Get reference genomes

# ** mm10 with alt

# from https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/seqs_for_alignment_pipelines/

proc = subprocess.run(
    f"wget -O {paths.mm10_full_fa_gz} {paths.mm10_full_fa_gz_url}",
    shell=True,
    check=True,
    executable="/bin/bash",
    capture_output=True,
    encoding="utf8",
)

# cannot run samtools faidx on the fa_gz file - the file is old, perhaps a version problem?
# there is an index available for download, but lets recreate it to be sure

proc = subprocess.run(
    f"zcat {paths.mm10_full_fa_gz} | bgzip > {paths.mm10_full_fa_re_gz}",
    shell=True,
    check=True,
    executable="/bin/bash",
    capture_output=True,
    encoding="utf8",
)

# !zcat {paths.mm10_full_fa_re_gz} | head


# *** Extract chrom 19 for testing

# first index the fa

proc = subprocess.run(
    ["samtools", "faidx", paths.mm10_full_fa_re_gz],
    capture_output=True,
    encoding="utf-8",
    check=True,
)

proc = subprocess.run(
    f"samtools faidx {paths.mm10_full_fa_re_gz} chr19 | bgzip > {paths.mm10_full_fa_re_gz_chr19_only}",
    shell=True,
    capture_output=True,
    encoding="utf-8",
    check=True,
)
# !zcat {paths.mm10_full_fa_re_gz_chr19_only} | head


# * Create indices

# ** mm10

# Test index for chrom19
proc = subprocess.run(
    ["biscuit", "index", paths.mm10_full_fa_re_gz_chr19_only],
    capture_output=True,
    encoding="utf-8",
    check=True,
)
jupyter_utils.send_desktop_notification("indexing chr19 done")
utils.run_shell_print(f"ls -l {os.path.dirname(paths.mm10_full_fa_re_gz_chr19_only)}")

proc = subprocess.run(
    ["biscuit", "index", paths.mm10_full_fa_re_gz],
    capture_output=True,
    encoding="utf-8",
    check=True,
)
jupyter_utils.send_desktop_notification("indexing mm10 done")
utils.run_shell_print(f"ls -l {os.path.dirname(paths.mm10_full_fa_re_gz_chr19_only)}")


# * Prepare probe sequences for specificity test

# paths.illumina_probes_csv
# paths.illumina_coordinate_bed
probe_sequences_fq = lib.create_probe_sequences_fa()

utils.run_shell_print(f'wc -l {probe_sequences_fq}')

# TODO
# - probe coordinates are 1-based, some dtypes may not be read correctly
# - transcribe into bed format
# - then extract all simple probe sequences
# - for first step, could just use type ii probes (which don't have Rs)
# - write to file, need probe id in read names

# * Alignment

proc3 = utils.run_shell(
    command=f"""
biscuit align \
        -t 16 \
        {paths.mm10_full_fa_re_gz} \
        {probe_sequences_fq} \
    > {paths.probe_alignment_bam}
    """,
    check = True
)
utils.run_shell_print(f'ls -l {os.path.dirname(paths.probe_alignment_bam)}')
jupyter_utils.send_desktop_notification('finished biscuit align mm10')

# print(proc3.returncode)
# print(proc3.stdout)
# print(proc3.stderr)


