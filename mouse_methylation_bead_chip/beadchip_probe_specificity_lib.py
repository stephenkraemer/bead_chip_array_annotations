import pandas as pd
import jupyter_utils
import subprocess
import os
import re
from pathlib import Path
from io import StringIO

import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
import mouse_methylation_bead_chip.utils as utils

import numpy as np



def scratch():
    probes_csv = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/epic-arrays-for-hematopoiesis/Infinium_20Mouse_20Methylation_20v1.0_20A1_20GS_20Manifest_20File.csv"
    # !head {probes_csv}
    probes_df = pd.read_csv(probes_csv, skiprows=7)
    # dtypewarning
    # comment on new illumin aid on github annotation page
    probes_df.iloc[10].to_frame()
    # AlleleA_ProbeSeq
    # Infinium_Design_Type
    probes_df.query('Infinium_Design_Type == "II"')["AlleleA_ProbeSeq"].head()
    # definitely contains Rs

    # no b seq for type ii, check
    assert (
        probes_df.query('Infinium_Design_Type == "II"')["AlleleB_ProbeSeq"].isna().all()
    )

    # both probes defined
    probes_df.query('Infinium_Design_Type == "I"')["AlleleA_ProbeSeq"].head()
    probes_df.query('Infinium_Design_Type == "I"')["AlleleB_ProbeSeq"].head()


def create_probe_sequences_fa():

    # exchange with curated and make arg
    # dtypewarning
    probes_df = pd.read_csv(paths.illumina_probes_csv, skiprows=7)

    # is illumina id unique?
    assert not probes_df.IlmnID.duplicated().any()

    # simple first: only type ii, a + b
    df1 = pd.melt(
        probes_df.query('Infinium_Design_Type == "I"')[
            ["IlmnID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq"]
        ],
        id_vars=["IlmnID"],
        var_name="allele",
        value_name="sequence",
    )
    # highest possible Q score here: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm is "I" for phred score 40
    # found no indication that biscuit is intended to work with fasta files directly (but it can take single reads)

    # pd.Series(['I', 'I', 'I']).str.repeat(pd.Series([1, 2, 3]))
    df2 = df1.assign(
        read_name=lambda df: "@" + df["IlmnID"] + "_" + df["allele"],
        plus="+",
        qual=lambda df: pd.Series(["I"] * df.shape[0]).str.repeat(
            df["sequence"].str.len()
        ),
    )

    probes_fq_content = (
        df2[["read_name", "sequence", "plus", "qual"]]
        .apply(lambda ser: ser.str.cat(sep="\n"), axis=1)
        .str.cat(sep="\n")
    )

    Path(paths.probe_sequences_fq).write_text(probes_fq_content)

    utils.run_shell_print(f"head {paths.probe_sequences_fq}")

    return paths.probe_sequences_fq


def bam_to_probe_df():
    bam_fp = paths.probe_alignment_bam

    # txt = utils.run_shell(f"samtools view {bam_fp} | head -n 2").stdout
    # re.sub(r'\t([^:]+?:[^:]+?:[^:]+?)', r'||\1', txt)
    # how to get tags: re.sub to remove first always tehre fields, get series of remainder, then extract

    df = pd.read_csv(
        StringIO(utils.run_shell(f"samtools view {bam_fp}").stdout),
        sep="\t",
        names=[
            "qname",
            "flag",
            "rname",
            "pos",
            "mapq",
            "cigar",
            "rnext",
            "pnext",
            "tlen",
            "seq",
            "qual",
        ],
        header=None,
        usecols=np.arange(11),
    )[["qname", "flag", "rname", "pos", "mapq", "cigar", "seq"]]


    df.mapq.plot.hist()
