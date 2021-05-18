# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.0
#   kernelspec:
#     display_name: Python [conda env:mouse_hema_meth_dev3] *
#     language: python
#     name: conda-env-mouse_hema_meth_dev3-xpython
# ---

# + [markdown] tags=[] heading_collapsed="true"
# # TODO
# -

# - check sorting order of chromosomes throughout analysis
# - there is a bug in gtfanno which makes it fail to read the appris principal score, does this influence the results?
# - there is a "transcript = 0" count in the gtfanno basic stat output. is this indicative of a problem?
# - double check the changes made to allow DCRD intervals enveloping the promoter interval
# - check annos in IGV
# - are feature coordinates really 0-based, right open?

# annos to add
# - motif: CG, CHH, CHG
# - strand, illumina strands are not cytosine strands
# - hematopoietic regions
#     - cis reg atlas
#     - vision
#     - amit enhancers
# - general regulatory regions
#     - ensembl reg regions
#     - chrom hmm - ask maxi again what he had in mind here - forgot which resource he mentioned
# - tfbs

# + [markdown] tags=[]
# # Setup

# + [markdown] tags=[] heading_collapsed="true"
# ## Resource parameters
# -

n_cores = 12

# + [markdown] tags=[] heading_collapsed="true"
# ## Imports

# + tags=[]
# isort: off
import os

num_threads = str(n_cores)

# these need to be set prior to numpy import
os.environ["OMP_NUM_THREADS"] = num_threads
os.environ["OPENBLAS_NUM_THREADS"] = num_threads
os.environ["MKL_NUM_THREADS"] = num_threads
os.environ["VECLIB_MAXIMUM_THREADS"] = num_threads
os.environ["NUMEXPR_NUM_THREADS"] = num_threads

import numpy as np

# isort: on

import subprocess
import tempfile

import gtfanno as ga
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from IPython.display import display

import mouse_hema_meth.utils as ut
# -

# %matplotlib inline

import mouse_hema_meth.methylome.annotation.epic_array_probe_annotation_lib as lib

# + [markdown] tags=[] heading_collapsed="true"
# ## Rerun flags
# -

recompute = True

# + [markdown] tags=[] heading_collapsed="true"
# ## Dtypes
# -

chrom_dtype_prefixed = pd.api.types.CategoricalDtype(
    categories=[
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chrX",
        "chrY",
        "chrMT",
    ],
    ordered=True,
)

# + [markdown] tags=[]
# # Paths

# + [markdown] tags=[] heading_collapsed="true"
# ## Project paths
# -

project_dir = (
    "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/epic-arrays-for-hematopoiesis"
)

temp_dir_obj = tempfile.TemporaryDirectory(dir=project_dir)
temp_dir_name = temp_dir_obj.name
temp_dir_name

# + [markdown] tags=[] heading_collapsed="true"
# ## Gencode
# -

gencode_download_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

# unmodified gencode gtf (note: with chr prefix)

gencode_gtf = project_dir + "/gencode.vM25.annotation.gtf.gz"
gencode_gtf

# gencode filtered for principal transcripts of protein coding genes, note that that chromosome prefix ('chr') is removed in this file

gencode_coding_canonical_gtf = (
    project_dir + "/gencode.vM25.annotation_coding_canonical.gtf.gz"
)
gencode_coding_canonical_gtf

# + [markdown] tags=[] heading_collapsed="true"
# ## Probes
# -

# path to original probe file

original_probes_bed = project_dir + "/2021-03-22_mmbc_probes.bed"

# !head {original_probes_bed}

# path to reformatted probe file

reformatted_probes_bed = project_dir + "/2021-03-22_mmbc_probes_reformatted.bed"

# full Illumina probe file

illumina_probes_url = "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/Infinium%20Mouse%20Methylation%20v1.0%20A1%20GS%20Manifest%20File.csv"
illumina_probes_csv = (
    project_dir
    + "/Infinium_20Mouse_20Methylation_20v1.0_20A1_20GS_20Manifest_20File.csv"
)
illumina_probes_csv

# full Illumina probe file BED coordinates

illumina_coordinate_bed = project_dir + "/illumina-all-probes.bed"
illumina_coordinate_bed

# + [markdown] tags=[]
# ## Gene annotation results

# + [markdown] heading_collapsed="true" tags=[]
# ### gtfanno results
# -

custom_intervals_results_dir = project_dir + "/custom-intervals_1500-500"
os.makedirs(custom_intervals_results_dir, exist_ok=True)

custom_intervals_trunk_path = (
    custom_intervals_results_dir + "/custom-intervals_1500-500"
)
custom_intervals_trunk_path

custom_intervals_results_paths_d = dict(
    primary_annos_bed=custom_intervals_trunk_path + "_primary-annotations.bed",
    primary_annos_p=custom_intervals_trunk_path + "_primary-annotations.p",
    all_annos_bed=custom_intervals_trunk_path + "_all-annotations.bed",
    all_annos_p=custom_intervals_trunk_path + "_all-annotations.p",
)
custom_intervals_results_paths_d

# + [markdown] tags=[] heading_collapsed="true"
# ### Final tables
# -

gene_annos_primary_one_row = project_dir + "/gene-annos_primary_one-row.bed"
print(gene_annos_primary_one_row)
gene_annos_primary_multi_row = project_dir + "/gene-annos_primary_multi-row.bed"
print(gene_annos_primary_multi_row)

# + [markdown] tags=[]
# # Analysis

# + [markdown] tags=[] heading_collapsed="true"
# ## Prepare input data

# + [markdown] tags=[] heading_collapsed="true"
# ### CpG island annos

# +
import mouse_hema_meth.genome_annotations.get_genome_annos_paths as get_genome_annos_paths

cpg_islands_pickle_d = get_genome_annos_paths.cpg_islands_shores_shelves_pickle_paths_d

# + [markdown] tags=[] heading_collapsed="true"
# ### Prepare gene annotation

# + [markdown] tags=[] heading_collapsed="true"
# #### download gencode

# + tags=[]
if recompute:
    subprocess.run(
        ["wget", "-O", gencode_gtf, gencode_download_url],
        check=True,
    )

# + tags=[]
# !zcat {gencode_gtf} | head -n 6

# + [markdown] tags=[] heading_collapsed="true"
# #### Filter and reformat gencode GTF
# -

# - restrict to canonical transcripts
# - restrict to coding transcripts
# - remove chr prefix
# - change M to MT

gencode_df = pr.read_gtf(gencode_gtf, as_df=True, duplicate_attr=True)

# extract appris principal score from tags
appris_principal_score = (
    gencode_df["tag"].str.extract(r"appris_principal_(\d)", expand=False).astype(float)
)

appris_principal_score.value_counts()

appris_principal_score.isnull().sum()

appris_principal_score.notnull().sum()

is_principal_transcript = appris_principal_score.notnull()

is_protein_coding = gencode_df["gene_type"].eq("protein_coding")

gencode_df_coding_canonical = gencode_df.loc[
    is_principal_transcript & is_protein_coding
].copy()

gencode_df_coding_canonical.head(3)

gencode_df_coding_canonical.shape

gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].str.replace("chr", "")
gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].replace("M", "MT")

gencode_pr = pr.PyRanges(gencode_df_coding_canonical)
gencode_pr.df.Chromosome.unique()

gencode_pr.to_gtf(gencode_coding_canonical_gtf)

# !zcat {gencode_coding_canonical_gtf} | head

# verify gtf

# + tags=[] jupyter={"outputs_hidden": true}
# !zcat {gencode_coding_canonical_gtf} | grep ^protein_coding

# + tags=[] jupyter={"outputs_hidden": true}
# !zcat {gencode_coding_canonical_gtf} | grep ^appris

# + [markdown] tags=[]
# ### Prepare and inspect probes files

# + [markdown] tags=[] heading_collapsed="true"
# #### Probe file from Maxi

# + [markdown] tags=[] heading_collapsed="true"
# ##### Inspect original probes file
# -

# - file has duplicates
# - file is not fully sorted

# + [markdown] tags=[]
# ###### General overview
# -

# !head -n 3 {original_probes_bed}

# !cut -f 1 < {original_probes_bed} | uniq

original_probes_df = pd.read_csv(
    original_probes_bed,
    sep="\t",
    header=None,
    names=["Chromosome", "Start", "End", "name"],
)
original_probes_df["Chromosome"] = pd.Categorical(
    original_probes_df["Chromosome"],
    categories=original_probes_df["Chromosome"].unique(),
    ordered=True,
)
original_probes_df

# + [markdown] tags=[]
# ###### File is not fully sorted
# -

# **Note that the original probes df is not completely sorted on Start/End**

original_probes_df_sorted = original_probes_df.sort_values(
    ["Chromosome", "Start", "End"]
).reset_index(drop=True)
original_probes_df_sorted

original_probes_df_sorted.Chromosome.dtype

# + [markdown] tags=[] heading_collapsed="true"
# ###### Several probes are present with the same coordinates, but different names
# -

original_probes_df.loc[
    original_probes_df.duplicated(["Chromosome", "Start", "End"], keep=False)
]

original_probes_df.loc[
    original_probes_df.duplicated(["Chromosome", "Start", "End", "name"], keep=False)
]

# + [markdown] tags=[] heading_collapsed="true"
# ##### Reformat probes file
# -

# - need to resort
# - need to remove chr prefix
# - drop duplicates

probes_df_no_prefix_sorted = (
    original_probes_df.assign(
        Chromosome=lambda df: df["Chromosome"].str.replace("chr", ""),
    )[["Chromosome", "Start", "End"]]
    .drop_duplicates()
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)

probes_df_no_prefix_sorted.to_csv(
    reformatted_probes_bed, sep="\t", header=False, index=False
)

# !head {reformatted_probes_bed}

# + [markdown] tags=[] heading_collapsed="true"
# #### Illumina probe file

# + [markdown] tags=[]
# ##### Schema
# -

# - MFG_CHANGE probes haben ein problem
# - there may be one row separating assay probes from controls somewhere in the dataframe? (info from Maxi)

# + [markdown] tags=[] heading_collapsed="true"
# ##### Download
# -

if recompute:
    subprocess.run(["wget", "-O", illumina_probes_csv, illumina_probes_url], check=True)

# !head {illumina_probes_csv}

# + [markdown] tags=[] heading_collapsed="true"
# ##### Get curated BED intervals for probes
# -

illumina_probes = pd.read_csv(
    illumina_probes_csv,
    skiprows=7,
    dtype={
        "AddressA_ID": str,
        "CHR": str,
        "MFG_Change_Flagged": "boolean",
        "MAPINFO": "Int64",
    },
)

# Fields, drop fields with longish sequence strings for display

# + tags=[]
illumina_probes.drop(["Forward_Sequence", "Top_Sequence"], axis=1).iloc[0].to_frame()
# -

# There are nan chromosomes entries, and also some entries for chromosome 0, just 410, so I assume this can just be discarded as controls or something similar

illumina_probes.CHR.value_counts()

# checked manually in my index files: 1-based Start info is in MAPINFO

# - for comparability with Maxis probes, also add 'chr' prefix and make Categorical
# - provide BED interval for cytosine

illumina_probes_curated_chrom_defined = (
    illumina_probes[["CHR", "MAPINFO", "IlmnID"]]
    .rename(columns={"CHR": "Chromosome", "MAPINFO": "Start", "IlmnID": "name"})
    .loc[lambda df: df.Chromosome.notnull() & df.Chromosome.ne("0")]
    .assign(
        Start=lambda df: df["Start"] - 1,
        End=lambda df: df["Start"] + 1,
        Chromosome=lambda df: ("chr" + df["Chromosome"]).astype(chrom_dtype_prefixed),
    )
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)[["Chromosome", "Start", "End", "name"]]
)
illumina_probes_curated_chrom_defined

# drop duplicate rows, remove prefix, change to alphabetic sorting order

illumina_probes_curated_chrom_defined.assign(
    Chromosome=lambda df: df.Chromosome.astype(str).str.replace("chr", "")
).iloc[:, 0:3].sort_values(["Chromosome", "Start", "End"]).drop_duplicates().to_csv(
    illumina_coordinate_bed, sep="\t", header=False, index=False
)

# !head {illumina_coordinate_bed}

# + [markdown] tags=[] heading_collapsed="true"
# ##### Check against Maxis probes to see whether I have correct manifest file
# -

# this is the correct manifest file - maxis coordinates are shifted when on minus strand

pd.merge(
    original_probes_df_sorted,
    illumina_probes_curated_chrom_defined,
    on=["Chromosome", "Start", "End", "name"],
    how="inner",
)

df = pd.merge(
    original_probes_df_sorted,
    illumina_probes_curated_chrom_defined,
    on=["name"],
    how="inner",
)
display(df)
assert df.shape[0] == original_probes_df_sorted.shape[0]

# + [markdown] tags=[] heading_collapsed="true"
# ##### Add motif and strand

# + [markdown] tags=[]
# ## Annotation

# + [markdown] tags=[]
# ### Gene annotation

# + [markdown] tags=[] heading_collapsed="true"
# #### Perform annotation

# + tags=[]
# %%time
ga.annotate(
    query_bed=illumina_coordinate_bed,
    gtf_fp=gencode_coding_canonical_gtf,
    trunk_path=custom_intervals_trunk_path,
    tmpdir=temp_dir_name,
    promoter=(-1500, 500),
    distant_cis_regulatory_domain=(-100_000, 100_000),
)

# + [markdown] tags=[]
# #### Inspect annotations
# -

primary_annos = pd.read_pickle(custom_intervals_results_paths_d["primary_annos_p"])

primary_annos.shape

# + [markdown] tags=[]
# ##### General checks
# -

primary_annos.query('feat_class == "Promoter"').head(3)

primary_annos.query('feat_class == "exon"').head(3)

# + [markdown] tags=[]
# ##### Multiple assignments per region

# + [markdown] tags=[]
# ###### How is this distributed across feature classes?
# -

multi_annos_crosstab = (
    primary_annos.groupby(["feat_class", "gtfanno_uid"], observed=True)
    .size()
    .groupby("feat_class")
    .value_counts()
    .unstack()
)
multi_annos_crosstab

# + active=""
# multi_annos_crosstab.to_clipboard()

# + [markdown] tags=[]
# ###### Example for Promoter multiple annotations - random samples indicate that these are indeed ambiguous sites
# -

primary_annos["is_duplicated"] = primary_annos.duplicated(
    subset=["Chromosome", "Start", "End"], keep=False
)

df = primary_annos.query('feat_class == "Promoter" & is_duplicated')[
    ["Chromosome", "Start", "End", "gtfanno_uid", "gene_name"]
]
display(df.head(20))
display(df.tail(20))

# Nsdhl
# http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000031349;r=X:71962163-72002120
#
# Rpl7
# http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000043716;r=1:16171519-16174886

# + [markdown] tags=[] heading_collapsed="true"
# #### merge annotations
# -

# Merging strategy: keep all
# - for Promoters, the window is relatively small. Ranking on TSS distance in such a small window seems arbitrary.
# - for enhancer candidates, a typical strategy would be to identify all TSS in +-100 kb window and try to find the target through correlation with gene expression, eg PMID: 30686579. So it also makes sense to indicate all genes in the window to give an impression of the number of possible target genes.

# %%time
merged_annos = lib.merge_annos(primary_annos=primary_annos)

merged_annos

merged_annos_new_chrom_dtype = merged_annos.copy()
merged_annos_new_chrom_dtype["Chromosome"] = (
    "chr" + merged_annos["Chromosome"].astype(str)
).astype(chrom_dtype_prefixed)
merged_annos_new_chrom_dtype = (
    merged_annos_new_chrom_dtype.sort_values(["Chromosome", "Start", "End"])
    .drop("gtfanno_uid", axis=1)
    .reset_index(drop=True)
)

merged_annos_final = pd.merge(
    merged_annos_new_chrom_dtype,
    illumina_probes_curated_chrom_defined,
    on=["Chromosome", "Start", "End"],
    how="left",
)

merged_annos_final.head(3)

merged_annos_final.shape

illumina_probes_curated_chrom_defined

assert merged_annos_final["name"].notnull().all()

pd.testing.assert_frame_equal(
    merged_annos_final[["Chromosome", "Start", "End"]],
    illumina_probes_curated_chrom_defined[["Chromosome", "Start", "End"]].astype(
        {"Start": "i8", "End": "i8"}
    ),
)

merged_annos_final.iloc[0]

# + [markdown] heading_collapsed="true" tags=[]
# #### Finalize annotation tables
# -

merged_annos_final.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    gene_annos_primary_one_row, sep="\t", header=True, index=False
)

primary_annos_final = (
    primary_annos.drop("gtfanno_uid", axis=1)
    .assign(
        Chromosome=lambda df: ("chr" + df["Chromosome"].astype(str)).astype(
            chrom_dtype_prefixed
        )
    )
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)

primary_annos_final.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    gene_annos_primary_multi_row, sep="\t", header=True, index=False
)

# !head {gene_annos_primary_multi_row}

# + [markdown] tags=[] heading_collapsed="true"
# ### CpG island annotations
# -

cpg_island_classif_df = lib.classify_cpg_island_overlap(
    granges_df=original_probes_df_sorted,
    cpg_islands_pickle_d=cpg_islands_pickle_d,
)
cpg_island_classif_df.head(3)

# ### Merge all annotations

cpg_island_classif_df
merged_annos_final

# # End


