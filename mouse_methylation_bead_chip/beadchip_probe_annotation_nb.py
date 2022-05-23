# ---
# jupyter:
#   jupytext:
#     formats: py:percent,ipynb,md:myst
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:mouse_hema_meth_dev3] *
#     language: python
#     name: conda-env-mouse_hema_meth_dev3-xpython
# ---

# %% [markdown]
# # TODO

# %% [markdown]
# - check sorting order of chromosomes throughout analysis
# - there is a bug in gtfanno which makes it fail to read the appris principal score, does this influence the results?
# - there is a "transcript = 0" count in the gtfanno basic stat output. is this indicative of a problem?
# - double check the changes made to allow DCRD intervals enveloping the promoter interval
# - check annos in IGV
# - are feature coordinates really 0-based, right open?

# %% [markdown]
# annos to add
# - motif: CG, CHH, CHG
# - strand, illumina strands are not cytosine strands

# %% [markdown]
# # Setup

# %% [markdown]
# ## Resource parameters

# %%
n_cores = 12

# %% [markdown]
# ## Imports

# %%
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

# %%
# %matplotlib inline

# %%
import mouse_methylation_bead_chip.beadchip_probe_annotation_lib as lib
import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths

# %% [markdown]
# ## Rerun flags

# %%
recompute = False

# %% [markdown]
# ## Dtypes

# %%
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

# %% [markdown]
# # Paths

# %%
temp_dir_obj = tempfile.TemporaryDirectory(dir=paths.project_dir)
temp_dir_name = temp_dir_obj.name
temp_dir_name

# %% [markdown]
# # Analysis

# %% [markdown]
# ## Prepare input data

# %% [markdown]
# ### CpG island annos

# %%
import mouse_hema_meth.genome_annotations.get_genome_annos_paths as get_genome_annos_paths

cpg_islands_pickle_d = get_genome_annos_paths.cpg_islands_shores_shelves_pickle_paths_d

# %% [markdown] tags=[] heading_collapsed="true"
#
# ### Prepare gene annotation

# %% [markdown] tags=[] heading_collapsed="true"
#
# #### download gencode

# %% tags=[]
if recompute:
    subprocess.run(
        ["wget", "-O", paths.gencode_gtf, paths.gencode_download_url],
        check=True,
    )

# %% tags=[]
# !zcat {paths.gencode_gtf} | head -n 6

# %% [markdown] tags=[] heading_collapsed="true"
#
# #### Filter and reformat gencode GTF
#

# %% [markdown]
# - restrict to canonical transcripts
# - restrict to coding transcripts
# - remove chr prefix
# - change M to MT

# %%
gencode_df = pr.read_gtf(paths.gencode_gtf, as_df=True, duplicate_attr=True)

# %%
# extract appris principal score from tags
appris_principal_score = (
    gencode_df["tag"].str.extract(r"appris_principal_(\d)", expand=False).astype(float)
)

# %%
appris_principal_score.value_counts()

# %%
appris_principal_score.isnull().sum()

# %%
appris_principal_score.notnull().sum()

# %%
is_principal_transcript = appris_principal_score.notnull()

# %%
is_protein_coding = gencode_df["gene_type"].eq("protein_coding")

# %%
gencode_df_coding_canonical = gencode_df.loc[
    is_principal_transcript & is_protein_coding
].copy()

# %%
gencode_df_coding_canonical.head(3)

# %%
gencode_df_coding_canonical.shape

# %%
gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].str.replace("chr", "")
gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].replace("M", "MT")

# %%
gencode_pr = pr.PyRanges(gencode_df_coding_canonical)
gencode_pr.df.Chromosome.unique()

# %%
gencode_pr.to_gtf(paths.gencode_coding_canonical_gtf)

# %%
# !zcat {paths.gencode_coding_canonical_gtf} | head

# %% [markdown]
# verify gtf

# %% tags=[] jupyter={"outputs_hidden": true}
# !zcat {paths.gencode_coding_canonical_gtf} | grep ^protein_coding

# %% tags=[] jupyter={"outputs_hidden": true}
# !zcat {paths.gencode_coding_canonical_gtf} | grep ^appris

# %% [markdown] tags=[]
# ### Prepare and inspect probes files

# %% [markdown] tags=[] heading_collapsed="true"
# #### Illumina probe file

# %% [markdown]
# ##### Download

# %%
if recompute:
    subprocess.run(
        ["wget", "-O", paths.illumina_probes_csv, paths.illumina_probes_url], check=True
    )

# %%
# !head {paths.illumina_probes_csv}

# %% [markdown]
# ##### Get curated BED intervals for probes

# %%
illumina_probes = pd.read_csv(
    paths.illumina_probes_csv,
    skiprows=7,
    dtype={
        "AddressA_ID": str,
        "CHR": str,
        "MFG_Change_Flagged": "boolean",
        "MAPINFO": "Int64",
    },
)

# %% [markdown]
# exemplary row, drop fields with longish sequence strings for display

# %%
illumina_probes.drop(["Forward_Sequence", "Top_Sequence"], axis=1).iloc[0].to_frame()

# %% [markdown]
# There are CHR in ['0', nan] probes

# %%
illumina_probes.CHR.value_counts(dropna=False)

# %% [markdown]
# example for control probe (everything NA)

# %%
illumina_probes.loc[illumina_probes.CHR.isnull()].drop(
    ["Forward_Sequence", "Top_Sequence"], axis=1
).iloc[0].to_frame()

# %% [markdown]
# example for chromosome 0 probe

# %%
illumina_probes.loc[illumina_probes.CHR.eq("0")].drop(
    ["Forward_Sequence", "Top_Sequence"], axis=1
).iloc[0].to_frame()

# %% [markdown]
# checked manually in my index files: 1-based Start info is in MAPINFO

# %% [markdown]
# - add 'chr' prefix and make Categorical
# - provide BED interval for the 1-bp long cytosine intervals
# - restrict to columns Chromosome, Start, End

# %%
illumina_probe_intervals_bed_convention = (
    illumina_probes[["CHR", "MAPINFO", "IlmnID"]]
    .rename(columns={"CHR": "Chromosome", "MAPINFO": "Start", "IlmnID": "name"})
    .assign(
        Start=lambda df: df["Start"] - 1,
        End=lambda df: df["Start"] + 1,
        Chromosome=lambda df: ("chr" + df["Chromosome"]).astype(chrom_dtype_prefixed),
    )[["Chromosome", "Start", "End", "name"]]
)


illumina_probes_curated_chrom_defined = (
    illumina_probe_intervals_bed_convention.loc[
        lambda df: df.Chromosome.notnull() & df.Chromosome.ne("0")
    ]
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)

# %% [markdown]
# assert that we have only lost the chromosome 0|na probes

n_probes_chr_defined = (
    illumina_probes.shape[0] - illumina_probes.CHR.isin(["0", np.nan]).sum()
)
assert n_probes_chr_defined == illumina_probes_curated_chrom_defined.shape[0]

# %% [markdown]
# There are duplicated intervals

# %%
illumina_probes_curated_chrom_defined[["Chromosome", "Start", "End"]].duplicated().sum()

# %% [markdown]
# to get unique intervals for annotation with gtfanno
# - drop duplicate rows
# - remove prefix
# - change to alphabetic sorting order
# - only save the grange columns

# %%
illumina_probes_curated_chrom_defined.assign(
    Chromosome=lambda df: df.Chromosome.astype(str).str.replace("chr", "")
).iloc[:, 0:3].sort_values(["Chromosome", "Start", "End"]).drop_duplicates().to_csv(
    paths.illumina_coordinate_bed, sep="\t", header=False, index=False
)

# %%
# !head {paths.illumina_coordinate_bed}

# %% [markdown]
# ## Annotation

# %% [markdown]
# ### Gene annotation

# %% [markdown]
# #### Perform annotation

# %%
# %%time
ga.annotate(
    query_bed=paths.illumina_coordinate_bed,
    gtf_fp=paths.gencode_coding_canonical_gtf,
    trunk_path=paths.custom_intervals_trunk_path,
    tmpdir=temp_dir_name,
    promoter=(-1500, 500),
    distant_cis_regulatory_domain=(-100_000, 100_000),
)

# %% [markdown]
# #### Inspect annotations

# %%
primary_annos = pd.read_pickle(
    paths.custom_intervals_results_paths_d["primary_annos_p"]
)

# %%
primary_annos.shape

# %% [markdown]
# ##### General checks
#

# %%
primary_annos.query('feat_class == "Promoter"').head(3)

# %%
primary_annos.query('feat_class == "exon"').head(3)

# %% [markdown]
# ##### Multiple assignments per region

# %% [markdown]
# ###### How is this distributed across feature classes?


# %%
multi_annos_crosstab = (
    primary_annos.groupby(["feat_class", "gtfanno_uid"], observed=True)
    .size()
    .groupby("feat_class")
    .value_counts()
    .unstack()
)
multi_annos_crosstab

# %% [markdown]
# ###### Example for Promoter multiple annotations

# %%
primary_annos["is_duplicated"] = primary_annos.duplicated(
    subset=["Chromosome", "Start", "End"], keep=False
)

# %%
df = primary_annos.query('feat_class == "Promoter" & is_duplicated')[
    ["Chromosome", "Start", "End", "gtfanno_uid", "gene_name"]
]
display(df.head(20))
display(df.tail(20))

# %% [markdown]
# TODO: exchange with appropriate ensembl version
# Nsdhl
# http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000031349;r=X:71962163-72002120
#
# Rpl7
# http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000043716;r=1:16171519-16174886

# %% [markdown]
# #### merge annotations

# %% [markdown]
# Merging strategy: keep all
# - for Promoters, the window is relatively small. Ranking on TSS distance in such a small window seems arbitrary.
# - for enhancer candidates, a typical strategy would be to identify all TSS in +-100 kb window and try to find the target through correlation with gene expression, eg PMID: 30686579. So it also makes sense to indicate all genes in the window to give an impression of the number of possible target genes.

# %%
# %%time
merged_annos = lib.merge_annos(primary_annos=primary_annos)

# %%
merged_annos

# %% [markdown]
# - change back to rnbeads chrom dtype (with prefix and natural sorting order)
# - sort according to new chrom dtype

# %%
merged_annos_new_chrom_dtype = merged_annos.copy()
merged_annos_new_chrom_dtype["Chromosome"] = (
    "chr" + merged_annos["Chromosome"].astype(str)
).astype(chrom_dtype_prefixed)
merged_annos_new_chrom_dtype = (
    merged_annos_new_chrom_dtype.sort_values(["Chromosome", "Start", "End"])
    .drop("gtfanno_uid", axis=1)
    .reset_index(drop=True)
)
assert merged_annos_new_chrom_dtype["Chromosome"].notnull().all()

# %% [markdown]
# merge gene annos into illumina probe df

# %%
merged_annos_final = pd.merge(
    illumina_probe_intervals_bed_convention,
    merged_annos_new_chrom_dtype.astype({"Start": "Int64", "End": "Int64"}),
    on=["Chromosome", "Start", "End"],
    how="left",
)

# %%
merged_annos_final.head(3)

# %%
assert merged_annos_final["name"].notnull().all()

# %%
pd.testing.assert_frame_equal(
    merged_annos_final[["Chromosome", "Start", "End", "name"]],
    illumina_probe_intervals_bed_convention[["Chromosome", "Start", "End", "name"]],
)

# %% [markdown]
# #### Finalize annotation tables

# %%
# assert against older version of annotation (which left out CHR in ['0', nan] probes)

# %%
gene_annos_primary_one_row_v1 = pd.read_csv(
    paths.gene_annos_primary_one_row,
    sep="\t",
    keep_default_na=False,
    na_values="",
    dtype={
        "#Chromosome": chrom_dtype_prefixed,
        "feat_class": merged_annos_final["feat_class"].dtype,
        "feat_chrom": str,
        "center": "f8",
        "Start": "Int64",
        "End": "Int64",
    },
)

# %%
cols = [
    "#Chromosome",
    "Start",
    "End",
    "center",
    "feat_class",
    "perc_feature",
    "perc_region",
    "distance",
    "has_center",
    "gene_name",
    "gene_id",
    "transcript_id",
    "appris_principal_score",
    "feat_chrom",
    "feat_start",
    "feat_end",
    "feat_center",
    "feat_strand",
    "feature_rank",
    "name",
]
pd.testing.assert_frame_equal(
    merged_annos_final.rename(columns={"Chromosome": "#Chromosome"})[cols]
    .sort_values(["#Chromosome", "Start", "End"])
    .loc[lambda df: df["#Chromosome"].notnull()]
    .reset_index(drop=True),
    gene_annos_primary_one_row_v1.drop("is_duplicated", axis=1)[cols]
    .sort_values(["#Chromosome", "Start", "End"])
    .reset_index(drop=True),
)

# %%
merged_annos_final.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    paths.gene_annos_primary_one_row_v2, sep="\t", header=True, index=False
)

# %% [markdown]
# ### CpG island annotations

# %% [markdown]
# This uses the mm10 cpg islands from UCSC


# %%
# two probes at same coordinate, different name
# remove first, merge again later
illumina_probes_curated_chrom_defined_no_names_no_dup_intervals = (
    illumina_probes_curated_chrom_defined[
        ["Chromosome", "Start", "End"]
    ].drop_duplicates(subset=["Chromosome", "Start", "End"])
)
granges_gr = pr.PyRanges(
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals
)

# %%
# join against cpg islands, shelves, shores
dfs = []
for cpg_island_region_name, cpg_island_region_p in cpg_islands_pickle_d.items():
    cpg_island_region_df = pd.read_pickle(cpg_island_region_p).assign(
        region_name=cpg_island_region_name,
    )
    cpg_island_region_gr = pr.PyRanges(cpg_island_region_df)
    join_gr = granges_gr.join(
        cpg_island_region_gr,
        how=None,
    )
    join_df = join_gr.df.rename(
        columns={"Start_b": "region_start", "End_b": "region_end"}
    )[
        [
            "Chromosome",
            "Start",
            "End",
            "region_name",
            "region_start",
            "region_end",
        ]
    ]
    dfs.append(join_df)
# NOTE: pyranges may change sorting order
cpg_island_anno_df = pd.concat(dfs, axis=0).reset_index(drop=True)

# %% [markdown]
# some cpgs are part of multiple features, e.g. cpg is located in a CpG island and in the shelve region of another CpG island; or a CpG is located in two shelve regions of different CpG islands

# %%
grange_cols = ["Chromosome", "Start", "End"]
cpg_island_anno_df.groupby(grange_cols, observed=True).size().value_counts()

cpg_island_anno_df.set_index(grange_cols).loc[
    cpg_island_anno_df.groupby(grange_cols, observed=True).size().gt(2)
].sort_values(grange_cols)

# %% [markdown]
# Sanity check: no CpG is annotated to two CpG islands

# %%
cpg_island_anno_df.query('region_name == "CpG islands"').groupby(
    grange_cols, observed=True
).size().value_counts()

# %% [markdown]
# pick unique classif according to precedence

# %%
cpg_island_anno_df_unique = (
    cpg_island_anno_df.assign(
        region_name=lambda df: pd.Categorical(
            df["region_name"],
            categories=["CpG islands", "CpG shelves", "CpG shores", "open sea"],
        )
    )
    .sort_values(["Chromosome", "Start", "End", "region_name"], ascending=True)
    .groupby(grange_cols, observed=True, as_index=False)
    .first()
    .assign(
        Chromosome=lambda df: df["Chromosome"].astype(
            illumina_probes_curated_chrom_defined.Chromosome.dtype
        )
    )
)

# %% [markdown]
# the grange join operation has discarded open sea cpgs, merge to get them back

# %%
full_cpg_island_anno_df = pd.merge(
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals,
    cpg_island_anno_df_unique,
    on=grange_cols,
    how="left",
)

# assert: the open sea cpgs have no region_name yet, all other features do have a region_name
assert (
    full_cpg_island_anno_df.region_name.notnull().sum()
    == cpg_island_anno_df_unique.shape[0]
)

full_cpg_island_anno_df["region_name"] = full_cpg_island_anno_df["region_name"].fillna(
    "open sea"
)

# add strand classification
nearest_cpg_island_grange = granges_gr.nearest(
    pr.PyRanges(pd.read_pickle(cpg_islands_pickle_d["CpG islands"])),
    strandedness=False,
    overlap=True,
    how=None,
)

# Distance computed by pyranges has no sign
nearest_cpg_island_df = nearest_cpg_island_grange.df.copy()
is_upstream = nearest_cpg_island_df.eval("Start < Start_b")
is_downstream = nearest_cpg_island_df.eval("Start >= End_b")
nearest_cpg_island_df["distance_signed"] = 0
nearest_cpg_island_df.loc[is_upstream, "distance_signed"] = nearest_cpg_island_df.loc[
    is_upstream
].eval("(End - 1) - Start_b")
nearest_cpg_island_df.loc[is_downstream, "distance_signed"] = nearest_cpg_island_df.loc[
    is_downstream
].eval("Start - (End_b - 1)")

pd.testing.assert_series_equal(
    nearest_cpg_island_df.distance_signed.abs(),
    nearest_cpg_island_df.Distance,
    check_names=False,
)

full_cpg_island_anno_df_with_dist = pd.merge(
    full_cpg_island_anno_df,
    nearest_cpg_island_df[grange_cols + ["distance_signed"]],
    on=grange_cols,
    how="left",
)
full_cpg_island_anno_df_with_dist.Chromosome = (
    full_cpg_island_anno_df_with_dist.Chromosome.astype(
        illumina_probes_curated_chrom_defined.Chromosome.dtype
    )
)
full_cpg_island_anno_df_with_dist["north_south_of_island"] = np.sign(
    full_cpg_island_anno_df_with_dist["distance_signed"]
).replace({0: "", 1: "N", -1: "S"})

full_cpg_island_anno_df_with_dist.loc[
    ~full_cpg_island_anno_df_with_dist.region_name.isin(["CpG shores", "CpG shelves"]),
    "north_south_of_island",
] = ""

# check that cpg islands and open sea have no north / south, but shelves and shores do have it
# check that strands are approx 50/50
pd.crosstab(
    full_cpg_island_anno_df_with_dist["region_name"],
    full_cpg_island_anno_df_with_dist["north_south_of_island"],
)

pd.testing.assert_frame_equal(
    full_cpg_island_anno_df_with_dist[grange_cols],
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals,
)


full_cpg_island_anno_df_with_dist["distance_signed"].plot.hist(
    bins=np.linspace(-1e5, 1e5, 100)
)


# %% [markdown]
# ### Add cytosine motif and strand

# %%

motif_index_grange_cols = ["Chromosome", "Start", "End", "motif", "strand"]
annotated_chrom_dfs_d = {}
motif_dtype = pd.CategoricalDtype(categories=["CG", "CHG", "CHH"], ordered=False)
strand_dtype = pd.CategoricalDtype(categories=["+", "-"], ordered=False)
for (
    chrom,
    chrom_df,
) in illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.groupby(
    "Chromosome"
):
    print(chrom)
    chrom_no_prefix = chrom.replace("chr", "")
    chrom_index_bed_gz = f"/omics/groups/OE0029/internal/kraemers/projects/mbias/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chrom_no_prefix}.bed.gz"
    chrom_index_df = pd.read_csv(
        chrom_index_bed_gz,
        sep="\t",
        usecols=[0, 1, 2, 3, 5],
        comment="#",
        header=None,
        names=motif_index_grange_cols,
        dtype={
            "Chromosome": str,
            "Start": "Int64",
            "End": "Int64",
            "motif": motif_dtype,
            "strand": strand_dtype,
        },
    )
    assert chrom_index_df.notnull().all().all()
    chrom_index_df["Chromosome"] = ("chr" + chrom_index_df["Chromosome"]).astype(
        chrom_df["Chromosome"].dtype
    )
    chrom_df_annotated = pd.merge(
        chrom_df,
        chrom_index_df,
        on=["Chromosome", "Start", "End"],
        how="left",
    )
    # chrom_df_annotated[["motif", "strand"]].notnull().all().all()
    # chrom_df_annotated.loc[chrom_df_annotated["motif"].isnull()]
    # chrom_df.shape
    pd.testing.assert_frame_equal(
        chrom_df.reset_index(drop=True), chrom_df_annotated[["Chromosome", "Start", "End"]].reset_index(drop=True)
    )
    annotated_chrom_dfs_d[chrom] = chrom_df_annotated

res = (
    pd.concat(annotated_chrom_dfs_d)
    # .assign(Chromosome=lambda df: df["Chromosome"].astype(chrom_dtype))
    # .sort_values(motif_index_grange_cols)
    # .reset_index(drop=True)
)
res.dtypes
res.Chromosome


# started: 19.40
if recompute:
    subprocess.run(["wget", "-O", paths.mm10_fa, paths.mm10_fa_url], check=True)
    subprocess.run(['samtools', 'faidx', paths.mm10_fa], check=True)

# add to env: samtools faidx
from pathlib import Path

regions_txt = temp_dir_name + "/regions.txt"
df = res.loc[res["motif"].isnull()]
Path(regions_txt).write_text((
    df["Chromosome"].astype(str).str.replace("chr", "")
    + ":"
    + df["Start"].add(1).astype(str)
    + "-"
    + df["End"].astype(str)
).str.cat(sep="\n"))
proc = subprocess.run(
    [
        "samtools",
        "faidx",
        '-r',
        regions_txt,
    ],
    check=True,
    capture_output=True,
)
assert 'C' not in ''.join(proc.stdout)
assert 'G' not in ''.join(proc.stdout)

# %%
print(
    df.shape[0],
    res.shape[0],
    )





# %% [markdown]
# ### Merge all annotations

# %%
cpg_island_classif_df
merged_annos_final

# %% [markdown]
# # End

# %%
