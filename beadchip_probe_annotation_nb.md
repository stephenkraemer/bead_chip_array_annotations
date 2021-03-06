---
jupytext:
  formats: py:percent,ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.1
kernelspec:
  display_name: Python [conda env:mouse_hema_meth_dev3] *
  language: python
  name: conda-env-mouse_hema_meth_dev3-xpython
---

# How to reproduce this analysis

+++

## Install the environment

+++

Environment yaml and lock file are distributed with the notebook

```
mamba env create -f /path/to/probe-anno_env_locked.yaml
or for development
mamba env create -f /path/to/probe-anno_env.yaml
```

+++

## Create a project directory and copy the manifest file into it

+++

- create a project directory
- copy the manifest file somewhere into this project directory
- specify these paths in the Config section

+++

## Specify the available number of cores

+++

The number of cores can be specified in the Config section as well

+++

# Setup

+++

## Config (ADAPT)

```{code-cell}
# TODO: set number of available cores as needed
n_cores = 12
```

```{code-cell}
# TODO: set basic paths as needed
project_dir = (
    "/omics/groups/OE0029/internal/kraemers/projects/epic-arrays-for-hematopoiesis"
)

# this file is distributed with the notebook, place it somewhere in project_dir
# and note the location here
illumina_probes_csv = (
    project_dir
    + "/Infinium_20Mouse_20Methylation_20v1.0_20A1_20GS_20Manifest_20File.csv"
)
```

## Imports

```{code-cell}
import os

num_threads = str(n_cores)

# these need to be set prior to numpy import
os.environ["OMP_NUM_THREADS"] = num_threads
os.environ["OPENBLAS_NUM_THREADS"] = num_threads
os.environ["MKL_NUM_THREADS"] = num_threads
os.environ["VECLIB_MAXIMUM_THREADS"] = num_threads
os.environ["NUMEXPR_NUM_THREADS"] = num_threads

import numpy as np

import subprocess
import tempfile
from pathlib import Path
from typing import Dict

import gtfanno as ga
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from IPython.display import display
```

## Rerun flags

```{code-cell}
recompute = True
```

## Paths and URLs

```{code-cell}
# ### Input data urls and filepaths
```

#### Reference genomes

```{code-cell}
mm10_fa_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz'
mm10_fa = project_dir + '/ucsc_mm10.fa.gz'
mm10_fa_bgz = project_dir + '/ucsc_mm10.fa.bgz'
```

```{code-cell}
ncbi_mm10_fa_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/all_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz'
ncbi_mm10_fa = project_dir + '/ncbi_mm10.fa.gz'
ncbi_mm10_fa_bgz = project_dir + '/ncbi_mm10.fa.bgz'
```

```{code-cell}
chrom_alias_url = 'https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromAlias.txt.gz'
chrom_alias_txt_gz = project_dir + '/chromAlias.txt.gz'
```

#### Gencode

```{code-cell}
gencode_download_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

# unmodified gencode gtf (note: with chr prefix)
gencode_gtf = project_dir + "/gencode.vM25.annotation.gtf.gz"

# gencode filtered for principal transcripts of protein coding genes, note that that chromosome prefix ('chr') is removed in this file
gencode_coding_canonical_gtf = (
    project_dir + "/gencode.vM25.annotation_coding_canonical.gtf.gz"
)
```

#### CpG islands

```{code-cell}
cpg_islands_ucsc_unmasked_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cpgIslandExtUnmasked.txt.gz'
cpg_islands_ucsc_unmasked_txt_gz = project_dir + "/cpg-islands_ucsc_unmasked.txt.gz"
```

#### Probes

+++

originally retrieved from:
illumina_probes_url = "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/Infinium%20Mouse%20Methylation%20v1.0%20A1%20GS%20Manifest%20File.csv"
saved to repo in case the link is not stable, use from here in the future (see Config section)

+++

### Output file paths

+++

#### gtfanno results

```{code-cell}
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
```

#### Final tables with CpG island and motif annotations

```{code-cell}
probe_annos_one_row_bed_csv = project_dir + "/gene-annos_primary_one-row_v2.bed"
probe_annos_one_row_bed_parquet = project_dir + "/gene-annos_primary_one-row_v2.parquet"
```

## TemporaryDirectory

```{code-cell}
temp_dir_obj = tempfile.TemporaryDirectory(dir=project_dir)
temp_dir_name = temp_dir_obj.name
temp_dir_name
```

# Functions and objects
## Dtypes

```{code-cell}
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
```

## Merge multiple gene annos per probe

```{code-cell}
def merge_annos(
    primary_annos: pd.DataFrame,
    #     chrom_dtype: pd.api.types.CategoricalDtype,
) -> pd.DataFrame:
    """Merge annotation format with multiple rows per DMR"""

    # TODO: could use some print statements

    # Merge information as CSV for fields with different values for different annotations
    fields_with_different_values = [
        "perc_feature",
        "perc_region",
        "distance",
        "has_center",
        "gene_name",
        "gene_id",
        "transcript_id",
        "appris_principal_score",
        "feat_start",
        "feat_end",
        "feat_center",
        "feat_strand",
    ]

    fields_with_unique_values = [
        f for f in primary_annos.columns if f not in fields_with_different_values
    ]

    # float fields require special str formatting
    float_fields = [
        f
        for f in fields_with_different_values
        if pd.api.types.is_float_dtype(primary_annos.dtypes[f])
    ]

    # When we concatenate, we need to convert all fields to str/object dtype
    primary_annos_new_dtypes = primary_annos.copy()
    primary_annos_new_dtypes[float_fields] = primary_annos_new_dtypes[
        float_fields
    ].round(2)
    primary_annos_new_dtypes[fields_with_different_values] = primary_annos_new_dtypes[
        fields_with_different_values
    ].astype(str)

    # computations are expensive, so only perform on features which are actually duplicated
    # always take gtfanno_uid into the mix so that we can avoid problems with duplicate regions
    # (which occur for example on the EPIC array)
    is_duplicated = primary_annos_new_dtypes.duplicated(
        subset=["Chromosome", "Start", "End", "gtfanno_uid"], keep=False
    )
    duplicated_df = primary_annos_new_dtypes.loc[is_duplicated].copy()

    #     primary_annos.query('feat_class == "Promoter"')

    print("merge unique value fields")
    merged_different_value_fields_df = (
        duplicated_df
        # for testing
        # .query('feat_class == "Promoter"')
        # .head(50)
        # / for testing
        .groupby(
            ["Chromosome", "Start", "End", "gtfanno_uid"], observed=True, as_index=False
        )[fields_with_different_values].aggregate(_agg_multi_value_columns)
    )
    merged_different_value_fields_df

    # Assert that unique value fields are really unique (1 or 0 unique values if nan)

    duplicated_df

    assert (
        duplicated_df
        # .head(100)
        .groupby(["Chromosome", "Start", "End", "gtfanno_uid"], observed=True)[
            fields_with_unique_values
        ]
        .agg(lambda ser: ser.nunique())
        # bug: nunique() does not respected observed atm
        # .nunique()
        .isin([0, 1])
        .all()
        .all()
    )

    merged_unique_value_fields_df = duplicated_df.groupby(
        ["Chromosome", "Start", "End", "gtfanno_uid"],
        observed=True,
        as_index=False,
    )[fields_with_unique_values].first()
    merged_unique_value_fields_df

    merged_duplicated_annos = pd.concat(
        [
            merged_different_value_fields_df.drop(
                ["Chromosome", "Start", "End", "gtfanno_uid"], axis=1
            ),
            merged_unique_value_fields_df,
        ],
        axis=1,
    )[primary_annos.columns]
    merged_duplicated_annos

    full_merged_annos = pd.concat(
        [primary_annos_new_dtypes.loc[~is_duplicated], merged_duplicated_annos], axis=0
    ).sort_values(["Chromosome", "Start", "End"])
    full_merged_annos

    assert (
        full_merged_annos["gtfanno_uid"] == np.arange(full_merged_annos.shape[0])
    ).all()

    return full_merged_annos.reset_index(drop=True)


def _agg_multi_value_columns(ser):
    if ser.eq("nan").all():
        return "nan"
    else:
        return ser.str.cat(sep=",")
```

## Classify cytosine motifs

```{code-cell}
def classify_motif(s):
    if s[2] == "C":
        if s[3] == "G":
            return "CG"
        elif s[3] == "N":
            return "CN"
        elif s[4] == "G":
            return "CHG"
        elif s[4] == "N":
            return "CHN"
        else:
            return "CHH"
    elif s[2] == "G":
        if s[1] == "C":
            return "CG"
        elif s[1] == "N":
            return "CN"
        elif s[0] == "C":
            return "CHG"
        elif s[0] == "N":
            return "CHN"
        else:
            return "CHH"
    else:
        return "D"
```

```{code-cell}
def find_cytosine_strand_and_motif(df, temp_dir_name, fa_bgz):
    """

    chrom names must match ref genome
    """

    regions_txt = temp_dir_name + "/regions.txt"
    Path(regions_txt).write_text(
        (
            df["Chromosome"].astype(str)
            + ":"
            + df["Start"].add(1 - 2).astype(str)
            + "-"
            + df["End"].add(2).astype(str)
        ).str.cat(sep="\n")
    )


    res_txt = temp_dir_name + "/res.txt"
    proc2 = subprocess.run(
        [
            "samtools",
            "faidx",
            "-r",
            regions_txt,
            "-o",
            res_txt,
            fa_bgz,
        ],
        check=True,
        capture_output=True,
        encoding="utf-8",
    )

    bases_ser = pd.Series(
        [s for s in Path(res_txt).read_text().upper().split() if not s.startswith(">")]
    )

    strand_ser = bases_ser.str.slice(2, 3).map({"C": "+", "G": "-"}).fillna("NA")
    motifs_ser = bases_ser.map(classify_motif)

    return strand_ser, motifs_ser
```

# Analysis

+++

## Prepare and inspect the probe manifest

```{code-cell}
!head {illumina_probes_csv}
```

### Get curated BED intervals for probes

```{code-cell}
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
```

exemplary row, drop fields with longish sequence strings for display

```{code-cell}
illumina_probes.drop(["Forward_Sequence", "Top_Sequence"], axis=1).iloc[0].to_frame()
```

There are CHR in ['0', nan] probes

```{code-cell}
illumina_probes.CHR.value_counts(dropna=False)
```

example for control probe (everything NA)

```{code-cell}
illumina_probes.loc[illumina_probes.CHR.isnull()].drop(
    ["Forward_Sequence", "Top_Sequence"], axis=1
).iloc[0].to_frame()
```

example for chromosome 0 probe

```{code-cell}
illumina_probes.loc[illumina_probes.CHR.eq("0")].drop(
    ["Forward_Sequence", "Top_Sequence"], axis=1
).iloc[0].to_frame()
```

checked manually in my index files: MAPINFO is 1-based Start info

+++

- add 'chr' prefix and make Categorical
- provide BED interval for the 1-bp long cytosine intervals
- restrict to columns Chromosome, Start, End

```{code-cell}
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
```

assert that we have only lost the chromosome 0|na probes

```{code-cell}
n_probes_chr_defined = (
    illumina_probes.shape[0] - illumina_probes.CHR.isin(["0", np.nan]).sum()
)
assert n_probes_chr_defined == illumina_probes_curated_chrom_defined.shape[0]
```

There are duplicated intervals, ie same interval, different probes, with different names

```{code-cell}
illumina_probes_curated_chrom_defined[["Chromosome", "Start", "End"]].duplicated().sum()
```

to get unique intervals for annotation with gtfanno
- drop duplicate rows
- remove prefix
- change to alphabetic sorting order
- only save the grange columns

```{code-cell}
illumina_coordinate_bed = temp_dir_name + '/illumina-coordinates.bed'
```

```{code-cell}
illumina_probes_curated_chrom_defined.assign(
    Chromosome=lambda df: df.Chromosome.astype(str).str.replace("chr", "")
).iloc[:, 0:3].sort_values(["Chromosome", "Start", "End"]).drop_duplicates().to_csv(
    illumina_coordinate_bed, sep="\t", header=False, index=False
)
```

```{code-cell}
!head {illumina_coordinate_bed}
```

```{code-cell}
illumina_probes_curated_chrom_defined_no_names_no_dup_intervals = (
    illumina_probes_curated_chrom_defined[
        ["Chromosome", "Start", "End"]
    ].drop_duplicates(subset=["Chromosome", "Start", "End"])
).reset_index(drop=True)
```

## Gene annotation

+++

### Prepare gene annotation

+++ {"tags": [], "heading_collapsed": "true"}

#### download gencode

```{code-cell}
:tags: []

if recompute:
    subprocess.run(
        ["wget", "-O", gencode_gtf, gencode_download_url],
        check=True,
        capture_output=True,
    )
```

#### Filter and reformat gencode GTF

+++

- restrict to canonical transcripts
- restrict to coding transcripts
- remove chr prefix
- change M to MT

```{code-cell}
gencode_df = pr.read_gtf(gencode_gtf, as_df=True, duplicate_attr=True)
```

```{code-cell}
# extract appris principal score from tags
appris_principal_score = (
    gencode_df["tag"].str.extract(r"appris_principal_(\d)", expand=False).astype(float)
)
```

```{code-cell}
appris_principal_score.value_counts()
```

```{code-cell}
appris_principal_score.isnull().sum()
```

```{code-cell}
appris_principal_score.notnull().sum()
```

```{code-cell}
is_principal_transcript = appris_principal_score.notnull()
```

```{code-cell}
is_protein_coding = gencode_df["gene_type"].eq("protein_coding")
```

```{code-cell}
gencode_df_coding_canonical = gencode_df.loc[
    is_principal_transcript & is_protein_coding
].copy()
```

```{code-cell}
gencode_df_coding_canonical.head(3)
```

```{code-cell}
gencode_df_coding_canonical.shape
```

```{code-cell}
gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].str.replace("chr", "")
gencode_df_coding_canonical["Chromosome"] = gencode_df_coding_canonical[
    "Chromosome"
].replace("M", "MT")
```

```{code-cell}
gencode_pr = pr.PyRanges(gencode_df_coding_canonical)
gencode_pr.df.Chromosome.unique()
```

```{code-cell}
gencode_pr.to_gtf(gencode_coding_canonical_gtf)
```

verify gtf

```{code-cell}
!zcat {gencode_coding_canonical_gtf} | grep ^protein_coding
```

```{code-cell}
!zcat {gencode_coding_canonical_gtf} | grep ^appris
```

### Perform gene annotations

```{code-cell}
if recompute:
    ga.annotate(
        query_bed=illumina_coordinate_bed,
        gtf_fp=gencode_coding_canonical_gtf,
        trunk_path=custom_intervals_trunk_path,
        tmpdir=temp_dir_name,
        promoter=(-1500, 500),
        distant_cis_regulatory_domain=(-100_000, 100_000),
    )
```

```{code-cell}
primary_annos = pd.read_pickle(
    custom_intervals_results_paths_d["primary_annos_p"]
)
```

#### Multiple assignments per region

+++

##### How is this distributed across feature classes?

```{code-cell}
multi_annos_crosstab = (
    primary_annos.groupby(["feat_class", "gtfanno_uid"], observed=True)
    .size()
    .groupby("feat_class")
    .value_counts()
    .unstack()
)
multi_annos_crosstab
```

##### Example for Promoter multiple annotations

```{code-cell}
primary_annos["is_duplicated"] = primary_annos.duplicated(
    subset=["Chromosome", "Start", "End"], keep=False
)
```

```{code-cell}
df = primary_annos.query('feat_class == "Promoter" & is_duplicated')[
    ["Chromosome", "Start", "End", "gtfanno_uid", "gene_name"]
]
display(df.head(20))
display(df.tail(20))
```

Nsdhl
http://nov2020.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000031349;r=X:71962163-72002120

Rpl7
http://nov2020.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000043716;r=1:16171519-16174886

+++

### Merge gene annotations (one row per unique probe interval)

+++

Merging strategy: keep all
- for Promoters, the window is relatively small. Ranking on TSS distance in such a small window seems arbitrary.
- for enhancer candidates, a typical strategy would be to identify all TSS in +-100 kb window and try to find the target through correlation with gene expression, eg PMID: 30686579. So it also makes sense to indicate all genes in the window to give an impression of the number of possible target genes.

```{code-cell}
merged_annos = merge_annos(primary_annos=primary_annos)
merged_annos
```

```{code-cell}
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
```

```{code-cell}
merged_annos_new_chrom_dtype = merged_annos_new_chrom_dtype.astype(
    {"Start": "Int64", "End": "Int64"}
)
```

## CpG island annotations

+++

### Download CpG island regions from UCSC

```{code-cell}
if recompute:
    subprocess.run(
        [
            "wget",
            "-O",
            cpg_islands_ucsc_unmasked_txt_gz,
            cpg_islands_ucsc_unmasked_url,
        ],
        check=True,
        capture_output=True,
    )
```

### Read in CpG islands

```{code-cell}
cpg_islands_df_prelim = pd.read_csv(  # type: ignore
        cpg_islands_ucsc_unmasked_txt_gz,
        sep="\t",
        header=None,
        names="bin Chromosome Start End name length cpgNum gcNum perCpg perGc obsExp".split(),
    )
```

no MT cpg islands

```{code-cell}
print(
    "chrMT" in cpg_islands_df_prelim.Chromosome.unique(),
    "chrM" in cpg_islands_df_prelim.Chromosome.unique(),
)
```

```{code-cell}
cpg_islands_df = (cpg_islands_df_prelim
    .astype({"Chromosome": chrom_dtype_prefixed})
    .dropna()
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)
```

### Compute shores and shelves

```{code-cell}
cpg_upstream_shores = cpg_islands_df[["Chromosome", "Start", "End"]].copy()
cpg_upstream_shores["End"] = cpg_islands_df["Start"]
cpg_upstream_shores["Start"] = cpg_islands_df["Start"] - 2000
cpg_downstream_shores = cpg_islands_df[["Chromosome", "Start", "End"]].copy()
cpg_downstream_shores["Start"] = cpg_islands_df["End"]
cpg_downstream_shores["End"] = cpg_islands_df["End"] + 2000
cpg_shores_df = (
    pd.concat([cpg_upstream_shores, cpg_downstream_shores], axis=0)
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
    .astype({"Chromosome": chrom_dtype_prefixed})
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)
```

```{code-cell}
cpg_upstream_shelves = cpg_islands_df[["Chromosome", "Start", "End"]].copy()
cpg_upstream_shelves["End"] = cpg_islands_df["Start"] - 2000
cpg_upstream_shelves["Start"] = cpg_islands_df["Start"] - 4000
cpg_downstream_shelves = cpg_islands_df[["Chromosome", "Start", "End"]].copy()
cpg_downstream_shelves["Start"] = cpg_islands_df["End"] + 2000
cpg_downstream_shelves["End"] = cpg_islands_df["End"] + 4000
cpg_shelves_df = (
    pd.concat([cpg_upstream_shelves, cpg_downstream_shelves], axis=0)
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
    .astype({"Chromosome": chrom_dtype_prefixed})
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)
```

```{code-cell}
cpg_islands_regions_gr_d = {'cpg_island': pr.PyRanges(cpg_islands_df[['Chromosome', 'Start', 'End']]),
                    'cpg_shore': pr.PyRanges(cpg_shores_df),
                    'cpg_shelve': pr.PyRanges(cpg_shelves_df),
                    }
```

### Join probe and cpg region intervals

```{code-cell}
illumina_probes_curated_chrom_defined_no_names_do_dup_intervals_gr = pr.PyRanges(
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals
)
```

```{code-cell}
# join against cpg islands, shelves, shores
dfs = []
for cpg_island_region_name, cpg_island_region_gr in cpg_islands_regions_gr_d.items():
    join_gr = illumina_probes_curated_chrom_defined_no_names_do_dup_intervals_gr.join(
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
            "region_start",
            "region_end",
        ]
    ].assign(region_name=cpg_island_region_name)
    dfs.append(join_df)
# NOTE: pyranges may change sorting order
cpg_island_anno_df = pd.concat(dfs, axis=0).reset_index(drop=True)
```

### Inspect multi-class assignments

+++

some cpgs are part of multiple features, e.g. cpg is located in a CpG island and in the shelve region of another CpG island; or a CpG is located in two shelve regions of different CpG islands

```{code-cell}
cpg_island_grange_cols = ["Chromosome", "Start", "End"]
cpg_island_anno_df.groupby(cpg_island_grange_cols, observed=True).size().value_counts()
```

```{code-cell}
cpg_island_anno_df.set_index(cpg_island_grange_cols).loc[
    cpg_island_anno_df.groupby(cpg_island_grange_cols, observed=True).size().gt(2)
].sort_values(cpg_island_grange_cols)
```

Sanity check: no CpG is annotated to two CpG islands

```{code-cell}
cpg_island_anno_df.query('region_name == "cpg_island"').groupby(
    cpg_island_grange_cols, observed=True
).size().value_counts()
```

### pick unique classif according to precedence

```{code-cell}
cpg_island_anno_df_unique = (
    cpg_island_anno_df.assign(
        region_name=lambda df: pd.Categorical(
            df["region_name"],
            categories=["cpg_island", "cpg_shelve", "cpg_shore", "open_sea"],
        )
    )
    .sort_values(["Chromosome", "Start", "End", "region_name"], ascending=True)
    .groupby(cpg_island_grange_cols, observed=True, as_index=False)
    .first()
    .assign(
        Chromosome=lambda df: df["Chromosome"].astype(
            illumina_probes_curated_chrom_defined.Chromosome.dtype
        )
    )
)
```

### add open sea probes

+++

the grange join operation has discarded open sea cpgs, merge to get them back

```{code-cell}
full_cpg_island_anno_df = pd.merge(
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals,
    cpg_island_anno_df_unique,
    on=cpg_island_grange_cols,
    how="left",
)
```

```{code-cell}
# assert: the open sea cpgs have no region_name yet, all other features do have a region_name
assert (
    full_cpg_island_anno_df.region_name.notnull().sum()
    == cpg_island_anno_df_unique.shape[0]
)
```

```{code-cell}
full_cpg_island_anno_df["region_name"] = full_cpg_island_anno_df["region_name"].fillna(
    "open_sea"
)
```

### add distance to nearest CpG island

```{code-cell}
nearest_cpg_island_grange = illumina_probes_curated_chrom_defined_no_names_do_dup_intervals_gr.nearest(
    cpg_islands_regions_gr_d['cpg_island'],
    strandedness=False,
    overlap=True,
    how=None,
)
```

```{code-cell}
# Distance computed by pyranges has no sign, add it
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
```

### add coords of nearest CpG island to each probe

```{code-cell}
nearest_cpg_island_df = nearest_cpg_island_df.rename(
    columns={"Start_b": "next_cpg_island_start", "End_b": "next_cpg_island_end"}
)

full_cpg_island_anno_df_with_dist = pd.merge(
    full_cpg_island_anno_df,
    nearest_cpg_island_df[
        cpg_island_grange_cols
        + ["distance_signed", "next_cpg_island_start", "next_cpg_island_end"]
    ],
    on=cpg_island_grange_cols,
    how="left",
)
# fix lost chrom dtype (pandas problem)
full_cpg_island_anno_df_with_dist.Chromosome = (
    full_cpg_island_anno_df_with_dist.Chromosome.astype(
        illumina_probes_curated_chrom_defined.Chromosome.dtype
    )
)
```

```{code-cell}
pd.testing.assert_frame_equal(
    full_cpg_island_anno_df_with_dist.query('region_name == "CpG islands"')[
        ["region_start", "region_end"]
    ].set_axis(["a", "b"], axis=1),
    full_cpg_island_anno_df_with_dist.query('region_name == "CpG islands"')[
        ["next_cpg_island_start", "next_cpg_island_end"]
    ].set_axis(["a", "b"], axis=1),
)
```

```{code-cell}
full_cpg_island_anno_df_with_dist = full_cpg_island_anno_df_with_dist.drop(
    ["region_start", "region_end"], axis=1
).rename(
    columns={
        "next_cpg_island_start": "region_start",
        "next_cpg_island_end": "region_end",
    }
)
```

```{code-cell}
# quick sanity check on distance
assert (
    full_cpg_island_anno_df_with_dist.query('region_name == "cpg_island"')[
        "distance_signed"
    ]
    .eq(0)
    .all()
)
```

### Illumina style north/south annos

```{code-cell}
full_cpg_island_anno_df_with_dist["north_south_of_island"] = np.sign(
    full_cpg_island_anno_df_with_dist["distance_signed"]
).replace({0: "", 1: "N", -1: "S"})

full_cpg_island_anno_df_with_dist.loc[
    ~full_cpg_island_anno_df_with_dist.region_name.isin(["cpg_shore", "cpg_shelve"]),
    "north_south_of_island",
] = ""
```

### Final CpG island annos

```{code-cell}
full_cpg_island_anno_df_with_dist = full_cpg_island_anno_df_with_dist.rename(
    columns={
        "region_name": "cpg_island_region_name",
        "region_start": "cpg_island_region_start",
        "region_end": "cpg_island_region_end",
        "distance_signed": "cpg_island_distance",
        "north_south_of_island": "north_south_of_cpg_island",
    }
)
```

### Inspect CpG anno results

```{code-cell}
full_cpg_island_anno_df_with_dist['cpg_island_region_name'].value_counts()
```

```{code-cell}
# check that cpg islands and open sea have no north / south, but shelves and shores do have it
# check that strands are approx 50/50
pd.crosstab(
    full_cpg_island_anno_df_with_dist["cpg_island_region_name"],
    full_cpg_island_anno_df_with_dist["north_south_of_cpg_island"],
)
```

```{code-cell}
pd.testing.assert_frame_equal(
    full_cpg_island_anno_df_with_dist[cpg_island_grange_cols],
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals,
)
```

```{code-cell}
full_cpg_island_anno_df_with_dist["cpg_island_distance"].plot.hist(
    bins=np.linspace(-1e5, 1e5, 100)
)
```

## Add probe cytosine motif and strand

+++

### Download ref genomes

+++

ucsc mm10

```{code-cell}
if recompute:
    subprocess.run(["wget", "-O", mm10_fa, mm10_fa_url], check=True)
    subprocess.run(
        f"zcat {mm10_fa} | bgzip > {mm10_fa_bgz}", shell=True, check=True
    )
    subprocess.run(["samtools", "faidx", mm10_fa_bgz], check=True)
```

NCBI GRCm38

```{code-cell}
if recompute:
    subprocess.run(
        ["wget", "-O", ncbi_mm10_fa, ncbi_mm10_fa_url], check=True,
        capture_output=True,
    )
    subprocess.run(
        f"zcat {ncbi_mm10_fa} | bgzip > {ncbi_mm10_fa_bgz}",
        shell=True,
        check=True,
        capture_output=True,
    )
    subprocess.run(["samtools", "faidx", ncbi_mm10_fa_bgz], check=True, capture_output=True)
```

### download seq id mapper

```{code-cell}
if recompute:
    subprocess.run(
        [
            "wget",
            "-O",
            chrom_alias_txt_gz,
            chrom_alias_url,
        ],
        check=True,
        capture_output=True,
    )
```

```{code-cell}
chrom_aliases = pd.read_csv(  # type_ignore
    chrom_alias_txt_gz,
    sep="\t",
    header=None,
    names=["other_db_id", "ucsc_id", "db"],
)

ucsc_to_refseq_chrom_name = chrom_aliases.query('db == "refseq"').set_index(
    "ucsc_id"
)["other_db_id"]
```

### Get strands and motifs

+++

- a small number of probes does not cover cytosines
- comparing against two separately downloaded ref genomes just to be overly cautious about getting this right

```{code-cell}
ucsc_strand_ser, ucsc_motifs_ser = find_cytosine_strand_and_motif(
    temp_dir_name=temp_dir_name,
    df=illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.assign(
        Chromosome=lambda df: df["Chromosome"].astype(str).replace({"chrMT": "chrM"})
    ),
    fa_bgz=mm10_fa_bgz,
)
```

```{code-cell}
df = illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.assign(
    Chromosome=lambda df: df["Chromosome"]
    .astype(str)
    .replace({"chrMT": "chrM"})
    .map(ucsc_to_refseq_chrom_name)
)
assert df.Chromosome.notnull().all()

ncbi_strand_ser, ncbi_motifs_ser = find_cytosine_strand_and_motif(
    df=df,
    temp_dir_name=temp_dir_name,
    fa_bgz=ncbi_mm10_fa_bgz,
)
```

```{code-cell}
pd.testing.assert_series_equal(ncbi_strand_ser, ucsc_strand_ser)
pd.testing.assert_series_equal(ncbi_motifs_ser, ucsc_motifs_ser)
```

```{code-cell}
ncbi_motifs_ser.value_counts(dropna=False)
ncbi_strand_ser.value_counts(dropna=False)
ncbi_motifs_ser.loc[ncbi_strand_ser.eq("-")].value_counts()
```

```{code-cell}
motifs_df_final = (
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.reset_index(
        drop=True
    ).assign(strand=ncbi_strand_ser, motif=ncbi_motifs_ser)
)
```

```{code-cell}
motifs_df_final["motif"].value_counts(dropna=False)
```

## Merge all annotations

```{code-cell}
pd.testing.assert_frame_equal(
    full_cpg_island_anno_df_with_dist[["Chromosome", "Start", "End"]],
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.reset_index(
        drop=True
    ),
)
pd.testing.assert_frame_equal(
    merged_annos_new_chrom_dtype[["Chromosome", "Start", "End"]],
    illumina_probes_curated_chrom_defined_no_names_no_dup_intervals.reset_index(
        drop=True
    ),
)
```

```{code-cell}
exp_number_of_nas = (
    illumina_probe_intervals_bed_convention["Chromosome"].isin([np.nan, "0"]).sum()
)
```

```{code-cell}
m1 = pd.merge(
    illumina_probe_intervals_bed_convention,
    motifs_df_final,
    on=["Chromosome", "Start", "End"],
    how="left",
)
assert m1["motif"].isnull().sum() == exp_number_of_nas
```

```{code-cell}
m2 = pd.merge(
    m1,
    full_cpg_island_anno_df_with_dist,
    on=["Chromosome", "Start", "End"],
    how="left",
)
assert m2["cpg_island_region_name"].isnull().sum() == exp_number_of_nas
```

```{code-cell}
m3 = pd.merge(
    m2,
    merged_annos_new_chrom_dtype,
    on=["Chromosome", "Start", "End"],
    how="left",
)
assert m3["feat_class"].isnull().sum() == exp_number_of_nas
m3["score"] = "."
```

```{code-cell}
final_cols_d = {
    "Chromosome": "Chromosome",
    "Start": "Start",
    "End": "End",
    "name": "IlmnID",
    "score": "Score",
    "strand": "Strand",
    "motif": "Motif",
    "feat_class": "Genomic_region_class",
    "distance": "Distance_to_genomic_region",
    "gene_name": "Gene_name",
    "gene_id": "Gene_id",
    "transcript_id": "Transcript_id",
    "feat_strand": "Gene_strand",
    "cpg_island_region_name": "Cpg_island_region_class",
    "cpg_island_region_start": "Cpg_island_region_start",
    "cpg_island_region_end": "Cpg_island_region_end",
    "cpg_island_distance": "Cpg_island_distance",
}

final_anno_table = m3[list(final_cols_d.keys())].rename(columns=final_cols_d)
```

```{code-cell}
final_anno_table.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    probe_annos_one_row_bed_csv, sep="\t", header=True, index=False
)
final_anno_table.to_parquet(probe_annos_one_row_bed_parquet)
```
