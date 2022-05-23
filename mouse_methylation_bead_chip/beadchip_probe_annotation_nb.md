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

+++ {"tags": [], "heading_collapsed": "true"}


# TODO

+++

- check sorting order of chromosomes throughout analysis
- there is a bug in gtfanno which makes it fail to read the appris principal score, does this influence the results?
- there is a "transcript = 0" count in the gtfanno basic stat output. is this indicative of a problem?
- double check the changes made to allow DCRD intervals enveloping the promoter interval
- check annos in IGV
- are feature coordinates really 0-based, right open?

+++

annos to add
- motif: CG, CHH, CHG
- strand, illumina strands are not cytosine strands
- hematopoietic regions
    - cis reg atlas
    - vision
    - amit enhancers
- general regulatory regions
    - ensembl reg regions
    - chrom hmm - ask maxi again what he had in mind here - forgot which resource he mentioned
- tfbs

+++ {"tags": []}


# Setup

+++ {"tags": [], "heading_collapsed": "true"}


## Resource parameters

```{code-cell}
n_cores = 12
```

+++ {"tags": [], "heading_collapsed": "true"}


## Imports

```{code-cell}
:tags: []

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
```

```{code-cell}
%matplotlib inline
```

```{code-cell}
import mouse_methylation_bead_chip.beadchip_probe_annotation_lib as lib
import mouse_methylation_bead_chip.beadchip_probe_annotation_paths as paths
```

+++ {"tags": [], "heading_collapsed": "true"}


## Rerun flags

```{code-cell}
recompute = True
```

+++ {"tags": [], "heading_collapsed": "true"}


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

+++ {"tags": []}


# Paths

temp_dir_obj = tempfile.TemporaryDirectory(dir=paths.project_dir)
temp_dir_name = temp_dir_obj.name
temp_dir_name

+++ {"tags": []}


# Analysis

+++ {"tags": [], "heading_collapsed": "true"}


## Prepare input data

+++ {"tags": [], "heading_collapsed": "true"}


### CpG island annos

```{code-cell}
import mouse_hema_meth.genome_annotations.get_genome_annos_paths as get_genome_annos_paths

cpg_islands_pickle_d = get_genome_annos_paths.cpg_islands_shores_shelves_pickle_paths_d
```

+++ {"tags": [], "heading_collapsed": "true"}


### Prepare gene annotation

+++ {"tags": [], "heading_collapsed": "true"}


#### download gencode

```{code-cell}
:tags: []

if recompute:
    subprocess.run(
        ["wget", "-O", paths.gencode_gtf, paths.gencode_download_url],
        check=True,
    )
```

```{code-cell}
:tags: []

!zcat {paths.gencode_gtf} | head -n 6
```

+++ {"tags": [], "heading_collapsed": "true"}


#### Filter and reformat gencode GTF

+++

- restrict to canonical transcripts
- restrict to coding transcripts
- remove chr prefix
- change M to MT

```{code-cell}
gencode_df = pr.read_gtf(paths.gencode_gtf, as_df=True, duplicate_attr=True)
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
gencode_pr.to_gtf(paths.gencode_coding_canonical_gtf)
```

```{code-cell}
!zcat {paths.gencode_coding_canonical_gtf} | head
```

verify gtf

```{code-cell}
---
jupyter:
  outputs_hidden: true
tags: []
---
!zcat {paths.gencode_coding_canonical_gtf} | grep ^protein_coding
```

```{code-cell}
---
jupyter:
  outputs_hidden: true
tags: []
---
!zcat {paths.gencode_coding_canonical_gtf} | grep ^appris
```

+++ {"tags": []}


### Prepare and inspect probes files

+++ {"tags": [], "heading_collapsed": "true"}


#### Probe file from Maxi

+++ {"tags": [], "heading_collapsed": "true"}


##### Inspect original probes file

+++

- file has duplicates
- file is not fully sorted

+++ {"tags": []}


###### General overview

```{code-cell}
!head -n 3 {paths.original_probes_bed}
```

```{code-cell}
!cut -f 1 < {paths.original_probes_bed} | uniq
```

```{code-cell}
original_probes_df = pd.read_csv(
    paths.original_probes_bed,
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
```

+++ {"tags": []}


###### File is not fully sorted

+++

**Note that the original probes df is not completely sorted on Start/End**

```{code-cell}
original_probes_df_sorted = original_probes_df.sort_values(
    ["Chromosome", "Start", "End"]
).reset_index(drop=True)
original_probes_df_sorted
```

```{code-cell}
original_probes_df_sorted.Chromosome.dtype
```

+++ {"tags": [], "heading_collapsed": "true"}


###### Several probes are present with the same coordinates, but different names

```{code-cell}
original_probes_df.loc[
    original_probes_df.duplicated(["Chromosome", "Start", "End"], keep=False)
]
```

```{code-cell}
original_probes_df.loc[
    original_probes_df.duplicated(["Chromosome", "Start", "End", "name"], keep=False)
]
```

+++ {"tags": [], "heading_collapsed": "true"}


##### Reformat probes file

+++

- need to resort
- need to remove chr prefix
- drop duplicates

```{code-cell}
probes_df_no_prefix_sorted = (
    original_probes_df.assign(
        Chromosome=lambda df: df["Chromosome"].str.replace("chr", ""),
    )[["Chromosome", "Start", "End"]]
    .drop_duplicates()
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)
```

```{code-cell}
probes_df_no_prefix_sorted.to_csv(
    paths.reformatted_probes_bed, sep="\t", header=False, index=False
)
```

```{code-cell}
!head {paths.reformatted_probes_bed}
```

+++ {"tags": [], "heading_collapsed": "true"}


#### Illumina probe file

+++ {"tags": []}


##### Schema

+++

- MFG_CHANGE probes haben ein problem
- there may be one row separating assay probes from controls somewhere in the dataframe? (info from Maxi)

+++ {"tags": [], "heading_collapsed": "true"}


##### Download

```{code-cell}
if recompute:
    subprocess.run(["wget", "-O", paths.illumina_probes_csv, paths.illumina_probes_url], check=True)
```

```{code-cell}
!head {paths.illumina_probes_csv}
```

+++ {"tags": [], "heading_collapsed": "true"}


##### Get curated BED intervals for probes

```{code-cell}
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
```

Fields, drop fields with longish sequence strings for display

```{code-cell}
:tags: []

illumina_probes.drop(["Forward_Sequence", "Top_Sequence"], axis=1).iloc[0].to_frame()
```

There are nan chromosomes entries, and also some entries for chromosome 0, just 410, so I assume this can just be discarded as controls or something similar

```{code-cell}
illumina_probes.CHR.value_counts()
```

checked manually in my index files: 1-based Start info is in MAPINFO

+++

- for comparability with Maxis probes, also add 'chr' prefix and make Categorical
- provide BED interval for cytosine

```{code-cell}
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
```

drop duplicate rows, remove prefix, change to alphabetic sorting order

```{code-cell}
illumina_probes_curated_chrom_defined.assign(
    Chromosome=lambda df: df.Chromosome.astype(str).str.replace("chr", "")
).iloc[:, 0:3].sort_values(["Chromosome", "Start", "End"]).drop_duplicates().to_csv(
    paths.illumina_coordinate_bed, sep="\t", header=False, index=False
)
```

```{code-cell}
!head {paths.illumina_coordinate_bed}
```

+++ {"tags": [], "heading_collapsed": "true"}


##### Check against Maxis probes to see whether I have correct manifest file

+++

this is the correct manifest file - maxis coordinates are shifted when on minus strand

```{code-cell}
pd.merge(
    original_probes_df_sorted,
    illumina_probes_curated_chrom_defined,
    on=["Chromosome", "Start", "End", "name"],
    how="inner",
)
```

```{code-cell}
df = pd.merge(
    original_probes_df_sorted,
    illumina_probes_curated_chrom_defined,
    on=["name"],
    how="inner",
)
display(df)
assert df.shape[0] == original_probes_df_sorted.shape[0]
```

+++ {"tags": [], "heading_collapsed": "true"}


##### Add motif and strand

+++ {"tags": []}


## Annotation

+++ {"tags": []}


### Gene annotation

+++ {"tags": [], "heading_collapsed": "true"}


#### Perform annotation

```{code-cell}
:tags: []

%%time
ga.annotate(
    query_bed=paths.illumina_coordinate_bed,
    gtf_fp=paths.gencode_coding_canonical_gtf,
    trunk_path=paths.custom_intervals_trunk_path,
    tmpdir=temp_dir_name,
    promoter=(-1500, 500),
    distant_cis_regulatory_domain=(-100_000, 100_000),
)
```

+++ {"tags": []}


#### Inspect annotations

```{code-cell}
primary_annos = pd.read_pickle(paths.custom_intervals_results_paths_d["primary_annos_p"])
```

```{code-cell}
primary_annos.shape
```

+++ {"tags": []}


##### General checks

```{code-cell}
primary_annos.query('feat_class == "Promoter"').head(3)
```

```{code-cell}
primary_annos.query('feat_class == "exon"').head(3)
```

+++ {"tags": []}


##### Multiple assignments per region

+++ {"tags": []}


###### How is this distributed across feature classes?

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

```{raw-cell}
multi_annos_crosstab.to_clipboard()
```

+++ {"tags": []}


###### Example for Promoter multiple annotations - random samples indicate that these are indeed ambiguous sites

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
http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000031349;r=X:71962163-72002120

Rpl7
http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000043716;r=1:16171519-16174886

+++ {"tags": [], "heading_collapsed": "true"}


#### merge annotations

+++

Merging strategy: keep all
- for Promoters, the window is relatively small. Ranking on TSS distance in such a small window seems arbitrary.
- for enhancer candidates, a typical strategy would be to identify all TSS in +-100 kb window and try to find the target through correlation with gene expression, eg PMID: 30686579. So it also makes sense to indicate all genes in the window to give an impression of the number of possible target genes.

```{code-cell}
%%time
merged_annos = lib.merge_annos(primary_annos=primary_annos)
```

```{code-cell}
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
```

```{code-cell}
merged_annos_final = pd.merge(
    merged_annos_new_chrom_dtype,
    illumina_probes_curated_chrom_defined,
    on=["Chromosome", "Start", "End"],
    how="left",
)
```

```{code-cell}
merged_annos_final.head(3)
```

```{code-cell}
merged_annos_final.shape
```

```{code-cell}
illumina_probes_curated_chrom_defined
```

```{code-cell}
assert merged_annos_final["name"].notnull().all()
```

```{code-cell}
pd.testing.assert_frame_equal(
    merged_annos_final[["Chromosome", "Start", "End"]],
    illumina_probes_curated_chrom_defined[["Chromosome", "Start", "End"]].astype(
        {"Start": "i8", "End": "i8"}
    ),
)
```

```{code-cell}
merged_annos_final.iloc[0]
```

+++ {"heading_collapsed": "true", "tags": []}


#### Finalize annotation tables

```{code-cell}
merged_annos_final.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    paths.gene_annos_primary_one_row, sep="\t", header=True, index=False
)
```

```{code-cell}
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
```

```{code-cell}
primary_annos_final.rename(columns={"Chromosome": "#Chromosome"}).to_csv(
    paths.gene_annos_primary_multi_row, sep="\t", header=True, index=False
)
```

```{code-cell}
!head {paths.gene_annos_primary_multi_row}
```

+++ {"tags": [], "heading_collapsed": "true"}


### CpG island annotations

```{code-cell}
cpg_island_classif_df = lib.classify_cpg_island_overlap(
    granges_df=original_probes_df_sorted,
    cpg_islands_pickle_d=cpg_islands_pickle_d,
)
cpg_island_classif_df.head(3)
```

### Merge all annotations

```{code-cell}
cpg_island_classif_df
merged_annos_final
```

# End

```{code-cell}

```