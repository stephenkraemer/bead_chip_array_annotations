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
print("reloaded")

# # Imports

# +
from typing import Dict

import numpy as np
import pandas as pd
import pyranges as pr


# -

# # Annotation mering

# - for intergenic, everything will be NaN anyway


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
        # # for testing
        # .query('feat_class == "Promoter"')
        # .head(50)
        # # / for testing
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

    # noinspection PyUnresolvedReferences
    assert (
        full_merged_annos["gtfanno_uid"] == np.arange(full_merged_annos.shape[0])
    ).all()

    return full_merged_annos.reset_index(drop=True)


def _agg_multi_value_columns(ser):
    if ser.eq("nan").all():
        return "nan"
    else:
        return ser.str.cat(sep=",")


# # CpG island classification

# Currently the classification picks the region class overlapping the largest fraction of the DMR.
#
# Given that CpG island hits are so rare, these hits may be so interesting that we want to be more sensitive.  Possibile alternative classification strategy:
#
# - filter overlaps, pick from
#   - \> p% overlap with region class
#   - \> q% overlap with DMR
#   - \> #bp overlap
# - among the remaining overlaps, pick according to precedence order islands > shores > shelves > open sea

# +
def classify_cpg_island_overlap(
    granges_df: pd.DataFrame, cpg_islands_pickle_d: Dict[str, str]
) -> pd.DataFrame:
    """Classify DMRs as cpg island, shelve, shore or open sea

    Parameters
    ----------
    granges_df
        region_id // ['Chromosome', 'Start', 'End'] ... optional other columns
    cpg_islands_pickle_d
        'CpG island', 'CpG shelve', ... -> df pickle
        df with no header, columns Chromosome Start End (?)

    Classification strategy and alternatives
    -----------------------------------------
    currently the classification picks the region class overlapping the largest fraction of the DMR
    given that CpG island hits are so rare, these hits may be so interesting that we want to be more sensitive
    Possibile alternative strategy
    - filter overlaps, pick from
      - > p% overlap with region class
      - > q% overlap with DMR
      - > #bp overlap
    - among the remaining overlaps, pick according to precedence order islands > shores > shelves > open sea

    Returns
    -------
    pd.DataFrame
        Chromosome Start End gain_loss_cassif
    """

    granges_gr = pr.PyRanges(granges_df)
    # need to consider double probes at same gcoordinate, different name
    grange_cols = ["Chromosome", "Start", "End", "name"]

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
        join_df = join_gr.df[["Chromosome", "Start", "End", "name", "region_name"]]
        dfs.append(join_df)

    cpg_island_anno_df = pd.concat(dfs, axis=0).reset_index(drop=True)

    # there are some cpgs right at the boundary
    cpg_island_anno_df.groupby(grange_cols, observed=True)[
        "region_name"
    ].nunique().value_counts()
    # there are also some CpGs overlapped by multiple CpG shores for example
    cpg_island_anno_df.groupby(grange_cols, observed=True).size().value_counts()
    cpg_island_anno_df.set_index(grange_cols).loc[
        cpg_island_anno_df.groupby(grange_cols, observed=True).size().gt(2)
    ].sort_values(grange_cols)
    # However, annotation to a CpG island is alway unique
    cpg_island_anno_df.query('region_name == "CpG islands"').groupby(
        grange_cols, observed=True
    ).size().value_counts()
    # so there are no example regions for multiple cpg islands
    cpg_island_anno_df.query('region_name == "CpG islands"').set_index(grange_cols).loc[
        cpg_island_anno_df.query('region_name == "CpG islands"')
        .groupby(grange_cols, observed=True)
        .size()
        .gt(1)
    ]

    # pick unique classif according to precedence
    cpg_island_anno_df_unique = (
        cpg_island_anno_df.assign(
            region_name=lambda df: pd.Categorical(
                df["region_name"],
                categories=["CpG islands", "CpG shelves", "CpG shores", "open sea"],
            )
        )
        .sort_values(
            ["Chromosome", "Start", "End", "name", "region_name"], ascending=True
        )
        .groupby(grange_cols, observed=True, as_index=False)
        .first()
        .assign(
            Chromosome=lambda df: df["Chromosome"].astype(granges_df.Chromosome.dtype)
        )
    )

    # detail all overlapping cpg island features
    (cpg_island_anno_df.groupby(grange_cols))

    # reindex to include open sea cpgs
    full_cpg_island_anno_df = pd.merge(
        granges_df,
        cpg_island_anno_df_unique,
        on=grange_cols,
        how="left",
    )

    assert (
        full_cpg_island_anno_df.region_name.notnull().sum()
        == cpg_island_anno_df_unique.shape[0]
    )

    full_cpg_island_anno_df["region_name"] = full_cpg_island_anno_df[
        "region_name"
    ].fillna("open sea")

    # pyranges may change sorting order, but I think pyranges using same sorting order as probes df,
    # so should be fine
    pd.testing.assert_frame_equal(
        full_cpg_island_anno_df[grange_cols],
        granges_df,
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
    nearest_cpg_island_df.loc[
        is_upstream, "distance_signed"
    ] = nearest_cpg_island_df.loc[is_upstream].eval("(End - 1) - Start_b")
    nearest_cpg_island_df.loc[
        is_downstream, "distance_signed"
    ] = nearest_cpg_island_df.loc[is_downstream].eval("Start - (End_b - 1)")

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
        full_cpg_island_anno_df_with_dist.Chromosome.astype(granges_df.Chromosome.dtype)
    )
    full_cpg_island_anno_df_with_dist["north_south_of_island"] = np.sign(
        full_cpg_island_anno_df_with_dist["distance_signed"]
    ).replace({0: "", 1: "N", -1: "S"})

    full_cpg_island_anno_df_with_dist.loc[
        ~full_cpg_island_anno_df_with_dist.region_name.isin(
            ["CpG shores", "CpG shelves"]
        ),
        "north_south_of_island",
    ] = ""

    # check that cpg islands and open sea have no north / south, but shelves and shores do have it
    # check that strands are approx 50/50
    pd.crosstab(
        full_cpg_island_anno_df_with_dist["region_name"],
        full_cpg_island_anno_df_with_dist["north_south_of_island"],
    )

    pd.testing.assert_frame_equal(
        full_cpg_island_anno_df_with_dist[grange_cols], granges_df
    )

    return full_cpg_island_anno_df_with_dist


def annotate_with_motif_and_strand(granges_df, chrom_dtype):
    """

    Parameters
    ----------
    granges_df
        BED coordinates

    Returns
    -------

    """
    grange_cols = ["Chromosome", "Start", "End", "motif", "strand"]
    annotated_chrom_dfs_d = {}
    for chrom, chrom_df in granges_df.groupby("Chromosome"):
        chrom_no_prefix = chrom.replace("chr", "")
        chrom_index_bed_gz = f"/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/genomes/GRCm38mm10_PhiX_Lambda/GRCm38mm10_PhiX_Lambda_CG-CHG-CHH_{chrom}.bed.gz"
        chrom_index_df = pd.read_csv(
            chrom_index_bed_gz,
            sep="\t",
            usecols=[0, 1, 2, 3, 5],
            comment="#",
            header=None,
            names=grange_cols,
        )
        chrom_df_annotated = pd.merge(
            chrom_df,
            chrom_index_df,
            on=["Chromosome", "Start", "End"],
            how="left",
        )
        pd.testing.assert_frame_equal(
            chrom_df[grange_cols], chrom_df_annotated[grange_cols]
        )
        annotated_chrom_dfs_d[chrom] = chrom_df_annotated
    res = (
        pd.concat(annotated_chrom_dfs_d)
        .assign(Chromosome=lambda df: df["Chromosome"].astype(chrom_dtype))
        .sort_values(grange_cols)
        .reset_index(drop=True)
    )
    return res

