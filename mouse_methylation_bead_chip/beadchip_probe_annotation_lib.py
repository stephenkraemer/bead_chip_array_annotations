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


