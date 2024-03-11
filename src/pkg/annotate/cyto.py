import enum
import logging

import duckdb
import numpy as np
import pandas as pd


LOG = logging.getLogger(__name__)


class CytoCat(enum.Enum):
    NONE = 0
    HET = 1
    CEN = 2
    NOR = 3


# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
# https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.sql
CYTO_BAND_FORMAT = dict(
    sep="\t", header=None, names=["chrom", "start", "end", "name", "stain"]
)


def load_cyto_band(filename: str) -> pd.DataFrame:
    return pd.read_csv(filename, **CYTO_BAND_FORMAT)


def attach_cyto_category(
    nci_table: pd.DataFrame,
    band_table: pd.DataFrame,
    extend_nor: bool = False,
) -> pd.DataFrame:
    nci_cat_table = nci_table.reset_index()

    with duckdb.connect(":memory:") as con:
        con.execute(
            """
            SELECT b.stain
            FROM nci_cat_table n
            LEFT JOIN band_table b ON
              n.chrom = b.chrom and
              n.start >= b.start and
              n.end <= b.end
            ORDER BY n.index
            """
        )
        stains = [value for value, in con.fetchall()]

    stain_to_cat = {
        "gpos25": CytoCat.HET,
        "gpos33": CytoCat.HET,
        "gpos50": CytoCat.HET,
        "gpos66": CytoCat.HET,
        "gpos75": CytoCat.HET,
        "gpos100": CytoCat.HET,
        "acen": CytoCat.CEN,
        "stalk": CytoCat.NOR,
    }

    nci_cat_table.drop("index", axis=1, inplace=True)
    nci_cat_table.insert(
        len(nci_cat_table.columns),
        "cat",
        [stain_to_cat.get(stain, CytoCat.NONE) for stain in stains],
    )

    for chrom, track in nci_cat_table.groupby("chrom", sort=False):
        if (track["cat"] == CytoCat.CEN).sum() == 0:
            LOG.warn("No centromere was identified on %s", chrom)

    if (nci_cat_table["cat"] == CytoCat.NOR).sum() == 0:
        LOG.warn("No NOR was identified")

    if extend_nor:
        nci_cat_table = do_extend_nor(nci_cat_table)

    return nci_cat_table


def do_extend_nor(table: pd.DataFrame) -> pd.DataFrame:
    new_cats = []

    for chrom, track in table.groupby("chrom", sort=False):
        cats = track["cat"].values
        seen_nor = False
        arm_end = len(cats)

        for i, cat in enumerate(cats):
            if cat == CytoCat.NOR:
                seen_nor = True
            if cat == CytoCat.CEN:
                arm_end = i
                break

        if seen_nor:
            cats = cats.copy()
            cats[:arm_end] = CytoCat.NOR

        new_cats.append(cats)

    return table.assign(cat=np.concatenate(new_cats))
