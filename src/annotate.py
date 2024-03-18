import argparse
import dataclasses
import enum
import logging
import os
import signal

import numpy as np
import pandas as pd

from pkg.common.args import remove_none
from pkg.common.cli import invoke_main
from pkg.annotate.cyto import CytoCat, load_cyto_band, attach_cyto_category


LOG = logging.getLogger()
NCI_FORMAT = dict(sep="\t")


def main(
    *,
    tristate: float = 0,
    activate_nor: list[str] = [],
    extend_nor: bool = False,
    smooth_window: int = 10,
    output_filename: str,
    nci_filename: str,
    band_filename: str,
):
    nci_table = pd.read_csv(nci_filename, **NCI_FORMAT)
    band_table = load_cyto_band(band_filename)
    nci_cat_table = attach_cyto_category(nci_table, band_table, extend_nor=extend_nor)

    chrom_names = list(nci_cat_table["chrom"].unique())
    chains = design_diploid_chains(chrom_names, activate_nor)

    exclude = ["chrX"]
    basic_center, basic_scale = compute_normalizer(
        nci_cat_table.query("chrom not in @exclude")["score"].values
    )

    output = open(output_filename, "w")
    need_header = True

    for chain in chains:
        chrom = chain.chrom
        nci_cat_track = nci_cat_table.query("chrom == @chrom")
        chain_length = len(nci_cat_track)

        cats = nci_cat_track["cat"].values
        tags = [[] for _ in range(chain_length)]

        for i in range(chain_length):
            match cats[i]:
                case CytoCat.NOR:
                    tags[i].append("anor" if chain.activate_nor else "bnor")
                case CytoCat.CEN:
                    tags[i].append("cen")
                case CytoCat.HET:
                    tags[i].append("het")

        scores = (
            nci_cat_track["score"]
            .rolling(window=smooth_window, center=True, min_periods=1)
            .mean()
            .values
        )
        parameters = np.zeros((chain_length, 2))

        match chain.annot_scheme:
            case AnnotScheme.GENOME_WIDE:
                center, scale = basic_center, basic_scale
            case AnnotScheme.SINGLE_CHROM:
                center, scale = compute_normalizer(scores)
            case AnnotScheme.ALL_B:
                center, scale = np.inf, 1

        for i in range(chain_length):
            z_score = (scores[i] - center) / scale
            chrom_type = determine_chromatin_type(z_score, tristate, tags[i])
            # chrom_type = determine_chromatin_type(z_score, tristate, tags[i])

            tags[i].append(CHROM_TYPE_TAGS[chrom_type])
            parameters[i] = CHROM_TYPE_PARAMETERS[chrom_type]
            #tags[i].append(CHROM_TYPE_TAGS[chrom_type])
            parameters[i] = CHROM_TYPE_PARAMETERS[determine_chromatin_type(z_score, tristate)]

            # Remove "het" as it is not used and also might be confusing.
            if "het" in tags[i]:
                tags[i].remove("het")

        track = pd.DataFrame({
            "chain": chain.name,
            "start": nci_cat_track["start"].values,
            "end": nci_cat_track["end"].values,
            "A": parameters[:, 0],
            "B": parameters[:, 1],
            "tags": map(",".join, tags),
        })

        track.to_csv(output, sep="\t", float_format="%g", index=False, header=need_header)
        need_header = False
        output.flush()


class AnnotScheme(enum.Enum):
    GENOME_WIDE = 0
    SINGLE_CHROM = 1
    ALL_B = 2


@dataclasses.dataclass
class Chain:
    name: str
    chrom: str
    annot_scheme: int
    activate_nor: bool = False


class ChromType(enum.Enum):
    A = 1
    B = 2
    U = 3


CHROM_TYPE_HEURISTICS = {
    "cen": ChromType.B,
    "anor": ChromType.A,
    "bnor": ChromType.B,
}

CHROM_TYPE_TAGS = {
    ChromType.A: "A",
    ChromType.B: "B",
    ChromType.U: "u",
}

CHROM_TYPE_PARAMETERS = {
    ChromType.A: (1.0, 0.0),
    ChromType.B: (0.0, 1.0),
    ChromType.U: (0.5, 0.5),
}


def design_diploid_chains(
    chrom_names: list[str],
    active_nor_patterns: list[str],
):
    chains = []

    normal_chroms = chrom_names.copy()
    use_xa = use_xb = use_y = False

    if "chrX" in chrom_names:
        use_xa = use_xb = True
        normal_chroms.remove("chrX")

    if "chrY" in chrom_names:
        use_xb = False
        use_y = True
        normal_chroms.remove("chrY")

    for chrom in normal_chroms:
        chain_name = chrom + ":a"

        chain = Chain(
            name=chain_name,
            chrom=chrom,
            annot_scheme=AnnotScheme.GENOME_WIDE,
            activate_nor=any(chain_name.endswith(pat) for pat in active_nor_patterns),
        )
        chains.append(chain)

    if use_xa:
        chain = Chain(
            name="chrX:a",
            chrom="chrX",
            annot_scheme=AnnotScheme.SINGLE_CHROM,
        )
        chains.append(chain)

    for chrom in normal_chroms:
        chain_name = chrom + ":b"

        chain = Chain(
            name=chain_name,
            chrom=chrom,
            annot_scheme=AnnotScheme.GENOME_WIDE,
            activate_nor=any(chain_name.endswith(pat) for pat in active_nor_patterns),
        )
        chains.append(chain)

    if use_xb:
        chain = Chain(
            name="chrX:b",
            chrom="chrX",
            annot_scheme=AnnotScheme.ALL_B,
        )
        chains.append(chain)

    if use_y:
        chain = Chain(
            name="chrY:b",
            chrom="chrY",
            annot_scheme=AnnotScheme.GENOME_WIDE,
        )
        chains.append(chain)

    return chains


def determine_chromatin_type(z_score: float, tristate: float, tags: list[str]) -> int:
    if np.isnan(z_score):
        return infer_chromatin_type(tags)

    if z_score > tristate:
        return ChromType.A
    if z_score < -tristate:
        return ChromType.B
    return ChromType.U


def infer_chromatin_type(tags: list[str]) -> int:
    for key, typ in CHROM_TYPE_HEURISTICS.items():
        if key in tags:
            return typ
    return ChromType.U


def compute_normalizer(values: np.ndarray) -> tuple[float, float]:
    MAD_FACTOR = 1.4826
    center = np.nanmedian(values)
    scale = np.nanmedian(np.abs(values - center)) * MAD_FACTOR
    return center, scale


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tristate", type=float)
    parser.add_argument("--activate-nor", type=str)
    parser.add_argument("--extend-nor", action="store_true")
    parser.add_argument("--smooth-window", type=int)
    parser.add_argument("--nci", dest="nci_filename", required=True)
    parser.add_argument("--band", dest="band_filename", required=True)
    parser.add_argument("--output", dest="output_filename", required=True)

    args = vars(parser.parse_args())
    reparse_list(args, "activate_nor")
    return remove_none(args)


def reparse_list(args: dict, key: str):
    if arg := args.get(key):
        args[key] = arg.split(",")


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
