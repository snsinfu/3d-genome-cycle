import argparse
import json
import logging
import typing

import cooler
import numpy as np
import pandas as pd
import tqdm


from pkg.common.args import remove_none
from pkg.common.cli import invoke_main
from pkg.pc1.math import power_svd


LOG = logging.getLogger()
MATRIX_FORMAT = dict(dtype=np.float32, order="C")
OUTPUT_FORMAT = dict(sep="\t", float_format="%g", na_rep="nan", index=False)


def main(
    *,
    cool: str,
    output: str,
    aux_output: str | None = None,
    mask_intra: bool = False,
    use_covariance: bool = False,
    use_raw: bool = False,
):
    LOG.info("Opening cooler dataset %s", cool)
    clr = cooler.Cooler(cool)
    bins = clr.bins()[:]
    matrix = clr.matrix(balance=(not use_raw))

    LOG.info("Loading contact matrix of shape %s", matrix.shape)
    contact_matrix = np.asarray(matrix[:, :], **MATRIX_FORMAT)
    chrom_ranges = {chrom: clr.extent(chrom) for chrom in clr.chromnames}

    # Transform to O/E matrix in place.
    LOG.info("Computing O/E matrix")
    expected = estimate_expected_contacts(contact_matrix, chrom_ranges)
    data_matrix = contact_matrix
    del contact_matrix

    for patch in scan_chrom_rows(data_matrix, chrom_ranges):
        patch.trans_1[...] /= expected.inter
        patch.trans_2[...] /= expected.inter

        if mask_intra:
            # Leaving bad entries as-is
            patch.cis[np.isfinite(patch.cis)] = 1
        else:
            cis_size = len(patch.cis)
            for k in range(-cis_size + 1, cis_size):
                diag = np.diagonal(patch.cis, k)
                diag.setflags(write=True) # https://github.com/numpy/numpy/issues/5407
                diag[:] /= expected.intra[abs(k)]

    # Center columns for PCA
    LOG.info("Standardizing columns")
    coverages = np.nansum(data_matrix, axis=0)
    selection = coverages > 0
    data_matrix = data_matrix[:, selection]
    data_dim = data_matrix.shape[1]

    data_matrix[...] -= np.nanmean(data_matrix, axis=0)
    if not use_covariance:
        data_matrix[...] /= np.nanstd(data_matrix, axis=0)

    LOG.info(">> Found %d valid bins out of %d", data_dim, data_matrix.shape[0])

    # SVD
    LOG.info("Running SVD")
    progress = tqdm.tqdm(desc="SVD")
    delta_tol = 1e-4

    for svd in power_svd(data_matrix):
        progress.update()
        progress.set_postfix(delta=svd.delta)
        if svd.delta < delta_tol:
            break

    progress.close()

    #
    LOG.info("Computing PC1")
    pc1 = data_matrix @ svd.vector
    ev1 = unselect_vector(svd.vector, selection)

    # pc1 is multiplied by the first singular value. Rescale PC1 so that the
    # variance of pc1 gives the EVR, for convenience.
    data_var = np.nansum(np.nanvar(data_matrix, axis=0))
    pc1[:] /= np.sqrt(data_var)
    evr = np.nanvar(pc1)

    LOG.info("Explained variance ratio: %.1f %%", evr * 100)

    aux_data = {
        # The default JSON encoder requires conversion from np.float32 to float
        "explained_variance_ratio": float(evr),
        "cis_decay_profile": [float(x) for x in expected.intra],
        "trans_contact": expected.inter,
    }

    table = pd.DataFrame.from_dict({
        "chrom": bins["chrom"].values,
        "start": bins["start"].values,
        "end": bins["end"].values,
        "ev1": ev1,
        "pc1": pc1,
    })
    table.to_csv(output, **OUTPUT_FORMAT)

    if aux_output:
        with open(aux_output, "w") as file:
            json.dump(aux_data, file)


def unselect_vector(
    vector: np.ndarray,
    selection: np.ndarray,
    placeholder: float = np.nan,
) -> np.ndarray:
    result = np.full(len(selection), placeholder, dtype=vector.dtype)
    indices = np.arange(len(result))[selection]
    result[indices] = vector
    return result


class ContactPatch(typing.NamedTuple):
    chrom: str
    rows: np.ndarray
    trans_1: np.ndarray
    trans_2: np.ndarray
    cis: np.ndarray


def scan_chrom_rows(contact_matrix: np.ndarray, chrom_ranges: dict):
    for chrom, (start, end) in chrom_ranges.items():
        rows = contact_matrix[start:end]
        yield ContactPatch(
            chrom=chrom,
            rows=rows,
            trans_1=rows[:, :start],
            trans_2=rows[:, end:],
            cis=rows[:, start:end],
        )


class ExpectedContacts(typing.NamedTuple):
    intra: np.ndarray
    inter: float


def estimate_expected_contacts(
    contact_matrix: np.ndarray,
    chrom_ranges: dict,
) -> ExpectedContacts:
    max_separation = max(end - start for start, end in chrom_ranges.values())

    inter_sum = 0.0
    inter_count = 0
    intra_sums = np.zeros(max_separation, dtype=contact_matrix.dtype)
    intra_counts = np.zeros(max_separation, dtype=contact_matrix.dtype)

    for patch in scan_chrom_rows(contact_matrix, chrom_ranges):
        sum_1, count_1 = valid_sum(patch.trans_1)
        sum_2, count_2 = valid_sum(patch.trans_2)
        inter_sum += sum_1 + sum_2
        inter_count += count_1 + count_2

        for s in range(len(patch.cis)):
            diag = np.diag(patch.cis, s)
            d_sum, d_count = valid_sum(diag)
            intra_sums[s] += d_sum
            intra_counts[s] += d_count

    with np.errstate(invalid="ignore", divide="ignore"):
        return ExpectedContacts(
            intra=(intra_sums / intra_counts),
            inter=float(inter_sum / inter_count),
        )


def valid_sum(vec: np.ndarray) -> (float, int):
    valid = np.isfinite(vec)
    return vec[valid].sum(), valid.sum()


def parse_args() -> dict:
    parser = argparse.ArgumentParser(
        description="Compute compartment signal",
    )
    arg = parser.add_argument

    arg(
        "--svd-tolerance",
        metavar="1e-4",
        type=float,
        default=None,
        help="Convergence threshold used in SVD calculation",
    )
    arg(
        "--use-raw",
        action="store_true",
        default=False,
        help="Use unnormalized contact matrix",
    )
    arg(
        "--use-covariance",
        action="store_true",
        default=False,
        help="Diagonalize O/E covariance, not correlation, matrix",
    )
    arg(
        "--mask-intra",
        action="store_true",
        default=False,
        help="Mask intra-chromosomal O/E contacts by setting 1",
    )
    arg(
        "--aux-output",
        metavar="aux.json",
        type=str,
        default=None,
        help="JSON file to store information output on PCA",
    )
    arg(
        "--output",
        metavar="pc1.tsv",
        type=str,
        help="TSV file to store PC1 values in",
    )
    arg(
        "cool",
        type=str,
        help="Path to a cool dataset",
    )

    return remove_none(vars(parser.parse_args()))


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
