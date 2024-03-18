import argparse
import logging

import cooler
import duckdb
import pandas as pd
import tqdm

from pkg.common.args import remove_none
from pkg.common.cli import invoke_main


LOG = logging.getLogger()
VIRTUAL_CHAINS = {"nucleoli"}
CHUNK_SIZE = 5_000_000


def main(
    *,
    output: str,
    input: str,
    no_balancing: bool = False,
):
    LOG.info("Opening cool dataset: %s", input)
    clr = cooler.Cooler(input)

    # Infer chromosome names from chain names.
    chrom_copies = infer_chromosome_copies(
        [name for name in clr.chromnames if name not in VIRTUAL_CHAINS]
    )
    LOG.info("Inferred chromosomes: %s", ", ".join(chrom_copies))

    # Define output bins.
    input_bins = clr.bins()[:]
    output_bins = make_output_bins(input_bins, chrom_copies)
    chrom_chain_mapping = make_mapping(input_bins, output_bins, chrom_copies)
    LOG.info("Bins reduced: %d -> %d", len(input_bins), len(output_bins))

    # Create a cool dataset from incrementally de-phased input contacts.
    cooler.create_cooler(
        output,
        output_bins,
        dephase_pixels(clr, chrom_chain_mapping, CHUNK_SIZE),
        assembly=clr.info.get("genome-assembly"),
    )

    output_clr = cooler.Cooler(output)
    LOG.info(">> Stored %d pixels", len(output_clr.pixels()))

    # Optional balancing. This may also take a while.
    if not no_balancing:
        LOG.info("Balancing contact matrix")
        cooler.balance_cooler(output_clr, store=True)


def infer_chromosome_copies(chain_names: list[str]) -> dict[str, list[str]]:
    chrom_copies: dict[str, list[str]] = {}

    for chain_name in chain_names:
        match chain_name.split(":"):
            case [chrom, suffix]:
                if chrom not in chrom_copies:
                    chrom_copies[chrom] = []
                chrom_copies[chrom].append(suffix)
            case _:
                LOG.warn("Skipping unrecognized chain: %s", chain_name)

    return chrom_copies


def make_output_bins(
    input_bins: pd.DataFrame,
    chrom_copies: dict[str, list[str]],
) -> pd.DataFrame:
    canon_renaming = {
        f"{chrom}:{suffixes[0]}": chrom
        for chrom, suffixes in chrom_copies.items()
    }
    output_bins = (
        input_bins
        .query("chrom in @canon_renaming")
        .reset_index(drop=True)
        .replace(canon_renaming)
    )
    return output_bins


def make_mapping(
    input_bins: pd.DataFrame,
    output_bins: pd.DataFrame,
    chrom_copies: dict[str, list[str]],
) -> pd.DataFrame:
    chrom_chain_mapping: list[dict] = []

    def chrom_range(bins: pd.DataFrame, name: str) -> tuple[int, int]:
        indices = bins.query("chrom == @name").index
        return indices[0], indices[-1] + 1

    for chrom, suffixes in chrom_copies.items():
        chrom_start, chrom_end = chrom_range(output_bins, chrom)

        for suffix in suffixes:
            chain = f"{chrom}:{suffix}"
            chain_start, chain_end = chrom_range(input_bins, chain)

            chrom_chain_mapping.append({
                "chrom_start": chrom_start,
                "chrom_end": chrom_end,
                "chain_start": chain_start,
                "chain_end": chain_end,
            })

    return pd.DataFrame.from_records(chrom_chain_mapping)


def dephase_pixels(
    clr: cooler.Cooler,
    chrom_chain_mapping: pd.DataFrame,
    chunk_size: int,
):
    source_pixels = clr.pixels()
    n_pixels = len(source_pixels)

    LOG.info("Dephasing %d pixels", n_pixels)

    for chunk_start in tqdm.trange(0, n_pixels, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_pixels)
        pixels = source_pixels[chunk_start:chunk_end]
        pixels = _map_pixels(pixels, chrom_chain_mapping)
        pixels = _dedupe_pixels(pixels)
        yield pixels

    # After consuming pixels Cooler merges them into a single cool.
    LOG.info("Merging pixels...")


def _map_pixels(
    pixels: pd.DataFrame,
    chrom_chain_mapping: pd.DataFrame,
) -> pd.DataFrame:
    # Map diploid pixel coordinates to haploid ones, using given mapping table.
    duckdb.execute(
        """
        SELECT
          (bin1_id - m1.chain_start + m1.chrom_start) as bin1_id,
          (bin2_id - m2.chain_start + m2.chrom_start) as bin2_id,
          count
        FROM pixels
        JOIN chrom_chain_mapping m1 ON
          bin1_id >= m1.chain_start AND bin1_id < m1.chain_end
        JOIN chrom_chain_mapping m2 ON
          bin2_id >= m2.chain_start AND bin2_id < m2.chain_end
        """
    )
    return duckdb.df()


def _dedupe_pixels(pixels: pd.DataFrame) -> pd.DataFrame:
    # Merge superposed pixels and output the results in the upper-triangular
    # format that Cooler expects.
    duckdb.execute(
        """
        WITH triu_pixels AS (
            SELECT
                 least(bin1_id, bin2_id) as bin1_id,
              greatest(bin1_id, bin2_id) as bin2_id,
              count
            FROM pixels
        )
        SELECT bin1_id, bin2_id, sum(count)::INTEGER as count
        FROM triu_pixels
        GROUP BY bin1_id, bin2_id
        ORDER BY bin1_id, bin2_id
        """
    )
    return duckdb.df()


def parse_args() -> dict:
    parser = argparse.ArgumentParser(
        description="Aggregate homologous contacts",
    )
    arg = parser.add_argument

    arg(
        "--no-balancing",
        action="store_true",
        help="Do not balance the resulting matrix",
    )
    arg(
        "--output",
        metavar="out.cool",
        type=str,
        help="Cool dataset to write the de-phased contact matrix to",
    )
    arg(
        "input",
        metavar="in.cool",
        type=str,
        help="Contact map created with the cool script",
    )

    return remove_none(vars(parser.parse_args()))


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
