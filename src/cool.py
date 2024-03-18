import argparse
import logging
import io
import json

import cooler
import h5py
import numpy as np
import pandas as pd
import tqdm


from pkg.common.args import remove_none
from pkg.common.cli import invoke_main


LOG = logging.getLogger()
NUCLEOLAR_CHAIN = "nucleoli"


def main(
    *,
    output: str,
    frames: list[slice] | None = None,
    input_sims: list[str],
):
    with h5py.File(input_sims[0], "r") as store:
        # ?
        assembly = None

        # Load chain configurations
        metadata = store["metadata"]

        n_bins = len(metadata["particle_types"])
        chain_names = [name.decode() for name in metadata["chain_names"]]
        chain_ranges = metadata["chain_ranges"][:]

        # Copy chromosomal bins from the simulation input (chains.tsv).
        chain_bins = pd.read_csv(
            io.StringIO(metadata["chains_source"][()].decode()),
            sep="\t",
        )

    bin_start_coords = chain_bins["start"].values
    bin_end_coords = chain_bins["end"].values
    bin_size = np.max(bin_end_coords - bin_start_coords)

    # Define bins. Homologous chromosomes are treated as distinct chains (as in
    # the chains.tsv simulation input), and the nucleolar particles are treated
    # as bins on a virtual nucleolar chain, giving NAD contacts. The NAD
    # contacts are trimmed out in the dephase.py script.
    bins_chrom = np.empty(n_bins, dtype=object)
    bins_start = np.empty(n_bins, dtype=int)
    bins_end = np.empty(n_bins, dtype=int)
    chains_end = chain_ranges.max()

    for name, (start, end) in zip(chain_names, chain_ranges):
        indices = np.arange(end - start)
        bins_chrom[start:end] = name
        bins_start[start:end] = bin_start_coords[start:end]
        bins_end[start:end] = bin_end_coords[start:end]

    indices = np.arange(n_bins - chains_end)
    bins_chrom[chains_end:] = NUCLEOLAR_CHAIN
    bins_start[chains_end:] = indices * bin_size
    bins_end[chains_end:] = (indices + 1) * bin_size

    bins = pd.DataFrame.from_dict({
        "chrom": bins_chrom,
        "start": bins_start,
        "end": bins_end,
    })

    LOG.info("Bins: %d", len(bins))

    # Load contact samples from input simulation data and merge into a single
    # cooler dataset.

    def scan_pixels():
        for input_sim in input_sims:
            LOG.info("Ingesting from %s", input_sim)
            try:
                with h5py.File(input_sim, "r") as store:
                    snapshots = store["stages"]["interphase"]
                    steps = [step.decode() for step in snapshots[".steps"]]

                    steps_to_use = steps
                    if frames is not None:
                        steps_to_use = np.concatenate(
                            [steps[frame_slice] for frame_slice in frames]
                        )

                    for i in tqdm.trange(len(steps_to_use), leave=False):
                        step = steps_to_use[i]

                        # Contact maps are thinned out in the simulation
                        # output file because of size.
                        sample = snapshots[step]
                        if "contacts" not in sample:
                            continue
                        contacts = sample["contacts"][:]

                        yield {
                            "bin1_id": contacts[:, 0],
                            "bin2_id": contacts[:, 1],
                            "count": contacts[:, 2],
                        }
            except Exception as ex:
                LOG.warning(">> Skipping: %s", ex)

    cooler.create_cooler(
        output, bins, scan_pixels(), assembly=assembly,
    )

    LOG.info("Balancing contact matrix")
    clr = cooler.Cooler(output)
    cooler.balance_cooler(clr, store=True)


def parse_args() -> dict:
    parser = argparse.ArgumentParser(
        description="Collect contact samples from simulation trajectories",
    )
    arg = parser.add_argument

    arg(
        "--output",
        metavar="sim.cool",
        type=str,
        help="Cool dataset to write the contact matrix to",
    )
    arg(
        "--frames",
        metavar="300-500,600",
        type=str,
        help="Frame range to collect contact samples from",
    )
    arg(
        "input_sims",
        metavar="sim.h5",
        type=str,
        nargs="+",
        help="Simulation trajectories",
    )

    args = vars(parser.parse_args())
    reparse_ranges(args, "frames")
    return remove_none(args)


def reparse_ranges(args: dict, key: str):
    # "100,200-300,400-" -> [slice(100, 101), slice(200, 301), slice(400, None)]
    spec = args.get(key)
    if spec:
        seq: list[slice] = []
        for sub_spec in spec.split(","):
            match sub_spec.split("-"):
                case [point]:
                    i = int(point)
                    seq.append(slice(i, i + 1))
                case [start, ""]:
                    i = int(start)
                    seq.append(slice(i, None))
                case [start, end]:
                    i = int(start)
                    j = int(end)
                    seq.append(slice(i, j + 1))
        args[key] = seq


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
