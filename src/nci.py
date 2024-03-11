import argparse
import logging
import os
import signal

import cooler

from pkg.common.args import remove_none
from pkg.common.cli import invoke_main
from pkg.nci import make_nci_track


LOG = logging.getLogger()


def main(
    *,
    input_filename: str,
    output_filename: str,
    binsize: int | None = None,
    exclude: list[str] = [],
):
    input_dataset = input_filename
    if binsize:
        input_binsize = binsize // 2
        LOG.info("Binsize %d => contact map resolution is %d", binsize, input_binsize)
        input_dataset += f"::/resolutions/{input_binsize}"

    LOG.info("Opening %s", input_dataset)
    cool = cooler.Cooler(input_dataset)

    LOG.info("Writing to %s", output_filename)
    with open(output_filename, "w") as output:
        output_format = determine_format(output_filename)
        need_header = True
        for chrom in cool.chromnames:
            if chrom in exclude:
                continue

            LOG.info("Computing NCI on %s", chrom)
            track = make_nci_track(cool, chrom, halve=True)
            track.to_csv(output, index=False, header=need_header, **output_format)
            need_header = False

            output.flush()


def determine_format(filename: str) -> dict:
    kwargs = dict(sep="\t", float_format="%g", na_rep="nan")
    if filename.endswith(".csv"):
        kwargs["sep"] = ","
    return kwargs


def parse_args() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binsize", type=int)
    parser.add_argument("--exclude", type=str)
    parser.add_argument("--output", dest="output_filename", type=str, required=True)
    parser.add_argument("input_filename", type=str)

    args = vars(parser.parse_args())
    reparse_list(args, "exclude")
    return remove_none(args)


def reparse_list(args: dict, key: str):
    if arg := args.get(key):
        args[key] = arg.split(",")


if __name__ == "__main__":
    invoke_main(main, parse_args(), LOG)
