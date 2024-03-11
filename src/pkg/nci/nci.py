import cooler
import numpy as np
import pandas as pd


DEFAULT_CHUNK_SIZE = 512


def compute_nci(
    matrix,
    start: int,
    end: int,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> np.ndarray:
    chunk_ncis = [np.empty(0)]

    for offset in range(start, end, chunk_size):
        # Must extend the end of the chunk (if it is not the last one) so that
        # the bin pair lying between the next chunk is handled.
        chunk_slice = slice(offset, min(offset + chunk_size + 1, end))
        chunk = matrix[chunk_slice, chunk_slice]
        diag = np.diag(chunk)
        sub = np.diag(chunk, 1)

        # Allow NaNs at zero-read sites.
        with np.errstate(divide="ignore", invalid="ignore"):
            nci = sub / np.sqrt(diag[1:] * diag[:-1])

        chunk_ncis.append(nci)

    return np.concatenate(chunk_ncis)


def make_nci_track(
    cool: cooler.Cooler,
    chrom: str,
    *,
    halve: bool = False,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> pd.DataFrame:
    def make_track(start, end, score):
        return pd.DataFrame(
            {"chrom": chrom, "start": start, "end": end, "score": score}
        )

    # No balancing is needed as NCI is invariant under multiplicative bias.
    cool_matrix = cool.matrix()
    cool_bins = cool.bins()
    chrom_start, chrom_end = cool.extent(chrom)

    if chrom_start == chrom_end:
        return make_track([], [], [])

    bins = cool_bins[chrom_start:chrom_end]
    nci = compute_nci(cool_matrix, chrom_start, chrom_end, chunk_size=chunk_size)
    assert len(nci) == len(bins) - 1

    # NCI at the i-th position characterizes the coalesced region spanning
    # across the i-th and the (i+1)-th input bins. Set the start-end pair
    # accordingly.
    #
    #     | * . . . . .
    #     |   * . . . .
    #  i  |     x s . . ] nci(i) = s / sqrt(x*y)
    # i+1 |       y . . ]
    #     |         * .
    #     V           *
    #
    if len(nci) > 0:
        track = make_track(
            start=bins["start"].values[:-1],
            end=bins["end"].values[1:],
            score=nci,
        )
    else:
        track = make_track(
            start=bins["start"],
            end=bins["end"],
            score=np.nan,
        )

    # Drop even rows so that the output has no overlapping bins.
    if halve:
        halved_track = track[::2]

        if len(track) % 2 == 0:
            remains = track[-1:]
            prev_end = track.at[track.index[-2], "end"]
            remains.at[remains.index[0], "start"] = prev_end
            halved_track = pd.concat([halved_track, remains], ignore_index=True)

        track = halved_track

    return track
