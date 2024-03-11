import typing

import numpy as np


class SVDStep(typing.NamedTuple):
    step: int
    vector: np.ndarray
    delta: float


def power_svd(
    data: np.ndarray,
    init: np.ndarray | None = None,
) -> typing.Iterator[SVDStep]:
    """
    Iteratively computes the first right singular vector of the given matrix.
    """
    step = 0
    prev_vec = init or _svd_init(data)
    while True:
        vec = _svd_iter(data, prev_vec)
        step += 1
        yield SVDStep(step=step, vector=vec, delta=np.abs(vec - prev_vec).max())
        prev_vec = vec


def _svd_init(data: np.ndarray) -> np.ndarray:
    data_dim = data.shape[1]
    return np.ones(data_dim, dtype=data.dtype) / np.sqrt(data_dim)


def _svd_iter(data: np.ndarray, vec: np.ndarray) -> np.ndarray:
    weights = data @ vec
    new_vec = np.nansum(weights[:, None] * data, axis=0)
    new_vec /= np.linalg.norm(new_vec)
    return new_vec
