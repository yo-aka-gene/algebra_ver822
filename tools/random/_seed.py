import numpy as np


def iterative_seed(
    n_iter: int,
    seed: int,
    min_val: int = 0,
    max_val: int = 10e6
) -> np.ndarray:
    np.random.seed(seed)
    return np.random.randint(min_val, max_val, n_iter)