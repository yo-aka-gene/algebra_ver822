import numpy as np
import pandas as pd

def random_permutation(
    data: pd.core.frame.DataFrame,
    random_state: int
) -> np.ndarray:
    np.random.seed(random_state)
    n_row, n_col = data.shape
    ret = data.values.ravel()
    np.random.shuffle(ret)
    return ret.reshape(n_row, n_col)

def null_model(
    data: pd.core.frame.DataFrame,
    random_state: int,
    mode: str = "norm"
) -> np.ndarray:
    assert mode in ["norm", "perm"], \
        f"mode expected 'norm' or 'perm', got {mode}"
    
    if mode == "null":
        np.random.seed(random_state)
        model = np.random.multivariate_normal(
            np.zeros(data.shape[1]),
            np.eye(data.shape[1]),
            data.shape[0]
        )
    
    else:
        model = random_permutation(data, random_state)
    
    return model