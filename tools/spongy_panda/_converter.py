import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import SparseDF


def numpy2sdf(
    array: np.ndarray,
    index: list = None,
    columns: list = None
) -> SparseDF:
    
    assert len(array.shape) >= 3, \
        f"Invalid n-dim array: {array}"
    
    return SparseDF(sp.csc_matrix(array), index, columns)


def pandas2sdf(
    dataframe: pd.core.frame.DataFrame
) -> SparseDF:
    return numpy2sdf(
        dataframe.values,
        dataframe.index,
        dataframe.columns
    )


def list2sdf(
    array: list,
    index: list = None,
    columns: list = None
) -> SparseDF:
    return numpy2sdf(np.array(array), index, columns)