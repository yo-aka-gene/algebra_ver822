import json
import pickle

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._sparse import SparseDF


def read_csv(
    filename: str,
    index_col: int = None,
    names: np.ndarray = None,
    usecols: list = None
) -> SparseDF:
    if names is None:
        df = pd.read_csv(
            filename,
            index_col=index_col,
            usecols=usecols
        )
    else:
        df = pd.read_csv(
            filename,
            index_col=index_col,
            names=names,
            usecols=usecols
        )
    
    return SparseDF(
        sp.csc_matrix(df.values),
        df.index, df.columns
    )


def read_pickle(
    filename: str,
    colidx_json: str = None
) -> SparseDF:
    with open(filename, "rb") as f:
        df = pickle.load(f)
    
    if isinstance(df, SparseDF):
        return df
    
    elif colidx_json is not None:
        with open(colidx_json, "r") as g:
            d = json.load(g)
    
    else:
        d = dict(index=None, columns=None)
    
    assert isinstance(
        df, (
            sp._csc.csc_matrix,
            np.ndarray,
            pd.core.frame.DataFrame,
            pd.core.series.Series
        )
    ), f"invalid dtype: {type(df)}"
    
    if isinstance(df, sp._csc.csc_matrix):
        return SparseDF(df, **d)
    
    elif isinstance(df, np.ndarray):
        return SparseDF(sp.csc_matrix(df), **d)
    
    elif isinstance(df, pd.core.series.Series):
        return SparseDF(df, df.index, df.name)
    
    else:
        return SparseDF(df, df.index, df.columns)


def load_npz(
    filename: str,
    colidx_json: str = None
) -> SparseDF:
    df = sp.load_npz(filename)
    
    if colidx_json is not None:
        with open(colidx_json, "r") as g:
            d = json.load(g)
    
    else:
        d = dict(index=None, columns=None)
        
    return SparseDF(df, **d)
    