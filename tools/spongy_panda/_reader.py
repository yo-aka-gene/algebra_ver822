import json
import pickle

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from ._sparse import SparseDF
from tools.r import read_json


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

def load_mtx(
    filename: str,
    colidx_json: str = None,
    from_r: bool = False,
    transpose: bool = None,
) -> SparseDF:
    transpose = from_r if transpose is None else transpose
    ret = SparseDF(
        mmread(filename).tocsc(),
        **read_json(colidx_json, from_r)
        )
    return ret.t() if transpose else ret


def _read_tsv(filename: str) -> list:
    data = pd.read_csv(filename, sep="\t")
    return data.columns.tolist() + data.values.ravel().tolist()


def load_10xdir(
    path: str,
    from_r: bool = False,
    transpose: bool = None
    ) -> SparseDF:
    transpose = from_r if transpose is None else transpose

    mtx = mmread(f"{path}/matrix.mtx").tocsc()

    return SparseDF(
        mtx.transpose() if transpose else mtx,
        index=_read_tsv(f"{path}/barcodes.tsv")
        columns=_read_tsv(f"{path}/features.tsv")
    )
