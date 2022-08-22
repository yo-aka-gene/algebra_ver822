import os
from typing import Callable, List, Union
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tqdm.notebook import tqdm

def submatrix(
    l_data: List[str],
    save_dir: str,
    filenames: str,
    rownames: List[str] = None,
    colnames: List[str] = None,
    index_col: int = 0,
    log2: bool = False,
) -> List[str]:
    
    assert isinstance(filenames, str), \
        f"filenames should be a string"
    
    header = open(l_data[0]).readline().split(",")[1:]

    if colnames is None:
        col = None
    else:
        col = np.array(
            [-1] + [int(np.where(np.array(header)==v)[0]) for v in colnames]
        ) + 1
    
    ret = []
    
    for i, v in tqdm(enumerate(l_data), desc='extraction', total=len(l_data)):
        temp = pd.read_csv(v, index_col=index_col, usecols=col)
        row = temp.index if rownames is None else np.intersect1d(rownames, temp.index)
        temp = temp.loc[row, :]

        if log2:
            temp = np.log2(temp + 1)
            temp.index.name = "$\log_2(" + f"{temp.index.name}" + "+1)$"
        
        temp.to_csv(f"{save_dir}/{filenames}.{v.split('.')[-1]}", index=True)
        ret += [f"{save_dir}/{filenames}.{v.split('.')[-1]}"]
        
    return ret


def clear_artifacts(l_data: List[str]) -> None:
     for i, v in tqdm(enumerate(l_data), desc='removing artifacts', total=len(l_data)):
            os.remove(v)


def fmt_table(
    l_data: List[str],
    save_dir: str,
    filenames: str,
    rownames: List[str] = None,
    colnames: List[str] = None,
    index_col: int = 0,
    log2: bool = False,
    axis: int = 0,
    export: bool = False,
    export_fmt: str = "pkl",
    save_as_csv: bool = False,
) -> pd.core.frame.DataFrame:
    l_path = submatrix(
        l_data,
        save_dir,
        filenames,
        rownames,
        colnames,
        index_col,
        log2
    )
    
    df = pd.concat(
        [pd.read_csv(v, index_col=index_col) for v in l_path],
        axis = axis
    )
    
    clear_artifacts(l_path)
    
    if save_as_csv:
        df.to_csv(f"{save_dir}/{filenames}.csv", index=True)
    
    elif export and export_fmt == "pkl":
        df.to_pickle(f"{save_dir}/{filenames}.pkl", index=True)
    
    elif export and export_fmt == "npz":
        sp.save_npz(f"{save_dir}/{filenames}.npz", sp.csc_matrix(df.values))

    return df