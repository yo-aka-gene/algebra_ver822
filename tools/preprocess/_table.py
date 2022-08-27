import os
from typing import List
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tqdm.notebook import tqdm

from  tools.r import list2tsv

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


def fmt_mtx(
    l_path: List[str],
    save_dir: str,
    fmt: str = "%.18e",
    axis: int = 0,
    mode: str = "py2r"
):
    assert mode in ["py2r", "py2py", "r2py", "r2r"], \
        f"invalid input: expected either of 'py2r', 'py2py', 'r2py', or 'r2r', got {mode}"
    n_r, n_c, n_nz = 0, 0, 0

    for i, v in tqdm(
        enumerate(l_path), desc="Concatenation", total=len(l_path)
        ):
        txt = open(f"{v}/matrix.mtx").readlines()
        info = txt[2] if " " not in txt[1] else txt[1]
        n_row, n_col, n_nonzero = map(lambda x: np.int32(x), info.split("\n")[0].split(" "))
        
        data = np.loadtxt(
            f"{v}/matrix.mtx", delimiter=" ", skiprows=3, dtype=np.int32
        )
        n_r += n_row if axis == 0 else 0
        n_c += n_col if axis == 1 else 0
        n_nz += n_nonzero
        archive = data if i == 0 else np.vstack([
            archive, data + np.tile([n_r, 0, 0] if axis == 0 else [0, n_c, 0], (n_nonzero, 1))
        ])
    
    if mode in ["py2py", "r2r"]:
        header = f"%%MatrixMarket matrix coodinate integer general\n{n_r} {n_c} {n_nz}"
    
    else:
        for i in tqdm([0], desc=f"Adjusting format", total=1):
            archive = archive[:, [1, 0, 2]]
            header = f"%%MatrixMarket matrix coordinate integer general\n{n_c} {n_r} {n_nz}"

    for i in tqdm([0], desc=f"Exporting log", total=1):
        np.savetxt(f"{save_dir}/matrix.mtx", archive, delimiter=" ", fmt=fmt, header=header, comments="")


def fmt_tsv(
    l_path: List[str],
    filenames: str,
    save_dir: str,
    unique: bool = True,
    concat: bool = True
):
    if concat:
        ret = []
        for v in tqdm(
            l_path, desc="Concatenation", total=len(l_path)
            ):
            ret += [
                np.loadtxt(
                    f"{v}/{filenames}.tsv", delimiter="\t", dtype=str
                ).ravel()
            ]
        ret = np.unique(np.concatenate(ret)) if unique else np.concatenate(ret)

    else:
        ret = np.loadtxt(
            f"{l_path[0]}/{filenames}.tsv", delimiter="\t", dtype=str
        ).ravel()
    for i in tqdm([0], desc=f"Exporting {filenames}.tsv", total=1):
        list2tsv(ret.tolist(), f"{save_dir}/{filenames}.tsv")
