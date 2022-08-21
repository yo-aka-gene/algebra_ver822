import numpy as np
import pandas as pd

import scipy.sparse as sp
from tqdm.notebook import tqdm
from typing import List

import tools.spongy_panda as spd

def count2rpm(data: pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    return 10e5 * (data.T / data.sum(axis=1)).T

def count2rpm_sdf(data: spd.SparseDF) -> spd.SparseDF:
    return 10e5 * (data.t() / data.sum(axis=1)).t()
    

def fmt_rpm(
    l_data: List[str],
    save_dir: str,
    index_col: int = 0,
    log2: bool = False,
    save_one_by_one: bool = False,
) -> None:
    head = "$\log_2(RPM+1)$" if log2 else "RPM"
    names = [head] + open(l_data[0]).readline().split(",")[1:]
    
    filename = "log2rpm" if log2 else "rpm"
    ret = []
    
    for i, v in tqdm(enumerate(l_data), desc='RPM normalization', total=len(l_data)):
        if i == 0:
            temp = pd.read_csv(v, index_col=index_col)
            temp.index.name = head
        
        else:
            temp = pd.read_csv(v, index_col=index_col, names=names)
        
        temp = np.log2(count2rpm(temp) + 1) if log2 else count2rpm(temp)
        
        if save_one_by_one:
            temp.to_pickle(f"{save_dir}/{filename}_{i}.pkl")
        
        else:
            ret += [temp]
            
    if not save_one_by_one:
        pd.concat(ret).to_pickle(f"{save_dir}/{filename}.pkl")

            
def fmt_rpm_sdf(
    l_data: List[str],
    save_dir: str,
    index_col: int = 0,
    log2: bool = False,
    save_one_by_one: bool = False,
) -> None:
    head = "$\log_2(RPM+1)$" if log2 else "RPM"
    names = [head] + open(l_data[0]).readline().split(",")[1:]
    
    filename = "log2rpm" if log2 else "rpm"
    ret = []
    
    for i, v in tqdm(enumerate(l_data), desc='RPM normalization', total=len(l_data)):
        if i == 0:
            temp = pd.read_csv(v, index_col=index_col)
            temp.index.name = head
        
        else:
            temp = pd.read_csv(v, index_col=index_col, names=names)
        
        temp = count2rpm_sdf(spd.pandas2sdf(temp))
        temp = spd.SparseDF(
            sp.csc_matrix(np.log2(temp + 1)), 
            index=temp.index, 
            columns=temp.columns
        ) if log2 else temp
        
        
        if save_one_by_one:
            temp.to_pickle(f"{save_dir}/{filename}_{i}.pkl")
        
        else:
            ret += [temp()]
        
    if not save_one_by_one:
        spd.concat(ret).to_pickle(f"{save_dir}/{filename}.pkl")