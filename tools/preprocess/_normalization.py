import numpy as np
import pandas as pd

from tqdm.notebook import tqdm
from typing import List

def count2rpm(data: pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    return 10e5 * (data.T / data.sum(axis=1)).T

def fmt_rpm(
    l_data: List[str],
    save_dir: str,
    index_col: int = 0,
    log2: bool = False,
) -> None:
    head = "$\log_2(RPM+1)$" if log2 else "RPM"
    names = [head] + open(l_data[0]).readline().split(",")[1:]
    
    filename = "log2rpm" if log2 else "rpm"
    
    for i, v in tqdm(enumerate(l_data), desc='RPM normalization', total=len(l_data)):
        if i == 0:
            temp = count2rpm(pd.read_csv(v, index_col=index_col))
            temp.index.name = head
        
        else:
            temp = count2rpm(pd.read_csv(v, index_col=index_col, names=names))
        
        temp = np.log2(temp + 1) if log2 else temp
        temp.to_csv(f"{save_dir}/{filename}.{v.split('.')[-1]}", index=True)