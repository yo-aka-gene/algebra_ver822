import numpy as np
import pandas as pd

from ._eigenvals import eigvals
from ._null_models import null_model

def parallel_analysis(
    data: pd.core.frame.DataFrame,
    random_state: int,
    nullmodel: str = "norm",
    use_smc: bool = False
) -> int:
    assert nullmodel in ["norm", "perm", "min(norm,perm)"], \
        f"nullmodel expected 'norm', 'perm', or 'min(norm,perm)', got {nullmodel}"
    
    ev = eigvals(data, use_smc)
    
    if nullmodel == "min(norm,perm)":
        nm = null_model(data, random_state, "norm")
        pm = null_model(data, random_state, "perm")
        
        evn = eigvals(nm, use_smc)
        evp = eigvals(pm, use_smc)
        
        n_f = min(len(ev[ev > evn]), len(ev[ev > evp]))
    
    else:
        model = null_model(data, random_state, nullmodel)
        evn = eigvals(model, use_smc)
        
        n_f = len(ev[ev > evn])
    
    return n_f