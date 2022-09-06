from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tools.figure import kwarg_mgr
from ._eigenvals import eigvals
from ._null_models import null_model
from ._pa import parallel_analysis

def plot_parallel_analysis(
    data: pd.core.frame.DataFrame,
    random_state: int,
    nullmodel: str = "norm",
    use_smc: bool = False,
    ax: Any = None,
    **kwargs
) -> None:
    assert nullmodel in ["norm", "perm", "min(norm,perm)"], \
        f"nullmodel expected 'norm', 'perm', or 'min(norm,perm)', got {nullmodel}"

    if ax is None:
        fig, ax = plt.subplots()
        
    assert issubclass(type(ax), mpl.axes.SubplotBase), \
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    ev = eigvals(data, use_smc)
    n_f = parallel_analysis(data, random_state, nullmodel, use_smc)
    
    c_factor = kwarg_mgr(kwargs, "c_factor", "r")
    c_resid = kwarg_mgr(kwargs, "c_resid", ".8")
    marker = kwarg_mgr(kwargs, "marker", "o")
    
    ax.plot(
        1 + np.arange(n_f), ev[:n_f],
        marker=marker, c=c_factor
    )
    ax.plot(
        1 + np.arange(n_f - 1, len(ev)), ev[n_f - 1:],
        marker=marker, c=c_resid, zorder=0
    )
    
    labels = {
        "norm": "null model$\sim \mathcal{N}(\mathbf{0}, I_n)$",
        "perm": "random permutation"
    }
    
    if nullmodel == "min(norm,perm)":
        evn = eigvals(
            null_model(data, random_state, "norm"),
            use_smc
        )
        evp = eigvals(
            null_model(data, random_state, "perm"),
            use_smc
        )
        
        c_norm = kwarg_mgr(kwargs, "c_norm", "k")
        c_perm = kwarg_mgr(kwargs, "c_perm", ".3")
        c_norm = kwarg_mgr(kwargs, "color", c_norm)
        c_perm = kwarg_mgr(kwargs, "color", c_perm)
        
        lw_norm = kwarg_mgr(kwargs, "lw_norm", 1)
        lw_perm = kwarg_mgr(kwargs, "lw_perm", 1)
        lw_norm = kwarg_mgr(kwargs, "linewidth", lw_norm)
        lw_perm = kwarg_mgr(kwargs, "linewidth", lw_perm)
        
        ax.plot(
            1 + np.arange(len(evn)), evn,
            color=c_norm, linewidth=lw_norm,
            label=labels["norm"]
        )
        ax.plot(
            1 + np.arange(len(evp)), evp,
            color=c_perm, linewidth=lw_perm,
            label=labels["perm"]
        )
    
    else:
        evn = eigvals(
            null_model(data, random_state, nullmodel),
            use_smc
        )

        c = kwarg_mgr(kwargs, "c_" + nullmodel, "k")
        c = kwarg_mgr(kwargs, "color", c)
        lw = kwarg_mgr(kwargs, "lw_" + nullmodel, 1)
        lw = kwarg_mgr(kwargs, "linewidth", lw)
        
        ax.plot(
            1 + np.arange(len(evn)), evn,
            color=c, linewidth=lw,
            label=labels[nullmodel]
        )
        
    ax.legend()
    
    xlabel = kwarg_mgr(kwargs, "xlabel", "Num. of factors")
    ylabel = kwarg_mgr(kwargs, "ylabel", "eigenvalues")
    title = kwarg_mgr(kwargs, "title", "Parallel Analysis ")

    ax.set(
        xlabel=xlabel, ylabel=ylabel,
        title=f"{title}({n_f} factors)"
    );
