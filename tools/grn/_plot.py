from typing import Any, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tools.figure import kwarg_mgr
from ._eval import d_asterisk

def annot_plot(
    key_sbj: Union[str, int],
    dict_sbj: Dict[Union[str, int], pd.core.frame.DataFrame],
    dict_obj: Dict[Union[str, int], pd.core.frame.DataFrame],
    dict_c_sbj: dict,
    dict_c_obj :dict,
    ax: Any = None,
    **kwargs
):
    
    if ax is None:
        fig, ax = plt.subplots()
    
    ax.scatter(0, 0, color=dict_c_sbj[key_sbj])
    ax.annotate(key_sbj, (0, 0))
    
    dist = np.array([d_asterisk(dict_sbj[key_sbj], dict_obj[k]) for k in dict_obj])
    
    linecolor = kwarg_mgr(kwargs, "linecolor", ".8")
    
    t = np.linspace(0, 2 * np.pi, num=len(dict_obj) + 2)
    x = np.cos(t)
    y = np.sin(t)
    u = np.linspace(0, 2 * np.pi, num=1000)
    
    for i, v in enumerate(dist):
        ax.scatter(v * x[i + 1], v * y[i + 1], color=dict_c_obj[list(dict_obj.keys())[i]])
        ax.annotate(list(dict_obj.keys())[i], (v * x[i  +1], v * y[i + 1]), ((v + .01) * x[i + 1], (v + .01) * y[i + 1]))
        ax.plot(v * np.cos(u), v * np.sin(u), color=linecolor, zorder=0)
    
    lim = np.abs([*ax.get_xlim()] + [*ax.get_ylim()]).max()
    
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
