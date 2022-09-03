from typing import Any

import factor_analyzer as fa
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from tools.figure import custom_bwr, kwarg_mgr, sns_color_mgr

def kmo_viz(
    data: pd.core.frame.DataFrame,
    ax: Any = None,
    criteria: float = 0.6,
    line: bool = False,
    landscape: bool = False,
    **kwargs
) -> None:
    if ax is None:
        fig, ax = plt.subplots()
        
    assert issubclass(type(ax), mpl.axes.SubplotBase), \
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
        
    vsep, vall = fa.calculate_kmo(data)
    
    rev = kwarg_mgr(kwargs, "reverse", False)
    sat = kwarg_mgr(kwargs, "saturation", 0.4)
    edgecolor = kwarg_mgr(kwargs, "edgecolor", ".5")
    digit = kwarg_mgr(kwargs, "round", 3)
    label = kwarg_mgr(kwargs, "xlabel", "MSA")
    title = kwarg_mgr(
        kwargs,
        "title", 
        "overall KMO index: "
    )
    
    color_kwargs = sns_color_mgr(
        kwargs,
        {"palette": custom_bwr(vsep, criteria, sat, rev)}
    )
    
    color, palette = color_kwargs["color"], color_kwargs["palette"]
    
    sns.barplot(
    x=data.columns if landscape else vsep,
    y=vsep if landscape else data.columns,
    orient='v' if landscape else 'h', ax=ax,
    palette=palette,
    color=color,
    edgecolor=edgecolor
    )

    vmin, vmax = ax.get_xlim() if landscape else ax.get_ylim()
    
    if not landscape:
        if line:
            ls = kwarg_mgr(kwargs, "linestyles", "--")
            lw = kwarg_mgr(kwargs, "linewidth", 1)
            lc = kwarg_mgr(kwargs, "linecolor", "k")
            ax.vlines(criteria, vmin, vmax, linestyles=ls, color=lc, linewidth=lw)
            
        ax.set(xlabel=label, ylim=[vmin, vmax], title=f"{title}{vall.round(digit)}")
    
    else:
        if line:
            ls = kwarg_mgr(kwargs, "linestyles", "--")
            lw = kwarg_mgr(kwargs, "linewidth", 1)
            lc = kwarg_mgr(kwargs, "linecolor", "k")
            ax.hlines(criteria, vmin, vmax, linestyles=ls, color=lc, linewidth=lw)
            
        ax.set(ylabel=label, xlim=[vmin, vmax], title=f"{title}{vall.round(digit)}")
