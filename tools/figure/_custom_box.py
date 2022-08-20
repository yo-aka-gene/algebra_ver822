from typing import Any, Union, List, Tuple
import matplotlib as mpl
import pandas as pd
import seaborn as sns

from tools.figure import kwarg_mgr


def boxplot(
    x: str = None,
    y: str = None,
    hue: str = None,
    data: pd.core.frame.DataFrame = None,
    order: List[str] = None,
    hue_order: List[str] = None,
    orient: str = None,
    color: Any = None,
    palette: Union[str, List[Tuple[float]]] = None,
    saturation: float = 0.75,
    width: float = 0.9,
    dodge: bool = True,
    fliersize: float = 1,
    linewidth: float = 1,
    whis: Union[float, Tuple[float]] = (0, 100),
    ax: Any = None,
    **kwargs
):
    alpha = kwarg_mgr(kwargs, "alpha", 1)
    edgecolor = kwarg_mgr(kwargs, "edgecolor", "k")
    color = edgecolor if color is None else color
    
    line = dict(color=edgecolor, linewidth=linewidth)
    flier = dict(color=edgecolor, markersize=fliersize)
    box = dict(edgecolor=edgecolor)
    
    args = dict(
        x=x, y=y, hue=hue, data=data, order=order,
        hue_order=hue_order, orient=orient, color=color,
        palette=palette, saturation=saturation, width=width,
        dodge=dodge, fliersize=fliersize, linewidth=linewidth,
        whis=whis, ax=ax, medianprops=line,
        whiskerprops=line, capprops=line,
        flierprops=flier,
    )
    
    sns.boxplot(**{**args, **dict(boxprops=dict(**box, facecolor="w"))})
    sns.boxplot(**{**args, **dict(boxprops=dict(**box, alpha=alpha))})