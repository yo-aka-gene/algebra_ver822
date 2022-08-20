from typing import Any, Dict, List, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.core.series import Series
from pandas.core.frame import DataFrame
from scipy.interpolate import PchipInterpolator

from ._utils import kwarg_mgr


def _fmt_for_pcp(
    data: DataFrame,
    twin_label: bool = True,
    index: list = None,
    label_width: int = 1
):
    
    col_0 = data.columns[0]
    col_last = data.columns[-1]
    
    scale_b = data.min().min()
    scale_a = 1 * (data.max().max() - scale_b) / (len(data) - 1)
    order = (scale_a * np.arange(len(data)) + scale_b).reshape(-1, 1)
    
    if label_width > 1:
        order = pd.concat(
            [pd.DataFrame(order)] * label_width,
            axis=1
        )
        order.columns = ["."] * (label_width - 1) + ["last"]
    
    if twin_label:
        if index is None:
            mat = data.sort_values(col_last, ascending=True)
            
        else:
            mat = data.loc[index, :]

        if label_width > 1:
            order.index = mat.index
            mat = pd.concat([mat, order], axis=1)

        else:
            mat = mat.assign(last=order)
    
    else:
        mat = data
    
    if index is None:
        mat = mat.sort_values(col_0, ascending=True)
    
    else:
        mat = mat.loc[index, :]
    
    if label_width > 1:
        order.index = mat.index
        order.columns = [""] + ["."] * (label_width - 1)
        mat = pd.concat([order, mat], axis=1)
    
    else:
        mat.insert(0, "", order)
    
    return mat


def ticks(
    data: Union[DataFrame, Series]
) -> Dict[str, np.ndarray]:
    
    assert isinstance(data, (DataFrame, Series)), \
        f"data expected DataFrame or Series, got {data}"
    
    rule = {
        0.005: 0.01, 0.01: 0.02,
        0.02: 0.025, 0.025: 0.05,
        0.05: 0.1, 0.1: 0.2, 0.2: 0.25,
        0.25: 0.4, 0.5: 10
    }
    
    if isinstance(data, DataFrame):
        vmax = data.max().max()
        vmin = data.min().min()
        
    else:
        vmax = data.max()
        vmin = data.min()
    
    dw = vmax - vmin
    check_val = (
        np.array(
            list(rule.keys())
        )[np.array(list(rule.keys())) <= 0.2 * dw]
    ).max()
    
    if rule[check_val] >= 0.2 * dw:
        interval = rule[check_val]
    
    elif 0.2 * dw > 0.5:
        interval = 1
    
    else:
        check_val = (
            np.array(
                list(rule.keys())
            )[np.array(list(rule.keys())) > 0.2 * dw]
        ).min()
        interval = rule[check_val]
        
    org_ticks = interval * np.arange(
        (vmin // interval),
        int(vmax // interval) + 2
    )
    
    scaled_ticks = (org_ticks - vmin) / dw
    
    return {
        "yticks": scaled_ticks,
        "yticklabels": org_ticks
    }


def pcp(
    data,
    ax: Any = None,
    cmap: Union[str, List[Tuple[float]], DataFrame, Series, tuple] = "hsv",
    adjust_scale: bool =True,
    curvate: bool = False,
    twin_label: bool =True,
    n_dot: int = 1000,
    remove_ticklabels: bool = False,
    index: list = None,
    label_width: int = 1,
    **kwargs
):
    if ax is None:
        fig, ax = plt.subplots()
        
    assert issubclass(type(ax), mpl.axes.SubplotBase), \
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    assert isinstance(cmap, (str, list, DataFrame, Series, tuple)), \
        f"cmap expected str or List[Tuple[float]] or DataFrame or Series, got {cmap}"
    
    if adjust_scale:
        scalemat_a = np.tile((data.max() - data.min()).values, (len(data), 1))
        scalemat_b = np.tile(data.min().values, (len(data), 1))
        
        mat = (data - scalemat_b) / scalemat_a
        
        tick_info = list(map(
            (lambda x: ticks(x)),
            [data.loc[:, v] for v in data]
        ))
    
    else:
        mat = data
    
    mat = _fmt_for_pcp(mat, twin_label, index, label_width)
    
    if isinstance(cmap, DataFrame):
        cmap = cmap.loc[mat.index, :]
    
    elif isinstance(cmap, Series):
        cmap = cmap.loc[mat.index]
    
    cm = lambda x: eval("plt.cm." + x)
    colors = [cm(cmap)(i/len(data)) for i in range(len(data))] if isinstance(cmap, str) else cmap
    
    alpha = kwarg_mgr(kwargs, "alpha", 1)
    lw = kwarg_mgr(kwargs, "linewidth", 1)
    lw = kwarg_mgr(kwargs, "lw", lw)
    
    if curvate:
        curve = lambda xrange: (lambda vec: PchipInterpolator(np.arange(len(vec)), vec)(xrange))
        linsp = lambda n_dot, matrix: np.linspace(0, matrix.shape[1] - 1, n_dot)

        fc = lambda n_dot, matrix: np.array(list(map(curve(linsp(n_dot, matrix)), matrix)))
        curve_mat = fc(n_dot, mat.values)
        
        plotter = lambda i: ax.plot(
            linsp(n_dot, mat), curve_mat[i, :], color=colors[i],
            alpha=alpha, linewidth=lw
        ) 
        
        ax.set_xticks(np.arange(label_width, mat.shape[1] - label_width))
        ax.set_xticklabels(mat.columns[label_width:-label_width])
        
    else:
        plotter = lambda i: ax.plot(
            mat.columns, mat.iloc[i, :], color=colors[i],
            alpha=alpha, linewidth=lw
        ) 
        
    list(map(plotter, np.arange(len(mat))))
    ax.set_xlim(0, mat.shape[1] - 1)
    ax.set_yticks(mat[""])
    ax.set_yticklabels(mat.index)
    init_ylims = ax.get_ylim()
    
    if remove_ticklabels:
        ax.axes.yaxis.set_ticks([])
    
    if twin_label:
        ax2 = ax.twinx()
        ax2.set_ylim(*init_ylims)
        ax2.set_yticks(mat[""])
        ax2.set_yticklabels(mat.sort_values("last", ascending=True).index)
        
        if remove_ticklabels:
            ax2.axes.yaxis.set_ticks([])
    
    f = lambda x: (x - ax.get_xlim()[0]) / np.diff(ax.get_xlim())[0]

    for i in np.arange(label_width, data.shape[1] + label_width):
        ax_i = ax.twinx()
        ax_i.spines["right"].set_position(("axes", f(i)))
        ax_i.set_ylim(*init_ylims)
        
        if adjust_scale:
            if i <= len(data):
                within_area = (-0.05 < tick_info[i - label_width]['yticks']
                              ) * (
                   tick_info[i - label_width]['yticks'] < 1.05)
                ax_i.set_yticks(tick_info[i - label_width]["yticks"][within_area])
                ax_i.set_yticklabels((
                    tick_info[i - label_width]["yticklabels"][within_area]
                ).round(2))
            
    ax.set_ylim(*init_ylims)