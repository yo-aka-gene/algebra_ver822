from typing import Any, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from tools.figure import kwarg_mgr


def silhouette_curve(
    l_silmean: list,
    l_silval: list,
    ax: Any,
    search_sp: list,
    **kwargs
):
    assert issubclass(type(ax), mpl.axes.SubplotBase),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    marker = kwarg_mgr(kwargs, "marker", "o")

    sns.lineplot(
        x=np.arange(search_sp[0], search_sp[1] + 1),
        y=l_silmean, marker="o", ax=ax,
        label="Mean"
    )

    sns.lineplot(
        x=np.arange(search_sp[0], search_sp[1] + 1),
        y=[v.std() for v in l_silval], marker="o", ax=ax,
        label="S.D."
    )
    
    ax.set(xlabel="Num. of $k$", ylabel="Silhouette coeff.")
    ax.set_xticks(np.arange(search_sp[0], search_sp[1] + 1, 2))
    ax.set_xticklabels(np.arange(search_sp[0], search_sp[1] + 1, 2));


def silhouette_plot(
    l_labels: list,
    ax: Any,
    cmap: list,
    silscore: float,
    margin: float = 0.1,
    **kwargs
) -> tuple:
    assert issubclass(type(ax), mpl.axes.SubplotBase),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    color = kwarg_mgr(kwargs, "color", ".2")
    linewidth = kwarg_mgr(kwargs, "linewidth", 1)
    linestyle = kwarg_mgr(kwargs, "linestyle", "--")
    
    l_y = np.array([
        len(v) for v in ([[]] + l_labels)
    ]).cumsum()
    
    margin = l_y[-1] * margin
    
    for i, v in enumerate(l_labels):
        ax.fill_betweenx(
            y=np.arange(l_y[i] + i * margin, l_y[i + 1] + i * margin)[:len(v)],
            x1=v.iloc[:, -1].values,
            x2=0,
            color=cmap[i]
        )

    ax.axes.yaxis.set_ticks([])
    ax.set_ylim([l_y[0] - margin, l_y[-1] + margin * (len(l_labels) + 1)])

    f = lambda x: (x - ax.get_xlim()[0]) / np.diff(ax.get_xlim())[0]

    ax2 = ax.twinx()
    ax2.spines["left"].set_position(("axes", f(0)))
    ax2.yaxis.set_ticks_position("left")
    ax2.set_ylim([*ax.get_ylim()])
    ax2.set_yticks(
        np.convolve(
            l_y, [0.5, 0.5]
        )[1:-1] + np.array([
            margin * i for i in range(len(l_labels))
        ])
    )
    ax2.set_yticklabels([f"cluster{i}" for i in range(len(l_labels))]);

    ax.vlines(
        silscore, *ax.get_ylim(),
        color=color, linewidth=linewidth, linestyle=linestyle, label="Mean",
    )

    ax.legend()
    
    return ax, ax2