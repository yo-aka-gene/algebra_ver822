from typing import Any, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from sklearn.cluster import KMeans

from tools.figure import kwarg_mgr


def elbowplot(
    l_km: List[KMeans],
    ax: Any,
    search_sp: list,
    show_diff: bool = False,
    **kwargs
):
    assert issubclass(type(ax), mpl.axes.SubplotBase) or isinstance(ax, np.ndarray),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    if isinstance(ax, np.ndarray):
        for v in ax.ravel():
            assert issubclass(type(v), mpl.axes.SubplotBase),\
                f"ax expected matplotlib.axes._subplots.AxesSubplot"

    elbow = np.array([model.inertia_ for model in l_km])
    
    marker = kwarg_mgr(kwargs, "marker", "o")

    ax = ax.ravel() if isinstance(ax, np.ndarray) else np.array([ax])

    sns.lineplot(
        x=np.arange(search_sp[0], search_sp[1] + 1),
        y=elbow,
        marker=marker, ax=ax[0],
        label="$Inertia$"
    )

    ax[0].set_xticks(np.arange(search_sp[0], search_sp[1], 2))
    ax[0].set_xticklabels(np.arange(search_sp[0], search_sp[1], 2))
    ax[0].set(ylabel="Inertia")

    if show_diff:
        sns.lineplot(
            x=np.arange(search_sp[0] + 1, search_sp[1] + 1),
            y=elbow[:-1] - elbow[1:],
            marker=marker, ax=ax[1], label="$Inertia_{k-1}-Inertia_{k}$"
        )

        ax[1].set_xticks(np.arange(search_sp[0], search_sp[1], 2))
        ax[1].set_xticklabels(np.arange(search_sp[0], search_sp[1], 2))
        ax[1].set(xlabel="Num. of $k$", ylabel="$\Delta$ Inertia");