from typing import Any, List, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from scipy.stats import entropy
from sklearn.cluster import KMeans

from tools.figure import kwarg_mgr


def entropy_plot(
    l_km: List[KMeans],
    ax: Any,
    search_sp: list,
    return_array: bool = False,
    **kwargs
) -> Union[np.ndarray, None]:
    assert issubclass(type(ax), mpl.axes.SubplotBase),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    marker = kwarg_mgr(kwargs, "marker", "o")
    color = kwarg_mgr(kwargs, "color", "C2")
    
    l_entropy = [
        entropy(
            np.bincount(v.labels_) / np.bincount(v.labels_).sum(),
        ) for v in l_km
    ]

    sns.lineplot(
        x=np.arange(search_sp[0], search_sp[1] + 1),
        y=l_entropy, 
        marker=marker, ax=ax, color=color, label="entropy"
    )

    ax.set(xlabel="Num. of $k$", ylabel="entropy");
    ax.set_xticks(np.arange(search_sp[0], search_sp[1] + 1, 2))
    ax.set_xticklabels(np.arange(search_sp[0], search_sp[1] + 1, 2))
    
    return l_entropy if return_array else None