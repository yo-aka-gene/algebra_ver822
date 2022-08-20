from typing import Any, List, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from tools.figure import kwarg_mgr


def countplot(
    data: pd.core.frame.DataFrame,
    label: pd.core.series.Series,
    cmap: Union[str, List[Tuple[float]]],
    ax: Any = None,
    sort: bool = False,
    **kwargs
):
    if ax is None:
        fig, ax = plt.subplots()
        
    assert issubclass(type(ax), mpl.axes.SubplotBase), \
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
    label = label.sort_values() if sort else label
    data = data.loc[label.index, :].assign(group=label).groupby("group").count()
    data = data if sort else data.loc[label.unique(), :]
    
    if isinstance(cmap, str):
        cmap = [
            eval("plt.cm" + cmap)(
                i/len(label.unique())
            ) for i in range(len(label.unique()))
        ]
    
    alpha = kwarg_mgr(kwargs, "alpha", 1)
    edgecolor = kwarg_mgr(kwargs, "edgecolor", ".2")
    
    sns.barplot(
        data=data, x=data.index, y=data.iloc[:, 0],
        palette=cmap, alpha=alpha, edgecolor=edgecolor,
        ax=ax
    )
    
    ax.set(ylabel="Counts", xlabel="")