from itertools import combinations
from typing import Any, List

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pgmpy.estimators import PC


class PCGRN():
    def __init__(
        self,
        data: pd.core.frame.DataFrame
    ):
        self.data = data
    
    def estimate(
        self,
        variant: str = "stable",
        ci_test: str = "pearsonr",
        max_cond_vars: int = None,
        return_type: str = "dag",
        significance_level: float = 0.01,
        n_jobs: int = -1,
        show_progress: bool = False
    ):
        max_cond_vars=self.data.shape[1] if max_cond_vars is None else max_cond_vars
        
        data = self.data.loc[:, self.data.columns[self.data.max() > 0]]

        model = PC(data=data).estimate(
            variant=variant,
            ci_test=ci_test,
            max_cond_vars=max_cond_vars,
            return_type=return_type,
            significance_level=significance_level,
            n_jobs=n_jobs,
            show_progress=show_progress
        )
        
        self.edges = list(model.edges)
    
    def get_matrix(self) -> pd.core.frame.DataFrame:
        df = pd.DataFrame(
            np.eye(self.data.shape[1]),
            index=self.data.columns,
            columns=self.data.columns
        )
        
        for i, v in combinations(df.columns, 2):
            df.loc[i, v] = 1 if (i, v) in self.edges else 0
        
        return df
    
    def plot(
        self,
        ax=None,
        c="C0"
    ):
        if ax is None:
            fig, ax = plt.subplots()
        
        feat = self.data.columns
        feat_dict = {v: i for i, v in enumerate(feat)}
        n = len(feat)

        t = np.linspace(0, 2 * np.pi, num=n + 1)
        x = 10 * np.cos(t)
        y = 10 * np.sin(t)

        ax.scatter(x, y, color=c)


        for i, v in enumerate(feat):
            ax.annotate(
                v,
                (1.1 * x[i] - 1, 1.1 * y[i] - .4),
                fontsize="x-small"
            )

        for v in self.edges:
            ax.plot([x[feat_dict[i]] for i in v], [y[feat_dict[i]] for i in v], color=c)

        ax.set_xlim([-12, 12])
        ax.set_ylim([-12, 12])

        ax.axis("off")
