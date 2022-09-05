import factor_analyzer as fa
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from tools.figure import kwarg_mgr
from ._pa import parallel_analysis

class FactorAnalyzer():
    def __init__(
        self,
        n_factor: int,
        data: pd.core.frame.DataFrame
    ):
        self.models = {}
        self.results = {}
        self.loadings = {}
        self.communalities = {}
        self.uniquenesses = {}
        self.figs = {}
        self.data = data
        self.n_f = n_factor
    
    def rotate(
        self,
        rotation: str,
        method: str = "minres",
        n_factor: int = None
    ):
        if n_factor is None:
            n_factor = self.n_f

        model = fa.FactorAnalyzer(
            n_factors=n_factor,
            rotation=rotation,
            method=method
        ).fit(self.data)
        
        result = pd.DataFrame(
            model.transform(self.data),
            index=self.data.index,
            columns=[f"factor{i+1}" for i in range(n_factor)]
        )
        
        loading = pd.DataFrame(
            model.loadings_,
            index=self.data.columns,
            columns=[f"factor{i+1}" for i in range(n_factor)],
        )
        
        com = model.get_communalities()
        uni = model.get_uniquenesses()
        
        condition = (n_factor, rotation, method)

        self.models = {
            **self.models,
            condition: model
        }
        self.model = model
        self.results = {
            **self.results,
            condition: result
        }
        self.result = result
        self.loadings = {
            **self.loadings,
            condition: loading
        }
        self.loading = loading
        self.communalities = {
            **self.communalities,
            condition: com
        }
        self.communality = com
        self.uniquenesses = {
            **self.uniquenesses,
            condition: uni
        }
        self.uniqueness = uni


    def heatmap(self, ax: None, condition: tuple = None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
            
            fig_id = (
                "heatmap",
                list(self.loadings.keys())[-1]
            ) if condition is None else ("heatmap", condition)
            
            self.figs = {**self.figs, fig_id: (fig, ax)}
        
        if (condition is not None) and (condition not in self.loadings):
            self.rotate(*list(condition))
            
        assert issubclass(type(ax), mpl.axes.SubplotBase),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
        cmap = kwarg_mgr(kwargs, "cmap", "bwr")
        annot = kwarg_mgr(kwargs, "annot", True)
        fmt = kwarg_mgr(kwargs, "fmt", ".2f")
        label = kwarg_mgr(kwargs, "label", "factor loadings")
        
        sns.heatmap(
            self.loading,
            vmax=1, vmin=-1, cmap=cmap,
            annot=annot, fmt=fmt,
            ax=ax
        )

        ax.collections[0].colorbar.set_label(label);
    

    def plot_var(self, ax: None, condition: tuple = None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
            
            fig_id = (
                "plot_var",
                list(self.communalities.keys())[-1]
            ) if condition is None else ("plot_var", condition)
            
            self.figs = {**self.figs, fig_id: (fig, ax)}
        
        if (condition is not None) and (condition not in self.communalities):
            self.rotate(*list(condition))
            
        assert issubclass(type(ax), mpl.axes.SubplotBase),\
        f"ax expected matplotlib.axes._subplots.AxesSubplot"
    
        color = kwarg_mgr(kwargs, "color", ["r", ".8"])
        ec = kwarg_mgr(kwargs, "annot", [".3", ".3"])
        orient = kwarg_mgr(kwargs, "orient", "h")
        
        sns.barplot(
            x=self.communality,
            y=self.data.columns,
            orient=orient, color=color[0],
            label="communalities",
            edgecolor=ec[0]
        )

        sns.barplot(
            x=self.communality + self.uniqueness,
            y=self.data.columns,
            orient=orient, color=color[1],
            zorder=0, label="uniquenesses",
            edgecolor=ec[1]
        )

        ax.legend();


class AutoFactorAnalyzer():
    def __init__(
        self,
        data: pd.core.frame.DataFrame,
        rotation: str,
        random_state: int = 0,
        nullmodel: str = "norm",
        use_smc: bool = False,
        method: str = "minres",
        threshold: float = 0.5,
        use_kaiser_criterion: bool = False,
    ):

        n_f = parallel_analysis(
            data,
            random_state=random_state,
            nullmodel=nullmodel,
            use_smc=use_smc
        )

        model = FactorAnalyzer(
            n_factor=n_f, data=data
        )

        model.rotate(
            rotation=rotation,
            method=method,
            n_factor=n_f
        )

        while not np.all(
            np.abs(model.loading).max(axis=0) > threshold
        ):
            if use_kaiser_criterion:
                n_f = len(
                    model.model.get_eigenvalues()[1][
                        model.model.get_eigenvalues()[1] > 1
                    ]
                )
            
            else:
                n_f = parallel_analysis(
                    model.result,
                    random_state=random_state,
                    nullmodel=nullmodel,
                    use_smc=use_smc
                )
            
            model.rotate(
                rotation=rotation,
                method=method,
                n_factor=n_f
            )
        
        self._model_ = model
        self.model = model.model
        self.models = model.models
        self.result = model.result
        self.results = model.results
        self.communality = model.communality
        self.communalities = model.communalities
        self.uniqueness = model.uniqueness
        self.uniquenesses = model.uniquenesses
        self.data = model.data
        self.n_f = model.n_f
        self.figs = model.figs
        
        return model
    

    def overwrite(self):
        self.model = self._model_.model
        self.models = self._model_.models
        self.result = self._model_.result
        self.results = self._model_.results
        self.communality = self._model_.communality
        self.communalities = self._model_.communalities
        self.uniqueness = self._model_.uniqueness
        self.uniquenesses = self._model_.uniquenesses
        self.data = self._model_.data
        self.n_f = self._model_.n_f
        self.figs = self._model_.figs

    def rotate(
        self,
        rotation: str,
        method: str = "minres",
        n_factor: int = None
    ):
        self._model_.rotate(
            rotation=rotation,
            method = method,
            n_factor=n_factor
        )
        
        self.overwrite()

    def heatmap(self, ax: None, condition: tuple = None, **kwargs):
        self._model_.heatmap(
            ax=ax,
            condition=condition,
            **kwargs
        )

        self.overwrite()
    

    def plot_var(self, ax: None, condition: tuple = None, **kwargs):
        self._model_.plot_var(
            ax=ax,
            condition=condition,
            **kwargs
        )

        self.overwrite()
