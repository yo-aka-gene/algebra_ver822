import factor_analyzer as fa
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from tools.figure import kwarg_mgr

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
