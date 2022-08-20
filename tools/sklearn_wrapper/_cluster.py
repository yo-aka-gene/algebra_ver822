from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples
from tqdm.notebook import tqdm

from ._elbow import elbowplot
from ._entropy import entropy_plot
from ._silhouette import silhouette_curve, silhouette_plot

class KMeansClustering():
    def __init__(
        self,
        minimum: int,
        maximum: int,
        data: pd.core.frame.DataFrame,
        random_state:int = 0
    ):
        self.objects = {}
        self.results = [
            KMeans(n_clusters=i, random_state=random_state).fit(data) for i in tqdm(
                np.arange(minimum, maximum + 1),
                desc='k-means'
            )
        ]
        
        self.searchspace = [minimum, maximum]
        self.data = data
        self.sil_val = None
        self.sil_mean = None
        self.k = None
        self.k_idx = None


    def elbow(self, ax: None, show_diff: bool = False, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(2, 1) if show_diff else plt.subplots()
            self.objects = {**self.objects, "elbow": (fig, ax)}
            
        elbowplot(self.results, ax, self.searchspace, show_diff, **kwargs)
    

    def silhouette(self):
        if self.sil_val is None:
            self.sil_val = [
                silhouette_samples(
                    self.data,
                    v.labels_
                ) for v in tqdm(self.results, desc='Calculating Silhouette Coeff.')
            ]
            self.sil_mean = [
                v.mean() for v in tqdm(self.sil_val, desc='Calculating the mean')
            ]
    

    def silhouette_curve(self, ax: None, **kwargs):
        if self.sil_mean is None:
            self.silhouette()

        if ax is None:
            fig, ax = plt.subplots()
            self.objects = {**self.objects, "silhouette_curve": (fig, ax)}
        
        silhouette_curve(self.sil_mean, self.sil_val, ax, self.searchspace, **kwargs)


    def entropy_curve(self, ax: None, return_array: bool = False, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
            self.objects = {**self.objects, "entropy_curve": (fig, ax)}
            
        ret = entropy_plot(
            self.results,
            ax,
            self.searchspace,
            return_array,
            **kwargs
        )
        
        if return_array:
            selt.objects = {**self.objects, "entropy": ret}

    
    def optimal_k(self, manual: bool = False, k: int = 0):
        if manual:
            self.k = k
            self.k_idx = k - self.searchspace[0]
        
        else:
            if self.sil_mean is None:
                self.silhouette()
            
            self.k_idx = np.argmax(self.sil_mean)
            self.k = self.k_idx + self.searchspace[0]
        
        print(f"Optimal Number of k: {self.k}")
        
        self.label = pd.DataFrame(
            {
                'cluster': self.results[self.k_idx].labels_,
                'silhouette coeff.': self.sil_val[self.k_idx]
            },
            index = self.data.index
        )

        self.label_sep = [
            self.label.loc[
                self.label.cluster == i,
                :
            ].sort_values(
                "silhouette coeff.",
                ascending=True
            ) for i in np.sort(self.label.cluster.unique())
        ]
        
        self.cmap = lambda v: [
            eval("plt.cm." + v)(
                i / len(self.label_sep)
            ) for i in range(len(self.label_sep))
        ]
    
    def labels_(self, k:int = None):
        if k is None:
            if self.k_idx is None:
                self.optimal_k()
            
            return self.results[self.k_idx].labels_
    
    def silhouette_plot(
        self,
        ax: None,
        margin: float = 0.1,
        cmap: str = "jet",
        k: int = None,
        **kwargs
    ):
        if ax is None:
            fig, ax = plt.subplots()

        if self.sil_val is None:
            self.silhouette()
            if k is None:
                self.optimal_k()
        
        if isinstance(k, (int, np.int_)):
            assert (k >= self.searchspace[0]) and (k <= self.searchspace[1]), \
                f"k should be between {self.searchspace}, got {k}"
            
            if k != self.k:
                k_idx = k - self.searchspace[0]
                df = pd.DataFrame(
                    {
                        'cluster': self.results[k_idx].labels_,
                        'silhouette coeff.': self.sil_val[k_idx]
                    },
                    index = self.data.index
                )
                l_labels = [
                    df.loc[
                        df.cluster == i,
                        :
                    ].sort_values(
                        "silhouette coeff.",
                        ascending=True
                    ) for i in np.sort(df.cluster.unique())
                ]
                l_cmap = (lambda v: [
                    eval("plt.cm." + v)(
                        i / k
                    ) for i in range(k)
                ])(cmap)

            else:
                l_labels = self.label_sep
                l_cmap = self.cmap(cmap)

        else:
            l_labels = self.label_sep
            l_cmap = self.cmap(cmap)
            k = self.k
        
        silscore = self.sil_mean[k - self.searchspace[0]]
        
        axes = silhouette_plot(l_labels, ax, l_cmap, silscore, margin)
        ax.set(title="Silhouette Plot $(k=" + f"{k})$");
        
        self.objects = {
            **self.objects,
            "silhouette plot": (fig, *axes) if ax is None else axes
        }