from typing import Any, Callable, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hdbscan import HDBSCAN
from joblib import Memory
from scipy.sparse import csr_matrix
from sklearn.metrics import silhouette_score

from ._hyperparams import _Gridsearch

class HdbscanClustering():
    def __init__(
        self,
        min_cluster_size: int = 5,
        min_samples: int = None,
        metric: Union[Callable, str] = 'euclidean',
        p: int = None,
        alpha: float = 1.0,
        cluster_selection_epsilon: float = 0.0,
        algorithm: str = "best",
        leaf_size: int = 40,
        memory: Union[Memory, str] = Memory(location=None),
        approx_min_span_tree: bool = True,
        gen_min_span_tree: bool = False,
        core_dist_n_jobs: int = 4,
        cluster_selection_method: str = "eom",
        allow_single_cluster: bool = False,
        prediction_data: bool = False,
        match_refernce_implementation: bool = False,
        init: bool = True,
        **kwargs
    ):
        if init:
            self.model = HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                metric=metric,
                p=p,
                alpha=alpha,
                cluster_selection_epsilon=cluster_selection_epsilon,
                algorithm=algorithm,
                leaf_size=leaf_size,
                memory=memory,
                approx_min_span_tree=approx_min_span_tree,
                gen_min_span_tree=gen_min_span_tree,
                core_dist_n_jobs=core_dist_n_jobs,
                cluster_selection_method=cluster_selection_method,
                allow_single_cluster=allow_single_cluster,
                prediction_data=prediction_data,
                match_refernce_implementation=match_refernce_implementation,
                **kwargs
            )
        
        self.searcher = None

    def fit(
        self,
        X: Union[csr_matrix, pd.core.frame.DataFrame, np.ndarray],
    ):
        return self.model.fit(X)


    def fit_predict(
        self,
        X: Union[csr_matrix, pd.core.frame.DataFrame, np.ndarray],
        y: np.ndarray
    ):
        return self.model.fit_predict(X, y)


    def generate_prediction_data(self):
        return self.model.generate_prediction_data()


    def grid_search(
        self,
        search_params: Dict[str, np.ndarray],
        X: pd.core.frame.DataFrame,
        objective: Callable = silhouette_score,
        fix_params: Dict[str, Any] = {},
        ignore_unclassifiable: bool = True,
        n_jobs: int = None,
        backend: Union[str, Any] = "loky",
        prefer: str = None,
        require: str = None,
        verbose: int = 0,
        timeout: float = None,
        pre_dispatch: Union[str, int] = "2 * n_jobs",
        batch_size: Union[int, str] = "auto",
        temp_folder: str = None,
        max_nbytes: Union[int, str] = "1M",
        mmap_mode: str = "r",
        eval_func: Callable = np.argmax
    ):
        if self.searcher is None:
            self.data = X
            self.searcher = _Gridsearch(
                search_params=search_params,
                X=X,
                objective=objective,
                fix_params=fix_params,
                ignore_unclassifiable=ignore_unclassifiable,
                n_jobs=n_jobs,
                backend=backend,
                prefer=prefer,
                require=require,
                verbose=verbose,
                timeout=timeout,
                pre_dispatch=pre_dispatch,
                batch_size=batch_size,
                temp_folder=temp_folder,
                max_nbytes=max_nbytes,
                mmap_mode=mmap_mode
            )
    
        return self.searcher.best_params(eval_func=eval_func)
    

    def tune_fit(self, eval_func: Callable = np.argmax):
        self.model = HDBSCAN(**self.searcher.best_params(eval_func))
        return self.fit(self.data)
    

    def tune_fit_predict(
        self,
        y: np.ndarray,
        eval_func: Callable = np.argmax
    ):
        self.model = HDBSCAN(**self.searcher.best_params(eval_func))
        return self.fit_predict(self.data, y)
