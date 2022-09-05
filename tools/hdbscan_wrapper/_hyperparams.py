from itertools import combinations
from typing import Any, Callable, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hdbscan import HDBSCAN
from joblib import delayed, Parallel
from sklearn.metrics import silhouette_score
from tqdm.notebook import tqdm

class _Gridsearch():
    def __init__(
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
        mmap_mode: str = "r"
    ):
        self.search_params = search_params.keys()
        self.space = np.meshgrid(
            *[v for v in search_params.values()]
        )

        def func(dat, obj, params, ignore):
            _cluster = pd.Series(
                HDBSCAN(**params).fit_predict(dat),
                index = dat.index,
                name = "cluster"
            )
            dat = dat.loc[_cluster[_cluster != -1].index, :] if ignore else dat
            _cluster = _cluster[_cluster != -1] if ignore else _cluster

            return obj(dat, _cluster)
        
        def functional(dat, obj, ignore):
            return lambda params: func(dat, obj, params, ignore)


        self.result = np.array(
            Parallel(
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
            )(
                delayed(
                    lambda tuple_of_array: list(
                        map(
                            functional(
                                X,
                                objective,
                                ignore_unclassifiable
                            ),
                            [
                                {
                                    **{k: int(v) for k, v in zip(
                                        search_params.keys(),
                                        list(tup_val)
                                    )},
                                    **fix_params
                                } for tup_val in zip(*tuple_of_array)
                            ]
                        )
                    )
                )(tup_array) for tup_array in tqdm(
                    zip(*self.space),
                    desc="Grid Search",
                    total=len(self.space[0])
                )
            )
        )

        return None
    

    def best_params(self, eval_func: Callable = np.argmax):
        lst_pos = [sp.ravel()[eval_func(self.result)] for sp in self.space]
        return {
            key: v for key, v in zip(
                self.search_params,
                lst_pos
            )
        }

    def coordinate(self):
        dict_vals = {
            k: v for k, v in zip(
                self.search_params,
                self.space
            )
        }
        return [
            dict(
                x=x,
                y=y,
                args=dict(
                    x=dict_vals[x],
                    y=dict_vals[y],
                    z=self.result
                )
            ) for x, y in combinations(self.search_params, 2)
        ]
