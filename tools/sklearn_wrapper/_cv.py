from typing import Union

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedGroupKFold as SGKF


class StratifiedGroupKFold():
    def __init__(
        self,
        n_splits: int = 5,
        shuffle: bool = False,
        random_state: int = None,
        groups: Union[np.ndarray, pd.core.frame.DataFrame, pd.core.series.Series] = None
    ):
        self.model = SGKF(n_splits=n_splits,shuffle=shuffle, random_state=random_state)
        self.groups = groups


    def get_n_splits(self, X: object = None, y: object = None, groups: object = None):
        return self.model.get_n_splits(X, y, groups)


    def split(
        self,
        X: Union[np.ndarray, pd.core.frame.DataFrame],
        y: Union[np.ndarray, pd.core.frame.DataFrame, pd.core.series.Series] = None,
        groups: Union[np.ndarray, pd.core.frame.DataFrame, pd.core.series.Series] = None
    ):
        return self.model.split(X, y, self.groups)
