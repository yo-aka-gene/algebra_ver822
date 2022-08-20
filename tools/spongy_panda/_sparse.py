import json
import pickle

import numpy as np
import pandas as pd
import scipy.sparse as sp


class SparseDF():
    def __init__(
        self,
        data: sp._csc.csc_matrix,
        index: list,
        columns: list
    ):
        self.values = data.toarray()
        self.index = index
        self.columns = columns
        self.shape = self.values.shape
        self.T = self.values.T
        
    def __call__(self):
        return sp.csc_matrix(self.values)
    
    def t(self):
        return SparseDF(
            sp.csc_matrix(self.values.T),
            index = self.columns,
            columns = self.index
        )
    
    def sum(self, axis: int = None):
        if axis is None:
            return self.values.sum()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.sum(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
        
    def min(self, axis: int = None):
        if axis is None:
            return self.values.min()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.min(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def max(self, axis: int = None):
        if axis is None:
            return self.values.max()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.max(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def mean(self, axis: int = None):
        if axis is None:
            return self.values.mean()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.mean(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def prod(self, axis: int = None):
        if axis is None:
            return self.values.prod()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.prod(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def var(self, axis: int = None):
        if axis is None:
            return self.values.var()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.var(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def std(self, axis: int = None):
        if axis is None:
            return self.values.std()
        else:
            return SparseDF(
                sp.csc_matrix(self.values.std(axis=axis)),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def median(self, axis: int = 0):
        return SparseDF(
                sp.csc_matrix(self.to_df().median(axis=axis).values),
                index=None if axis == 0 else self.index,
                columns=None if axis == 1 else self.columns
            )
    
    def describe(self):
        df = self.to_df().describe()
        return SparseDF(
            sp.csc_matrix(df.values),
            index=df.index,
            columns=df.columns
        )
    
    def __len__(self):
        return self.values.__len__()
    
    def __add__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__add__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __sub__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__sub__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __mul__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__mul__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __matmul__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__matmul__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __pow__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__pow__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __truediv__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__truediv__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __floordiv__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__floordiv__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __mod__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__mod__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __radd__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__radd__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rsub__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rsub__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rmul__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rmul__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rmatmul__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rmatmul__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rpow__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rpow__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rtruediv__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rtruediv__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rfloordiv__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rfloordiv__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __rmod__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__rmod__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    
    def __pos__(self):
        return SparseDF(
            sp.csc_matrix(self.values.__pos__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __neg__(self):
        return SparseDF(
            sp.csc_matrix(self.values.__neg__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __eq__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__eq__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __lt__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__lt__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __le__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__le__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __ne__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__ne__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __ge__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__ge__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __gt__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__gt__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __cmp__(self, other):
        return SparseDF(
            sp.csc_matrix(self.values.__cmp__(other)), 
            index=self.index,
            columns=self.columns
        )
    
    def __contains__(self, other):
        return self.values.__contains__(other)
    
    def to_df(self):
        return pd.DataFrame(
            self.values,
            index=self.index,
            columns=self.columns
        )
    
    def __getitem__(self, item):
        return self.to_df()[item] if isinstance(item, str) else self.values[item]
    
    def to_pickle(self, filename: str):
        with open(filename, "wb") as f:
            pickle.dump(self, f)
    
    def save_npz(self, filename):
        sp.save_npz(filename, self())
    
    def to_csv(self, filename, index: bool = True):
        self.to_df().to_csv(filename, index=index)
        
    def colidx2json(self, filename):
        with open(filename, "w") as f:
            d = dict(
                index=self.index.to_list(),
                columns=self.columns.to_list()
            )
            json.dump(d, f)