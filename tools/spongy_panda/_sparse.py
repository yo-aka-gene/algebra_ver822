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
        
    def __call__(self):
        return sp.csc_matrix(self.values)
    
    def __len__(self):
        return self.values.__len__()
    
    def __add__(self, other):
        return SparseDF(
            self.values.__add__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __sub__(self, other):
        return SparseDF(
            self.values.__sub__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __mul__(self, other):
        return SparseDF(
            self.values.__mul__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __matmul__(self, other):
        return SparseDF(
            self.values.__matmul__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __pow__(self, other):
        return SparseDF(
            self.values.__pow__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __truediv__(self, other):
        return SparseDF(
            self.values.__truediv__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __floordiv__(self, other):
        return SparseDF(
            self.values.__floordiv__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __mod__(self, other):
        return SparseDF(
            self.values.__mod__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __pos__(self):
        return SparseDF(
            self.values.__pos__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __neg__(self):
        return SparseDF(
            self.values.__neg__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __eq__(self, other):
        return SparseDF(
            self.values.__eq__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __lt__(self, other):
        return SparseDF(
            self.values.__lt__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __le__(self, other):
        return SparseDF(
            self.values.__le__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __ne__(self, other):
        return SparseDF(
            self.values.__ne__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __ge__(self, other):
        return SparseDF(
            self.values.__ge__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __gt__(self, other):
        return SparseDF(
            self.values.__gt__(other), 
            index=self.index,
            columns=self.columns
        )
    
    def __cmp__(self, other):
        return SparseDF(
            self.values.__cmp__(other), 
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