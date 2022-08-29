from typing import Union

import numpy as np

from ._sparse import SparseDF

def log2normalize(
    sdf: SparseDF
) -> SparseDF:
    data = sdf().tocoo()
    data.data = np.log2(data.data + 1)
    return SparseDF(data=data.tocsc(), index=sdf.index, columns=sdf.columns)


def log10normalize(
    sdf: SparseDF
) -> SparseDF:
    data = sdf().tocoo()
    data.data = np.log10(data.data + 1)
    return SparseDF(data=data.tocsc(), index=sdf.index, columns=sdf.columns)


def lnnormalize(
    sdf: SparseDF
) -> SparseDF:
    data = sdf().tocoo()
    data.data = np.log(data.data + 1)
    return SparseDF(data=data.tocsc(), index=sdf.index, columns=sdf.columns)
