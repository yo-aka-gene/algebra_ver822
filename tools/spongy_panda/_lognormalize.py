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


def change_base(
    sdf: SparseDF,
    current_base: Union[int, float],
    target_base: Union[int, float]
) -> SparseDF:
    base = "" if target_base == np.e else target_base
    data = sdf().tocoo()
    data.data = eval(f"np.log{base}")(current_base ** data.data)
    return SparseDF(data=data.tocsc(), index=sdf.index, columns=sdf.columns)


def reverse_transform(
    sdf: SparseDF,
    current_base: Union[int, float]
) -> SparseDF:
    data = sdf().tocoo()
    data.data = current_base ** data.data
    return SparseDF(data=data.tocsc(), index=sdf.index, columns=sdf.columns)
