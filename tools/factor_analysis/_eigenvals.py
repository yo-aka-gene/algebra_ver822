import factor_analyzer as fa
import numpy as np
import pandas as pd

def eigvals(
    data: pd.core.frame.DataFrame,
    use_smc: bool = False
) -> np.ndarray:
    cormat = np.corrcoef(data.T)
    
    if use_smc:
        np.fill_diagonal(cormat, fa.utils.smc(cormat))
    
    return np.sort(np.linalg.eigvals(cormat))[::-1]