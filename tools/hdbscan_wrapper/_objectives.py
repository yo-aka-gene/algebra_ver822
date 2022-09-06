import pandas as pd

def n_cluster(
    X: pd.core.frame.DataFrame,
    cluster: pd.core.series.Series
) -> int:
    X = None
    return len(cluster.unique())
