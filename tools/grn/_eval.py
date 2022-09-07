import pandas as pd

def d_asterisk(
    subjective: pd.core.frame.DataFrame,
    objective: pd.core.frame.DataFrame
) -> float:
    s = subjective.loc[subjective.index, subjective.index]
    o = objective.loc[s.index, s.columns]
    return 1 - (((s.values == 1) * (o.values == 1)).sum() / (s.values == 1).sum())
