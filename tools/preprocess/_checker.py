from typing import List

def find(
    whole_list: List[str], 
    feat_list: List[str],
    show_missing: bool = True
) -> List[str]:
    ret = []
    for f in feat_list:
        ret += [f] if (f in whole_list) != show_missing else []
    return ret
