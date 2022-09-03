import json

def read_json(path: str, from_r: bool = False) -> dict:
    with open(path, "r") as f:
        ret = json.load(f)
    
    if from_r:
        lst = [
            (
                v.split(": ")[0][1:-1],
                list(
                    w[1:-1] for w in v.split(": ")[1][
                        1:(lambda x: -3 if x == len(ret[0].split("\n  ")[1:]) else -2)(i + 1)
                    ].split(", ")
                )
            )
            for i, v in enumerate(ret[0].split("\n  ")[1:])
        ]

        ret = {v[0]: v[1] for v in lst}
    
    return ret
