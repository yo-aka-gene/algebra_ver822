def list2tsv(array: list, filename: str):
    with open(filename, "w") as f:
        f.write("\n".join([v for v in array[:-1]] + [array[-1] + "\n"]))
        f.close()
