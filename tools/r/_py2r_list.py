def list2tsv(array: list, filename: str):
    with open(filename, "w") as f:
        f.write("\n".join(list(array)))
        f.close()
