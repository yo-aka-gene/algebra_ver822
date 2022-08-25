# Schematic of this directory
---
```mermaid
graph LR
    A("GW9") =="preprocess"==> E("GW9")
    B("GW10")=="preprocess"==> F("GW10")
    C("GW11")=="preprocess"==> G("GW11")
    D("GW12")=="preprocess"==> H("GW12")

    E-."conjugated".->I(("GSE165388"))
    F-."conjugated".->I
    G-."conjugated".->I
    H-."conjugated".->I
```
## Overview of scripts
- GW9: [01_preprocess_GSE165388_gw9.md](./01_preprocess_GSE165388_gw9.md)
- GW10: [01_preprocess_GSE165388_gw10.md](./01_preprocess_GSE165388_gw10.md)
- GW11: [01_preprocess_GSE165388_gw11.md](./01_preprocess_GSE165388_gw11.md)
- GW12: [01_preprocess_GSE165388_gw12.md](./01_preprocess_GSE165388_gw12.md)

## Raw scripts
- GW9: [01_preprocess_GSE165388_gw9.Rmd](./01_preprocess_GSE165388_gw9.Rmd)
- GW10: [01_preprocess_GSE165388_gw10.Rmd](./01_preprocess_GSE165388_gw10.Rmd)
- GW11: [01_preprocess_GSE165388_gw11.Rmd](./01_preprocess_GSE165388_gw11.Rmd)
- GW12: [01_preprocess_GSE165388_gw12.Rmd](./01_preprocess_GSE165388_gw12.Rmd)


---
# Whole Picture of the Analysis
```mermaid
graph LR
I(("GSE165388"))-."Factor Analysis".->J("Cell Classes")-."GRN".->K{"comparison"}

    E("GW9")-."DEG".->L("annotation")
    F("GW10")-."DEG".->M("annotation")
    G("GW11")-."DEG".->N("annotation")
    H("GW12")-."DEG".->O("annotation")

    L-.->P{"comparison"}
    M-.->P
    N-.->P
    O-.->P

    K-."Calcualte Simirality".->Q(("m1_10x"))
    P-."check consistency".->Q
```
