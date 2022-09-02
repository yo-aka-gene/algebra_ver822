from typing import Union

import numpy as np
import matplotlib.pyplot as plt

fig_kwarg = {
    'facecolor': 'white',
    'dpi': 300,
    'bbox_inches': 'tight',
    'pad_inches': 0.05
}
    
cns_markers = [
    # astrocyte markers
    "GFAP", "S100B", "SLC1A3", "ALDH1L1",
    
    # oligodendrocyte markers
    "OLIG2", "CNP", "CA2",
    
    # neuron markers
    'TUBB3', 'RBFOX3', 'NEUROD1', 'DCX', 'SOX5',
    'FOXO3',
    
    # neural progenitor markers
    'CALB1', 'CALB2', 'CD24', 'PROX1', 'NCAM1',
    'REST','CSPG4', 'PAX6', 'ASCL1', 
    
    # neural stem markers
    'FABP7', 'NES', 'EOMES', 'MCM2', 
    
    # vascular marker
    'PECAM1'
]

cns_subtype_markers = [
    # glia markers
    "SLC1A2", "VIM", "AQP4",

    # astrocyte markers
    "GFAP", "S100B", "SLC1A3", "ALDH1L1", "BYSL",
    "GJA1", "GLUL", "PYGB", 

    # oligodendrocyte markers
    "OLIG2", "CNP", "CA2", "NFIA", "NFIB",

    # OPC mekrers
    "PDGFRA",

    # Microglial markers
    "SPP1",

    # neuron markers
    "RBFOX3", "TBR1", "SOX5", "DCX", "BCL11B",
    "FEZF2", "SATB2", "CUX1", "CUX2", "POU3F2",
    "POU3F3", "TUBB3",

    # CGE markers
    "CALB2", "RELN", "NR2F1", "NR2F2", "VIP",

    # Excitatory Neuron markers
    "NEUROD1", "NEUROD2", "NEUROD6", "GRIN2B",

    # LGE markers
    "SIX3", "ISL1", "EBF1",

    # MGE
    "LHX6", "LHX8", "MAF", "SST", "ERBB4",
    "SOX6", "NKX2-1",

    # NPC
    "HES1", "HES5", "TYMS", "FABP7", "EOMES",
    "NEUROG1", "NCAM1", "TTYH1", "DLX2", "GAD2",
    "ASCL1", "NEUROG1", "NEUROG2", "PROX1", "TOP2A",
    "NUSAP1", "NR2E1", "CD24",

    # NSC
    "PAX6", "NES", "SOX1", "SOX2", "MCM2",
    "PCNA", "MKI67", "FOXO3", "BHLHE22", "REST",

    # RBC
    "HBB", "HBM", "HBA1", "HBA2",

    # Endothelial
    "IGFBP7", "PECAM1",

    # Epithelial
    "KRT14", "KRT16", "KRT17", "EMX2"

]

def custom_bwr(
    array: np.ndarray,
    centor: Union[int, float],
    saturation: float,
    reverse: bool = False
):
    assert 0 < saturation <= 1, \
        f"saturation expected within range [0, 1], got {saturation}"
    
    ret = [
        plt.cm.bwr(0.5 + (v - centor) / np.abs((v - centor))**(saturation)) for v in array
    ]
    
    if reverse:
        ret = [
        plt.cm.bwr(0.5 + (centor - v) / np.abs((centor - v))**(saturation)) for v in array
    ]
    
    return ret
    