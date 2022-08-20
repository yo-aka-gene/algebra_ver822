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
    