import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import LinearSegmentedColormap

mixture = tuple((3 * np.array(plt.cm.Purples(.6)) + 2.5 * np.array([0, 0, 1, 1])) / 5)

cdict = {
    "red": [
        (x, y0, y1) for x, y0, y1 in zip(
            np.linspace(0, 1, 256),
            np.linspace(plt.cm.Greys(.2)[0], mixture[0], 257)[:-1],
            np.linspace(plt.cm.Greys(.2)[0], mixture[0], 257)[1:],  
        )
    ],
    "green": [
        (x, y0, y1) for x, y0, y1 in zip(
            np.linspace(0, 1, 256),
            np.linspace(plt.cm.Greys(.2)[1], mixture[1], 257)[:-1],
            np.linspace(plt.cm.Greys(.2)[1], mixture[1], 257)[1:],  
        )
    ],
    "blue": [
        (x, y0, y1) for x, y0, y1 in zip(
            np.linspace(0, 1, 256),
            np.linspace(plt.cm.Greys(.2)[2], mixture[2], 257)[:-1],
            np.linspace(plt.cm.Greys(.2)[2], mixture[2], 257)[1:],  
        )
    ],
}

seurat = LinearSegmentedColormap("seurat", cdict, gamma=1)