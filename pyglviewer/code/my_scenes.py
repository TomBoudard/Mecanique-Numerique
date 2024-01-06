#!/usr/bin/env python3

import numpy as np
import random as rd

from dynamics.dummy_dynamic_system import *
from dynamics.curling_system import *
from graphics import *
from geom import *

def curling(viewer):
    """

    """

    palets = []
    for i in range(16):
        color = np.array([1, 0, 0]) if i % 2 else np.array([0, 0, 1])
        # velocity = np.array([rd.uniform(-0.15, 0.15), rd.uniform(1.7, 2.1)]) * 3
        velocity = np.array([0.5, 1.5]) * (1. + 0.2 * i) * 3
        palet = Palet2D(np.array([0.0, -10]), velocity, radius=0.15, color=color, visible=False)
        paletRenderable = Palet2DRenderable(palet, scale=0.19)
        viewer.addRenderable(paletRenderable)

        palets.append(palet)

    dyn = CurlingDynamic(palets)
    viewer.addDynamicSystem(dyn)

    return dyn