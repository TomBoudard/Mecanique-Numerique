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
    for _ in range(20):
        palet = Palet2D(np.array([0.0, -1.5]), np.array([rd.uniform(-1, 1),rd.uniform(1, 2)]) * 3, radius=0.1, visible=False)
        paletRenderable = Palet2DRenderable(palet)
        viewer.addRenderable(paletRenderable)

        palets.append(palet)

    dyn = CurlingDynamic(palets)
    viewer.addDynamicSystem(dyn)

    return dyn