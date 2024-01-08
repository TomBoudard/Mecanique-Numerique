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
        velocity = np.array([0.0 + i * 0.01, 1.0]) * (1. + 0.2 * i) * 3
        palet = Palet2D(np.array([0.0, -10]), velocity, radius=0.2, color=color, visible=False)
        paletRenderable = Palet2DRenderable(palet, scale=0.19)
        viewer.addRenderable(paletRenderable)

        palets.append(palet)

    dyn = CurlingDynamic(palets)
    viewer.addDynamicSystem(dyn)

    disque4 = Palet2D(np.array([0.0, 5.0]), np.array([0.0, 0.0, 0.0]), radius=0.5, color=np.array([0, 0, 0]), visible=True)
    disque4Renderable = Palet2DRenderable(disque4, scale=0.19, mainColor=np.array((0.8745, 0.9961, 1)), isPalet=False)
    viewer.addRenderable(disque4Renderable)
    disque3 = Palet2D(np.array([0.0, 5.0]), np.array([0.0, 0.0, 0.0]), radius=1.5, color=np.array([0, 0, 0]), visible=True)
    disque3Renderable = Palet2DRenderable(disque3, scale=0.19, mainColor=np.array((1.0, 0.0, 0.0)), isPalet=False)
    viewer.addRenderable(disque3Renderable)
    disque2 = Palet2D(np.array([0.0, 5.0]), np.array([0.0, 0.0, 0.0]), radius=2.5, color=np.array([0, 0, 0]), visible=True)
    disque2Renderable = Palet2DRenderable(disque2, scale=0.19, mainColor=np.array((0.8745, 0.9961, 1)), isPalet=False)
    viewer.addRenderable(disque2Renderable)
    disque1 = Palet2D(np.array([0.0, 5.0]), np.array([0.0, 0.0, 0.0]), radius=3.5, color=np.array([0, 0, 0]), visible=True)
    disque1Renderable = Palet2DRenderable(disque1, scale=0.19, mainColor=np.array((0.0, 0.0, 1.0)), isPalet=False)
    viewer.addRenderable(disque1Renderable)

    return dyn