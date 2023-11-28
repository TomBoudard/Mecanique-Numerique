#!/usr/bin/env python3

import numpy as np

import matplotlib.pyplot as plt

from .abstract_dynamic_system import AbstractDynamicSystem

#Global Variable(s)
g = 9.81
alpha = 150

def eulerExplicite(X, h, m, force):
    a = np.array([force[0], force[1]])/m
    V = np.array([X[2], X[3], a[0], a[1]])
    X_1 = X + h*V

    return X_1

class CurlingDynamic(AbstractDynamicSystem):

    def __init__(self, palets):
        ## Constructor
        # @param self
        # @param mesh
        super().__init__()
        self.palets = palets

        # Parameters
        self.h = .02
        self.time = 0.

    def step(self):
    
        #Compute the step n+1 of each palet
        activePalets = [elt for elt in self.palets if elt.visible]

        for palet in activePalets:

            X = np.array([*palet.position, *palet.velocity])
            #Forces
            
            resForce = -alpha*palet.velocity

            #Step n+1
            X_1 = eulerExplicite(X, self.h, palet.mass, resForce)


            palet.position = X_1[:2]
            palet.velocity = X_1[2:]

            #FIXME To Change
            if np.linalg.norm(palet.velocity) < 0.1:
                palet.velocity *= 0

        for i, palet in enumerate(activePalets):
            for paletCompare in activePalets[i+1:]:
                if palet == paletCompare:
                    continue

                if np.linalg.norm(palet.position - paletCompare.position) < palet.radius + paletCompare.radius:
                    u1 = palet.velocity.copy()
                    u2 = paletCompare.velocity.copy()
                    dist = palet.position-paletCompare.position
                    palet.velocity -= 2*paletCompare.mass/(palet.mass+paletCompare.mass) * np.dot(u1-u2, dist)/np.dot(dist, dist)*dist
                    paletCompare.velocity -= 2*palet.mass/(palet.mass+paletCompare.mass) * np.dot(u2-u1, -dist)/np.dot(-dist, -dist)*(-dist)

        #Update next palet if possible
        allStopped = True
        for palet in activePalets:
            if np.linalg.norm(palet.velocity) != 0:
                allStopped = False
                break
        if allStopped:

            try:
                [elt for elt in self.palets if not(elt.visible)][0].visible = True
            except:
                pass

        self.time += self.h