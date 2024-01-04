#!/usr/bin/env python3

import numpy as np

import matplotlib.pyplot as plt

from .fischerburmeister import FischerBurmeister
from .nsnewton import NonSmoothNewton

from .abstract_dynamic_system import AbstractDynamicSystem

#Global Variable(s)
g = 9.81
alpha = 6*np.pi*0.14 * 0.018e-3
mu = 0.1
mu2 = 0.1

def eulerExplicite(X, h, m, force):
    a = np.array([force[0], force[1]])/m
    V = np.array([X[2], X[3], a[0], a[1]])
    X_1 = X + h*V

    return X_1

def eulerSemiImplicite(X, h, m, force): # X = [x, y, theta, vx, vy, omega]
    Vel = X[2:]
    Vel_1 = Vel + h*(1/m)*force

    Pos = X[:2]
    Pos_1 = Pos + h * Vel_1

    X_1 = np.array([Pos_1[0], Pos_1[1], Vel_1[0], Vel_1[1]])
    
    return X_1

class CurlingDynamic(AbstractDynamicSystem):

    def __init__(self, palets):
        ## Constructor
        # @param self
        # @param mesh
        super().__init__()
        self.palets = palets

        # Parameters
        self.h = .005
        self.time = 0.

    def step(self):
        for _ in range(3):
            self.compute()

    def compute(self):

    
        #Compute the step n+1 of each palet
        activePalets = [elt for elt in self.palets if elt.visible]

        nbContacts = 0
        dico = {}
        pair = {i: [] for i, _ in enumerate(activePalets)}
        H = np.zeros((0, 0))
        M_1 = np.zeros((0, 0))
        V = np.zeros((0))
        #Liste de palets en contacts
        for i, palet in enumerate(activePalets):
            for jindex, paletCompare in enumerate(activePalets[i+1:]):
                j = jindex + i + 1
                if palet == paletCompare:
                    continue                
                if np.linalg.norm(palet.position - paletCompare.position) < palet.radius + paletCompare.radius:

                    pair[i].append((nbContacts, -1, j))
                    pair[j].append((nbContacts, 1, i))
                    nbContacts += 1

                    # Calc Pc
                    nc = (paletCompare.position - palet.position) / np.linalg.norm(paletCompare.position - palet.position)
                    tc = np.array([-nc[1], nc[0]])
                    Pc = np.array([nc, tc]).T
                    
                    partialJPalet = np.array((paletCompare.position + palet.position)/2 - palet.position) #All palets are the same size
                    partialJPalet[0], partialJPalet[1] = partialJPalet[1], -partialJPalet[0]

                    JPalet = np.array([
                        [1, 0, partialJPalet[0]],
                        [0, 1, partialJPalet[1]]
                    ])

                    HcPalet = np.matmul(Pc, JPalet)

                    partialJPaletCompare = np.array((palet.position + paletCompare.position)/2 - paletCompare.position) #All palets are the same size
                    partialJPaletCompare[0], partialJPaletCompare[1] = partialJPaletCompare[1], -partialJPaletCompare[0]

                    JPaletCompare = np.array([
                        [1, 0, partialJPaletCompare[0]],
                        [0, 1, partialJPaletCompare[1]]
                    ])

                    HcPaletCompare = np.matmul(Pc, JPaletCompare)

                    # H
                    H = np.vstack((H, np.zeros((2, H.shape[1]))))

                    if i not in dico.keys():
                        dico[i] = len(dico)
                        H = np.hstack((H, np.zeros((H.shape[0], 3))))
                        M_1 = np.hstack((M_1, np.zeros((M_1.shape[0], 3))))
                        M_1 = np.vstack((M_1, np.zeros((3, M_1.shape[1]))))
                        M_1[M_1.shape[0]-3: M_1.shape[0], M_1.shape[1]-3: M_1.shape[1]] = np.array([[1/palet.mass, 0, 0], [0, 1/palet.mass, 0], [0, 0, 1/palet.mInertie]])
                        V = np.hstack((V, np.array([palet.velocity[0], palet.velocity[1], palet.theta])))

                    Hindex = dico[i]
                    H[H.shape[0]-2:H.shape[0], 3*Hindex:3*Hindex+3] = -HcPalet

                    if j not in dico.keys():
                        dico[j] = len(dico)
                        H = np.hstack((H, np.zeros((H.shape[0], 3))))
                        M_1 = np.hstack((M_1, np.zeros((M_1.shape[0], 3))))
                        M_1 = np.vstack((M_1, np.zeros((3, M_1.shape[1]))))
                        M_1[M_1.shape[0]-3: M_1.shape[0], M_1.shape[1]-3: M_1.shape[1]] = np.array([[1/paletCompare.mass, 0, 0], [0, 1/paletCompare.mass, 0], [0, 0, 1/paletCompare.mInertie]])
                        V = np.hstack((V, np.array([paletCompare.velocity[0], paletCompare.velocity[1], paletCompare.theta])))

                    Hindex = dico[j]
                    H[H.shape[0]-2:H.shape[0], 3*Hindex:3*Hindex+3] = HcPaletCompare
        

        if (nbContacts != 0):
            A = np.matmul(np.matmul(H, M_1), np.transpose(H)) 

            b = np.matmul(H, V)

            # Create fischerburmeister object

            fish = FischerBurmeister(nbContacts, mu2, A, b)
            
            # Call nsnewton to solve and have a resForce

            tolerance = 1e-6
            maxIter = 100
            nnsm = NonSmoothNewton(tolerance, maxIter)

            x0 = np.zeros((2*nbContacts))
            xBest, _ = nnsm.solve(fish, x0)

            forceContact = H.T @ xBest

        for i, palet in enumerate(activePalets):

            X = np.array([*palet.position, *palet.velocity])
            #Forces
            resForce = np.zeros(2, dtype=np.float32)
            resCouple = 0

            for index, sens, paletCompareIndex in pair[i]:
                paletCompare = activePalets[paletCompareIndex]
                nc = (paletCompare.position - palet.position) / np.linalg.norm(paletCompare.position - palet.position)
                tc = np.array([-nc[1], nc[0]])
                Pc = np.array([nc, tc]).T
                print(" -------> xBest: ", xBest[2*index:2*index+2])
                print(" -------> forceContact: ", forceContact[2*index:2*index+2])

                resForce = resForce + np.matmul(np.linalg.inv(Pc), forceContact[2*index:2*index+2]) / self.h
                # print(" -------> resForce: ", resForce)

                # Frottement fluide
            resForce += -alpha*palet.velocity

                # Frottement solide
            paletNorm = np.linalg.norm(palet.velocity)
            resNorm = np.linalg.norm(resForce)
            #FIXME Integrer la vitesse pour utiliser la bonne velocity au bon pas de temps =)
            if (paletNorm != 0):
                # Vel + h*(1/m)*force
                # m/h > mu * m * g / paletNorm
                resForce += - min(mu * palet.mass * g * (1/paletNorm), palet.mass/self.h) * palet.velocity
            elif (resNorm != 0):
                resForce += - min(mu * palet.mass * g * (1/resNorm), 1) * resForce


            #Step n+1
            X_1 = eulerSemiImplicite(X, self.h, palet.mass, resForce)


            palet.position = X_1[:2]
            palet.velocity = X_1[2:]

            # #FIXME To Change
            # if np.linalg.norm(palet.velocity) < 0.1:
            #     palet.velocity *= 0

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