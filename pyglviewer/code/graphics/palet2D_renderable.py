#!/usr/bin/env python3

import OpenGL.GL as GL
import numpy as np
import numpy.linalg as la

from geom.mesh2D import Mesh2D
from .mesh2D_renderable import Mesh2DRenderable


## Class rendering a rod
## Creates a mesh from a rod to give it a thickness
## as GLLineWidth has no standard implementation...
class Palet2DRenderable(Mesh2DRenderable):


    def __init__(self, palet, scale=1):
        ## Constructor
        # Generates a mesh around the rod and
        # initialized the GPU buffers
        # @param self
        # @param palet

        self.glId = None
        
        self.palet = palet
        self.scale = scale

        # Init mesh
        nbVerticesPalet = palet.nbVertices
        nbVerticesMesh = nbVerticesPalet + 1
        positions = np.zeros(2 * nbVerticesMesh, np.float64)
        colours = np.zeros(3 * nbVerticesMesh, np.float32)
        for i in range(nbVerticesPalet//8 + 1, 3*nbVerticesPalet//8 + 1):
            colours[3*i:3*i + 3] = palet.color

        ## Indices
        indices = []
        for i in range(1, nbVerticesMesh):
            indices.append(0)
            indices.append(i)
            indices.append(i+1)
        indices[-1] = 1
        

        self.mesh = Mesh2D(positions, indices, colours)
        self.updateMeshPositions()
        self.updateMeshColours()
        
        # Init mesh renderable
        super().__init__(self.mesh)



    def updateMeshPositions(self):
        ## Compute the positions of the mesh from the palet
        # Just adds a thickness to the segments
        # @param self

        # Data
        position = self.palet.position
        meshPositions = self.mesh.positions

        meshPositions[0] = self.scale * position[0]
        meshPositions[1] = self.scale * position[1]
        for vId in range(self.palet.nbVertices):
            meshPositions[2*(vId+1)] = self.scale * (position[0] + np.cos(2*vId*np.pi/self.palet.nbVertices + self.palet.theta)*self.palet.radius)
            meshPositions[2*(vId+1)+1] = self.scale * (position[1] + np.sin(2*vId*np.pi/self.palet.nbVertices + self.palet.theta)*self.palet.radius)
            
    def updateMeshColours(self):
        ## Compute the colours of the mesh from the rod
        #  Simply report them to the extended vertices
        # @param self

        # Data
        pass
        # colours = self.palet.color
        # meshColours = self.mesh.colours

        # for vId in range(self.palet.nbVertices):
        #     meshColours[3*vId:3*(vId+1)] = colours
        # meshColours[3*(vId+1):3*(vId+2)] = (1, 0, 0)
                
        
            
    def updatePositionsBuffer(self):
        ## Update the GPU colour buffer
        # @param self
        if (not self.palet.positionsUpdated):
            return
        self.updateMeshPositions()
        super().updatePositionsBuffer()
        self.palet.positionsUpdated = self.mesh.positionsUpdated
        
    def updateColoursBuffer(self):
        ## Update the GPU colour buffer
        # @param self
        if (not self.palet.coloursUpdated):
            return
        self.updateMeshColours()
        super().updateColoursBuffer()
        self.palet.coloursUpdated = self.mesh.coloursUpdated
        

    def draw(self, *args):
        if self.palet.visible:
            super().draw(*args)
        return

    def __del__(self):
        super().__del__()
