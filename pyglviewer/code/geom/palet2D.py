#!/usr/bin/env python3

import numpy as np



## Class defining a 2D palet
class Palet2D(object):

    def __init__(self, position, velocity, radius = 1, theta = 0, omega = 0, nbVertices = 40, color = np.array((0.1, 0.1, 0.1), np.float32), mass=18, visible=True):
        ## Constructor
        # @param positions  1-D Numpy array concatenating
        #                   the 2D positions of the mesh
        # @param indices    1-D Numpy array for the triangles indices (triplets)
        # @param colours    1-D Numpy array concatenating the
        #                   vertices colours (r, g, b)

        # Register the data and make sure of the types
        self.position = position
        self.radius = radius
        self.nbVertices = nbVertices
        self.velocity = velocity
        self.color = color
        self.mass = mass
        self.mInertie = 0.5*self.mass*self.radius*self.radius
        self.theta = theta
        self.omega = omega
        self.visible = visible

        if (color.size != 3):
            raise Exception("Wrong buffer size")

        # Fields to lighten the redraw
        self.positionsUpdated = True
        self.coloursUpdated = True


    def __getattribute__(self, name):
        ## Attribute accessor
        # Overload it to have "const" accessor to the posititions and colours :
        # * positions or colours : "non-const" accessor,
        #     trigger the buffer update
        # * constPositions or constColours  : "const" accessor,
        #     does not trigger the buffer update

        if (name == "position"):
            self.positionsUpdated = True
        elif (name == "constPosition"):
            return object.__getattribute__(self, "position")
        elif ((name == "colour") or (name == "color")):
            self.coloursUpdated = True
            return object.__getattribute__(self, "colour")
        elif ((name == "constColour") or (name == "constColor")):
            return object.__getattribute__(self, "colour")
        
        return object.__getattribute__(self, name)
            


    def __setattr__(self, name, value):
        ## Attribute setter
        # Overload it to have "const" accessor to the posititions and colours :
        # * positions or colours : "non-const" accessor,
        #     set trigger the buffer update
        # * constPositions or constColours  : "const" accessor,
        #     fails

        if (name == "position"):
            self.positionsUpdated = True
        elif ((name == "colour") or (name == "color")):
            self.coloursUpdated = True
            object.__setattr__(self, "colour", value)
            return
        elif (name == "constPosition") or \
             (name == "constColour") or (name == "constColor"):
            raise Exception("Tried to set a const field")
        
        object.__setattr__(self, name, value)
