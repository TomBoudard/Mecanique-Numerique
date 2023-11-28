#!/usr/bin/env python3

import glfw
import numpy as np

from my_scenes import *
from viewer import Viewer


# Separates the main to make sure all objects are deleted
# before glfw.terminate is called
def main():
    viewer=Viewer(height=900, width=550,
                  bgColor = np.array([0.8745, 0.9961, 1]))

    # Loading the scene
    system=curling(viewer)

    # Main loop
    viewer.run()


if __name__ == '__main__':

    glfw.init()
    main()    
    glfw.terminate()
