#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# This file is part of SimulationTeachingElan, a python code used for teaching at Elan Inria.
#
# Copyright 2020 Mickael Ly <mickael.ly@inria.fr> (Elan / Inria - Université Grenoble Alpes)
# SimulationTeachingElan is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# SimulationTeachingElan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SimulationTeachingElan.  If not, see <http://www.gnu.org/licenses/>.
#

import glfw
import numpy as np

from scenes import *
from viewer import Viewer


# Separates the main to make sure all objects are deleted
# before glfw.terminate is called
def main():
    viewer=Viewer(height=960, width=1280,
                  bgColor = np.array([0.4, 0.4, 0.4]))

    # Loading the scene
    
    #indexedTest(viewer)
    #dynamicTest(viewer)
    #rodTest(viewer)
    system=pendule(viewer)

    # Main loop
    viewer.run()

    system.plot()



if __name__ == '__main__':

    # Initialization
    glfw.init()

    main()
    
    # End
    glfw.terminate()
