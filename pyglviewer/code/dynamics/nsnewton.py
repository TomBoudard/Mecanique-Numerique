#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# This file is part of SimulationTeachingElan (pyDFCP), a python framework used for teaching at Elan Inria.
#
# Copyright 2020 Thibaut METIVET <thibaut.metivet@inria.fr> (Elan / Inria)
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

import numpy as np

class NonSmoothNewton:
    def __init__( self, tolerance, maxIter, verbose=False ):
        self.tolerance = tolerance
        self.maxIter = maxIter

        self.sigma2 = 1.e-4
        self.alpha = 0.5

        self.verbose = verbose

    def solve( self, F, x0 ):
        x = x0.copy()
        f = F.compute( x0 )
        phiInit = 0.5 * np.inner(f, f)
        
        if phiInit < self.tolerance:
            return x0, phiInit

        xBest = np.zeros( x0.shape )
        f = F.compute( xBest )
        phiZero = 0.5 * np.inner(f, f)
        phiBest = phiZero

        if phiZero < self.tolerance:
            return xBest, phiZero

        if phiZero < phiInit:
            phiBest = phiZero
            x[:] = xBest[:]
        else:
            phiBest = phiInit
            xBest[:] = x[:]

        for iter in range(self.maxIter):
            f = F.compute(x)
            dF_dx = F.computeJacobian(x)
            phi = 0.5 * np.inner(f, f)
            
            if phi < self.tolerance:
                phiBest = phi
                xBest[:] = x
                break
            if phi < phiBest:
                phiBest = phi
                xBest[:] = x[:]

            dphi_dx = dF_dx.transpose().dot(f)
            dx = np.linalg.solve( dF_dx, -f )
            if np.any( np.isnan( dx ) ):
                break
            
            proj = dx.dot(dphi_dx)

            if( proj > 0. or proj * proj < self.sigma2 * np.inner(dx, dx) * np.inner(dphi_dx, dphi_dx) ):
                dx *= self.alpha

            x += dx
            if self.verbose:
                print(f"x{iter} = ", x, "phi = ", phi)
        
        if self.verbose:
            print( "xBest =", xBest, "phiBest =", phiBest)
            print( "FBest =", F.compute(xBest) )

        return xBest, phiBest
