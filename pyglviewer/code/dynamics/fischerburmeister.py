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


import math
import numpy as np

class FischerBurmeister:
    def __init__(self, nContacts, mu, A, b, rho = 1):
        self.nContacts = nContacts
        self.mu = mu
        self.A = A
        self.b = b
        self.rho = rho

    def compute(self, x):
        d = x.shape[0] // self.nContacts
        fb = np.empty( x.shape )

        # compute FB(x, Ax+b) for all contacts
        y = self.A.dot(x) + self.b
 
        for i in range(self.nContacts):
            id = i*d
            fb[id:id+d] = self.computeLocalFBMuDeSaxce(self.mu, x[id:id+d], y[id:id+d])

        return fb

    def computeJacobian(self, x):
        d = x.shape[0] // self.nContacts
        dfb_dx = np.zeros( (x.shape[0], x.shape[0]) )
        
        # compute FB Jacobian(x, Ax+b) for all contacts
        y = self.A.dot(x) + self.b

        for i in range(self.nContacts):
            id = i*d
            xi = x[id:id+d]
            yi = y[id:id+d]
            dfbi_dx, dfbi_dy = self.computeLocalFBJacobianMuDeSaxce(self.mu, xi, yi)
            dfb_dx[id:id+d, id:id+d] = dfbi_dx
            for j in range(self.nContacts):
                jd = j*d
                dfb_dx[id:id+d,jd:jd+d] += dfbi_dy.dot(self.A[id:id+d,jd:jd+d])
    
        return dfb_dx
        
    def computeLocalFBMuDeSaxce(self, mu, x, y):
        d = x.shape[0]
        yt = y.copy()
        yt[0] += mu * np.linalg.norm(y[1:d])
        return self.computeLocalFBMu(mu, x, yt)

    def computeLocalFBJacobianMuDeSaxce(self, mu, x, y):
        d = x.shape[0]
        yt = y.copy()
        s = np.linalg.norm(y[1:d])
        yt[0] += mu * s

        dfb_dx, dfb_dy = self.computeLocalFBJacobianMu(mu, x, yt)
        if s != 0.:
            dfb_dy[:,1:] += mu/s * np.outer(dfb_dy[:,0], y[1:d])

        return dfb_dx, dfb_dy

    def computeLocalFBMu(self, mu, x, y):
        d = x.shape[0]

        if self.mu == 0.:
            fb = np.empty(d)
            fb[0] = x[0] + y[0] - math.sqrt( x[0]**2 + y[0]**2 )
            fb[1:d] = x[1:d]
            return fb
        else:
            xh = x.copy()
            yh = y.copy()
            xh[0] *= mu
            yh[1:d] *= mu
            return self.computeLocalFB(xh, yh)

    def computeLocalFBJacobianMu(self, mu, x, y):
        d = x.shape[0]
        if self.mu == 0.:
            z = math.sqrt( x[0]**2 + y[0]**2 )
            dfb_dx = np.eye(d)
            dfb_dy = np.zeros((d,d))
            if z != 0.:
                dfb_dx[0,0] = 1. - x[0] / z
                dfb_dy[0,0] = 1. - y[0] / z
            return dfb_dx, dfb_dy
        else:
            xh = x.copy()
            yh = y.copy()
            xh[0] *= mu
            yh[1:d] *= mu
            dfb_dx, dfb_dy = self.computeLocalFBJacobian(xh, yh)
            dfb_dx[:,0] *= mu
            dfb_dy[:,1:] *= mu
            return dfb_dx, dfb_dy
        

    def computeLocalFB(self, x, y):
        d = x.shape[0]
        z2 = np.empty( d )
        z2[0] = np.inner(x,x) + np.inner(y,y)
        z2[1:d] =  x[0] * x[1:d] + y[0] * y[1:d]
        nz2t = np.inner(z2[1:d],z2[1:d])

        w1 = np.zeros( d )
        w2 = np.zeros( d )
        w1[0] = 0.5
        w2[0] = 0.5

        if nz2t == 0.:
            w1[1] = -0.5
            w2[1] = 0.5
        else:
            w1[1:d] = -0.5 * z2[1:d] / nz2t
            w2[1:d] = 0.5 * z2[1:d] / nz2t

        rlambda1 = math.sqrt( max( 0., z2[0] - 2.*nz2t ) )
        rlambda2 = math.sqrt( z2[0] + 2.*nz2t )

        z = rlambda1 * w1 + rlambda2 * w2
        fb = x + y - z

        return fb

    def computeLocalFBJacobian(self, x, y):
        d = x.shape[0]
        z2 = np.empty( d )
        z2[0] = np.inner(x,x) + np.inner(y,y)
        z2[1:d] =  x[0] * x[1:d] + y[0] * y[1:d]
        nz2t = np.inner(z2[1:d],z2[1:d])

        w1 = np.zeros( d )
        w2 = np.zeros( d )
        w1[0] = 0.5
        w2[0] = 0.5

        if nz2t == 0.:
            w1[1] = -0.5
            w2[1] = 0.5
        else:
            w1[1:d] = -0.5 * z2[1:d] / nz2t
            w2[1:d] = 0.5 * z2[1:d] / nz2t

        rlambda1 = math.sqrt( max( 0., z2[0] - 2.*nz2t ) )
        rlambda2 = math.sqrt( z2[0] + 2.*nz2t )

        z = rlambda1 * w1 + rlambda2 * w2

        dfb_dx = np.zeros( (d,d) )
        dfb_dy = np.zeros( (d,d) )
        if rlambda2 == 0.:
            return dfb_dx, dfb_dy
        elif rlambda1 == 0.:
            izn = 1. / ( x[0]*x[0] + y[0]*y[0] )
            dfb_dx[:] = ( 1. - x[0]*izn ) * np.eye(d)
            dfb_dy[:] = ( 1. - y[0]*izn ) * np.eye(d)
        else:
            det = rlambda1 * rlambda2
                                
            invLz = np.empty( (d, d) )

            invLz[0, 0] = z[0]
            invLz[0,1:d] = -z[1:d]
            invLz[1:d,0] = -z[1:d]
            invLz[1:d,1:d] = ( det * np.eye(d-1) + np.outer(z[1:d], z[1:d]) ) / z[0]
            invLz /= det

            L = x[0] * np.eye(d)
            L[0,1:d] = x[1:d]
            L[1:d,0] = x[1:d]
            dfb_dx = np.eye(d)
            dfb_dx -= invLz.dot(L)

            L = y[0] * np.eye(d)
            L[0,1:d]= y[1:d]
            L[1:d,0] = y[1:d]
            dfb_dy[:] = np.eye(d)
            dfb_dy[:] -= invLz.dot(L)
        
        # print("x=", x)
        # print("y=", y)
        # print("rlambda1 = ", rlambda1)
        # print("rlambda2 = ", rlambda2)
        # print( "dfb_dx = ", dfb_dx )
        # print( "dfb_dy = ", dfb_dy )
        
        return dfb_dx, dfb_dy

