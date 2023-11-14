#!/usr/bin/env python3

import numpy as np

import matplotlib.pyplot as plt

from .abstract_dynamic_system import AbstractDynamicSystem

class Pendule(AbstractDynamicSystem):

    def __init__(self, mesh):
        ## Constructor
        # @param self
        # @param mesh
        super().__init__()
        self.mesh = mesh

        # Parameters
        self.h = .005
        initial_theta_0 = np.pi/5
        initial_theta_1 = 0.
        self.theta = np.array([initial_theta_0, initial_theta_1])
        self.g = 9.81
        self.l = 1.
        self.m = 1.

        self.time = 0.
        self.time_plot = []
        self.theta_plot = []
        self.potential_plot = []
        self.cinetic_plot = []
        self.mecanic_plot = []

    # def step(self):
    #     """
    #         Euler explicite
    #     """
    #
    #     #Compute the step n+1
    #     self.theta[0] = self.theta[0] + self.h * self.theta[1]
    #     self.theta[1] = self.theta[1] + self.h * (-(self.g / self.l) * self.theta[0])
    #
    #     self.mesh.positions[2] = self.l*np.sin(self.theta[0])
    #     self.mesh.positions[3] = -self.l*np.cos(self.theta[0])
    #
    #     self.time += self.h
    #     self.time_plot.append(self.time)
    #     self.theta_plot.append(self.theta[0])
    #     self.cinetic_plot.append(self.m*(self.l*self.theta[1]*self.l*self.theta[1])/2)
    #     self.potential_plot.append(self.m*self.g*(self.l-self.l*np.cos(self.theta[0])))
    #     self.mecanic_plot.append(self.cinetic_plot[-1]+self.potential_plot[-1])

    def step(self):
        """
        Euler semi-implicite
        """

        #Compute the step n+1
        self.theta[1] = self.theta[1] + self.h*(-(self.g/self.l)*self.theta[0])
        self.theta[0] = self.theta[0] + self.h*self.theta[1]

        self.mesh.positions[2] = self.l*np.sin(self.theta[0])
        self.mesh.positions[3] = -self.l*np.cos(self.theta[0])

        self.time += self.h
        self.time_plot.append(self.time)
        self.theta_plot.append(self.theta[0])
        self.cinetic_plot.append(self.m*(self.l*self.theta[1]*self.l*self.theta[1])/2)
        self.potential_plot.append(self.m*self.g*(self.l-self.l*np.cos(self.theta[0])))
        self.mecanic_plot.append(self.cinetic_plot[-1]+self.potential_plot[-1])

        # def step(self):
        #     """
        #     Euler implicite
        #     """
        #
        #     #Méthode de Newton

    def plot(self):

        plt.plot(self.time_plot, self.theta_plot, label='Theta')
        plt.plot(self.time_plot, self.potential_plot, label='Potential')
        plt.plot(self.time_plot, self.cinetic_plot, label='Cinétique')
        plt.plot(self.time_plot, self.mecanic_plot, label='Mécanique')

        plt.show()

