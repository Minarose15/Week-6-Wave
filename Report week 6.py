# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:23:22 2021

@author: Sabrina
"""

from numpy import exp, linspace, sqrt, copy, pi, zeros_like, abs
from matplotlib import pyplot as plt

class Wave:
    
    def __init__(self, L, sigma, speed):
        self.L = L #1m
        self.x0 = 0.4*self.L
        self.sigma = sigma #0.1
        self.v = speed
        self.dt = 0.0099/self.v  #auto-adjusts to avoid r > 1
        self.runTime = 0.5

    #This function is to help with the plotting of the Gaussian wave
    def Gauss(self, x):
        return exp(-(x-self.x0)**2/(2*self.sigma**2))/(self.sigma*sqrt(2*pi))
    
    #this function creates the movie. It also keeps track of the point of the original peak over time
    def propagateWave(self):
        self.t = 0
        steps = 100
        x, dx = linspace(0, self.L, steps, retstep = True)
        print(dx)
        y = self.Gauss(x)
        yNew = zeros_like(y)
        yOld = copy(y)
        con = self.v**2*self.dt**2/dx**2
        self.f = []
        lim = max(y)

        while(self.t<self.runTime):
            yNew[1:-1] = 2*y[1:-1]-yOld[1:-1]+con*(y[2:]-2*y[1:-1]+y[:-2])
            yNew[-1] = yNew[-2]
            yOld = copy(y)
            y = copy(yNew)
            self.t += self.dt
            plt.plot(x,y, "fuchsia")
            ax = plt.axes()
            ax.set_facecolor("k")
            plt.ylim(-lim,lim)
            plt.draw()
            plt.pause(0.01)
            plt.clf()
            #this is for the Fourier analysis. This will be stored in the member variable for later. 
            self.f.append(y[int(0.40*steps)-1])

    #coding lesson's function, finds Fourier coefficients
    def DFT(self, samples):
        N = len(samples)
        gamma = []
        for k in range(N//2+1): 
            gammaK = 0
            for n,yn in enumerate(samples):
                gammaK += yn * exp(-2j * pi * k * n/N )
            gamma.append(gammaK/N)

        return gamma

    #finds frequency spectrum of one single point over time on the string
    def spectrumAnalysis(self):
        k = linspace(0,1/(self.runTime*2),len(self.f)//2+1) #the 1/2*time is frequency/2
        self.gamma = self.DFT(self.f)
        plot2 = plt.plot(k, abs(self.gamma))
        plt.show(plot2)
        
wave1 = Wave(1,0.1,3)
wave1.propagateWave()
wave1.spectrumAnalysis()
