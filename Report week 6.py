# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:23:22 2021

@author: Sabrina
"""

from numpy import exp, linspace, sqrt, copy, pi, zeros_like, abs
from numpy.fft import fft
from matplotlib import pyplot as plt

L = 1.0
x0 = 0.4*L
sigma = 0.1
t = 0
v = 3
dt = 0.0099/v  #auto-adjusts to avoid r > 1
runTime = 0.5

def Gauss(x):
    return exp(-(x-x0)**2/(2*sigma**2))/(sigma*sqrt(2*pi))

x, dx = linspace(0, L, 100, retstep = True)
print(dx)
y = Gauss(x)
yNew = zeros_like(y)
yOld = copy(y)

con = v**2*dt**2/dx**2
print(con)
plt.plot(x,y)
plt.show()

#ax = plt.axes()
f = []

while(t<runTime):
    yNew[1:-1] = 2*y[1:-1]-yOld[1:-1]+con*(y[2:]-2*y[1:-1]+y[:-2])
    yNew[-1] = y[-2]
    yOld = copy(y)
    y = copy(yNew)
    t += dt
    plt.plot(x,y, "fuchsia")
    plt.ylim(-4,4)
    plt.draw()
    plt.pause(0.01)
    plt.clf()
    f.append(y[39])
     
fS = 10000  #Rate at which I sample the function
dt = 1/fS  # time between adjacent samples
samplingTime = 1.0  # How long do I sample my function
nSamples = int(samplingTime/dt)  # How many samples am I going to take

def DFT(samples):
    N = len(samples) #length of dataset
    gamma = [] #create new array for all the coefficients
    for k in range(N//2+1): #I'm not sure why we only go halfway. I feel like I remember something about it double 
        #counting something
        gammaK = 0
        for n,yn in enumerate(samples):
            gammaK += yn * exp(-2j * pi * k * n/N )
        gamma.append(gammaK/N)

    return gamma

k = linspace(0,runTime//2,len(f)//2+1)
#t = linspace(0,samplingTime,nSamples) 
#f = Gauss(x)
gamma = DFT(f)
print(gamma)
plt.plot(k, abs(gamma))
plt.show()