#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 17:56:50 2021

@author: zow
"""
#%%
import numpy as np
import scipy

import scipy.fftpack as syfp
scipy.__version__
import scipy.signal
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from numpy import linspace, loadtxt
import pylab as pyl


#from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
from matplotlib.font_manager import FontProperties

#%%
#var declarations
#var declarations
Vtlr = 10
Ktlr = 300
Otlr = 0.001
Vnkc = 10000
Knkc = 500
Cni =0.001
Vik = Vnkc
Kik = Knkc
Oik = 10
Inkn = 0.1
Onkn = 0.5
Knkn = 10
scale = 1e22

#second part

Vmcx= 1 
Kmcx =0.05
Omcx =1
Icx = 0.9
Ocx = 0.5
Vpg = 0.3
Kpg = 0.1




parameters = [Vtlr, Ktlr, Otlr, Vnkc, Knkc, Cni, Vik, Kik, Oik, Inkn, Onkn, Knkn]
#param_labels = ["Vtlr", "Ktlr", "Otlr", "Vnkc", "Knkc", "Cni", "Vik", "Kik", "Oik", "Inkn", "Onkn", "Knkn"]

parameters2 = [Vtlr, Ktlr, Otlr, Vnkc, Knkc, Cni, Vik, Kik, Oik, Inkn, 
               Onkn, Knkn, Vmcx, Kmcx, Omcx, Vpg, Kpg]

#param_labels2 = ["Vtlr", "Ktlr", "Otlr", "Vnkc", "Knkc", "Cni", "Vik", "Kik", "Oik", "Inkn", "Onkn", "Knkn", "Vmcx", "Kmcx", "Omcx", "Vpg", "Kpg"]

#Initial conditions
aTLR0 = 1
NKc0 = 6
IK0 = 10
NKn0= 1 
mCX0 = 1
CX0 = 1
PG0 = 1

t0 = -1000
tf = 20000

# time points
t = np.arange(t0,tf,1)

initial_conditions = [aTLR0, NKc0, IK0, NKn0]

var_labels= ["aTLR", "NKc", "IK", "NKn"]
initial_conditions2 = [aTLR0, NKc0, IK0, NKn0, mCX0, CX0, PG0]
var_labels2= ["aTLR", "NKc", "IK", "NKn", "mCX", "CX", "PG"]


#%%
#EQs

#LPS [ t _] : = If [ 1000 ≤ t ≤ 1500, 500, 0 ] ;
def Lps(t):
    
    if (1000 <= t) and (t <= 1500):
        dydt = 500
    else:
        dydt = 0
    
    return dydt

#funcion oscilatoria
#ys = 7 * np.sin(15 * 2 * np.pi * xs) + 3 * np.sin(13 * 2 * np.pi * xs)

#aTLR ' [ t ] == Vtlr * LPS [ t ] / ( Ktlr + LPS [ t ]) - Otlr * aTLR [ t ]
def ddtaTLR(aTLR,t):
    num = Vtlr * Lps(t) 
    den = abs(Ktlr + Lps(t)) 
    subst=  Otlr * aTLR  
    dydt = (num/den - subst) 
    return dydt

#NKc ' [ t ] == Vnkc * aTLR [ t ] / ( Knkc + aTLR [ t ]) - Cni * IK [ t ] * NKc [ t ] - Inkn * NKc [ t ] ,

def ddtNKc(aTLR, IK, NKc):
    num = Vtlr * aTLR  
    den = Knkc + aTLR
    subst=  Cni * IK * NKc - Inkn * NKc
    dydt = num/den - subst
    return dydt


#IK ' [ t ] == Vik * aTLR [ t ] / ( Kik + aTLR [ t ]) - Oik * IK [ t ] - Cni * IK [ t ] * NKc [ t ] ,

def ddtIK(aTLR, IK, NKc):
    num = Vik * aTLR
    den = Kik + aTLR
    subst=  Oik * IK - Cni * IK * NKc
    dydt = num/den - subst
    return dydt


#NKn ' [ t ] == Inkn * NKc [ t ] / ( Knkn + NKc [ t ]) - Onkn * NKn [ t ] ,

def ddtNKn(NKc, NKn):
    num = Inkn * NKc
    den = Knkn + NKc
    subst=   Onkn * NKn 
    dydt = num/den - subst
    return dydt

def ddtmCX(NKn, mCX):
    
    num = Vmcx * NKn
    den = Kmcx + NKn
    if den == 0:
        
        den = 1
    subst = Omcx *mCX
    dydt = num/den - subst
    return dydt

def ddtCX(mCX, CX):
    
    dydt = Icx * mCX - Ocx * CX
    return dydt

def ddtPG(CX):
    num = Vpg * CX
    den = Kpg + CX
    dydt = num/den 
    return dydt


#%%
#Source https://scipy-cookbook.readthedocs.io/items/CoupledSpringMassSystem.html
#def EQsystem(state_vars, t, parameters):
def EQsystem( t,state_vars, parameters):
    
    aTLR, NKc, IK, NKn = state_vars
    Vtlr, Ktlr, Otlr, Vnkc, Knkc, Cni, Vik, Kik, Oik, Inkn, Onkn, Knkn = parameters

    # Create f = (x1',y1',x2',y2'):
    
    f = [ddtaTLR(aTLR,t),
         ddtNKc(aTLR,IK, NKc),
         ddtIK(aTLR, IK, NKc),
         ddtNKn(NKc, NKn)]
    f = np.divide(f, scale)
    return f

#%%
    
def EQsystem2(t, state_vars, parameters):
    
    aTLR, NKc, IK, NKn, mCX, CX, PG = state_vars
    Vtlr, Ktlr, Otlr, Vnkc, Knkc, Cni, Vik, Kik, Oik, Inkn, Onkn, Knkn, Vmcx, Kmcx, Omcx, Vpg, Kpg = parameters
 
    
    f = [ddtaTLR(aTLR,t),
         ddtNKc(aTLR,IK, NKc),
         ddtIK(aTLR, IK, NKc),
         ddtNKn(NKc, NKn),
         ddtmCX(NKn, mCX),
         ddtCX(mCX, CX),
         ddtPG(CX)]
    return f

#%%
    
# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6


#%%

#original problem     
#sol = solve_ivp(fun = EQsystem, y0= initial_conditions,  t_span = (t0, tf), args=(parameters,), method="LSODA",  atol=abserr, rtol=relerr)
#print("success: ", sol.success)




#%%
        
#extended version
sol2 = solve_ivp(fun = EQsystem2, y0= initial_conditions2,  t_span = (t0, tf), args=(parameters2,), method="LSODA",  atol=abserr, rtol=relerr)
print("success: ", sol2.success)



#%%
#PLOTS



def plotConcentrations(sol, labels):
    
    plt.figure(figsize=(5, 4))
    
    plt.xlabel('$t$') # the horizontal axis represents the time 
    plt.ylabel('Concentration')
    plt.title('Complete model')
    #plt.yscale('log')
    
    for i in range(0,3):#range(sol.y.shape[0]):
                
        lab = labels[i]        
        plt.plot(sol.t, sol.y[i], label=lab)
        
    
    #plt.axis([-1000, 6000, 0, 5000])

    plt.legend() # show how the colors correspond to the components of X
    plt.show()
    
    
    
    
 #incomplete model  
#plotConcentrations(sol, var_labels)    
    
#complete model
plotConcentrations(sol2, var_labels2)


#%%
# Do FFT analysis of array
aTLRcurve = sol2.y[0]
NKccurve = sol2.y[1]
IKcurve = sol2.y[2]
wave_array = NKccurve
FFT = scipy.fft.fft(wave_array)
up_limit = min(len(t), len(wave_array))
# Getting the related frequencies
freqs = syfp.fftfreq(len(wave_array), 1) 


# Create subplot windows and show plot
pyl.subplot(211)

pyl.plot(t[0:up_limit], wave_array[0:up_limit], color='y')
pyl.xlabel('Time')
pyl.ylabel('Amplitude')
pyl.subplot(212)

pyl.plot(freqs, np.lib.scimath.log10(abs(FFT)), '.', color='y')  
pyl.xlim(-.05, .05)                       
pyl.show()