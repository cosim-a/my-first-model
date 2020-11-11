# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 09:30:32 2020

@author: Cosima
"""
# this is just testing around
t = 1
e = 2
s = 3
t = 4
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.integrate import odeint
time = pd.read_csv("D:\Cosima\programming\R\lena_time.csv",header=None)
co2_measured = pd.read_csv("D:\Cosima\programming\R\lena_delta_co2.csv",header=None)
probe = pd.read_csv("D:\Cosima\programming\R\probes.csv",header=None)
c_initial = pd.read_csv("D:\Cosima\programming\R\c_init.csv",header=None)


# paramter: 
# k1 0.0155    max specific assimilation rate
# k2 0.0000    may specific depolymerization rate
# k3 0.0209    half saturation constant for the assimilation of C_doc by C_bio
# k4 3.7918    half saturation constant for the depolymerization of C_pom
# k5 0.0039    mortality rate
# k6 0.6212    unsoluble fraction of dead microbes  (1-k6 = soluble fraction)
# k7 0.6371    CUE
# k8 0.0353    initial fraction of C_bio
# k9 0.4267    initial fraction of C_pom ( (initial_C - initial_C*k8) * k9)
#[#k10 0.0001    k adsorption
#k11 0.0005    k desorption
#k12 0.1900    may DOC adsorption capacity]  # k10-k12 not used

initial_C = 0.0388
# set parameters
maxSpecAssimilRate              = 0.0155
maxSpecDepolymRate              = 0.0000
halfSaturationConstantAssimil   = 0.0209
halfSaturationConstantDepolym   = 3.7918
mortalityRate                   = 0.0039
unsolubleFractionOfDeadMicrobes = 0.6212
CUE                             = 0.6371
initialFractionCbio             = 0.0353
#initialFractionCpom             = initial_C - (initial_C * initialFractionCbio) * 0.4267
#initialFractionCdoc             = initial_C - (initial_C * initialFractionCbio) * (1-initialFractionCpom)

# initial conditions:
C_bio = initial_C * initialFractionCbio    
C_pom = initial_C - (initial_C * initialFractionCbio) * 0.4267
C_doc = initial_C - (initial_C * initialFractionCbio) * (1-0.4267)
C_atm = 0

# tendencies:
depolymerization = maxSpecDepolymRate * C_pom * C_bio / (C_bio + halfSaturationConstantDepolym)
assimilation     = maxSpecAssimilRate * C_bio * C_doc / (C_doc + halfSaturationConstantAssimil)

deltaC_pom       = unsolubleFractionOfDeadMicrobes * mortalityRate * C_bio - depolymerization
deltaC_doc       = (1-unsolubleFractionOfDeadMicrobes) * mortalityRate * C_bio + depolymerization - assimilation
deltaC_bio       = assimilation * CUE - mortalityRate * C_bio

deltaC_atm = assimilation * (1-CUE)
ks = (maxSpecAssimilRate, maxSpecDepolymRate, halfSaturationConstantAssimil, halfSaturationConstantDepolym, mortalityRate, unsolubleFractionOfDeadMicrobes, CUE, initialFractionCbio)
Cs = (C_bio, C_pom, C_doc, C_atm)

# Euler forward scheme:
# set run parameters
startTime = 0.0
nTimeSteps = 1492
    
myTime = startTime
deltaT = 1
nTimeSteps = 1492
# time loop
#x0 = 0  startTime
#y0 = 1  initialState (C_atm  0)
#xf = 10  max x = stopTime
# n = 100 number of timestamps
# xf = 10
# y0 = 1
# n = 100
# deltax = (xf-y0)/(n-1)

x = np.linspace(myTime, nTimeSteps-1, nTimeSteps)
y = np.zeros([nTimeSteps])

initial_C = 0.0388
#y[0] = y0

pom = list()
doc = list()
bio = list()
bio = list()
co2 = list()
hlp_co2 = 0
hlp_pom = C_pom
hlp_doc = C_doc
hlp_bio = C_bio

#y = odeint(compute_tendencies, y0, t)

for i in range(0, nTimeSteps):   
   
    for j in range(0, deltaT):
        depolymerization = maxSpecDepolymRate * hlp_pom * hlp_bio / (hlp_bio + halfSaturationConstantDepolym)
        assimilation     = maxSpecAssimilRate * hlp_bio * hlp_doc / (hlp_doc + halfSaturationConstantAssimil)


        deltaC_pom       = unsolubleFractionOfDeadMicrobes * mortalityRate * hlp_bio - depolymerization
        deltaC_doc       = (1-unsolubleFractionOfDeadMicrobes) * mortalityRate * hlp_bio + depolymerization - assimilation
        deltaC_bio       = assimilation * CUE - mortalityRate * hlp_bio
        deltaC_atm       = assimilation * (1-CUE)

        hlp_pom = hlp_pom + deltaC_pom 
        hlp_doc = hlp_doc + deltaC_doc
        hlp_bio = hlp_bio + deltaC_bio
        hlp_co2 = hlp_co2 + deltaC_atm
        
    pom.append(hlp_pom) 
    doc.append(hlp_doc)
    bio.append(hlp_bio)
    co2.append(hlp_co2)
    myTime = myTime + deltaT
        
plt.plot(x, co2, 'k-.',   label = 'co2')
plt.plot(x, pom, 'r-',  label = 'pom')
plt.plot(x, doc, 'b--', label = 'doc')
plt.plot(x, bio, 'g:',  label = 'bio')
plt.xlabel("Value of x")
plt.ylabel("Value of y")
plt.legend()
plt.title("Approximate Solution with Forward Euler's Method")
plt.show()

p1 = pd.DataFrame(data = time[1:])
#p1 = p1[1:0]

