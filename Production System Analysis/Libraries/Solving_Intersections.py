# coding=utf-8
from __future__ import division 
import math
import FluidProps
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import Vogel
import BeggsandBrill as BB
from scipy.optimize import fsolve
import Hagendornandbrown as HB

Oil_Rate=100
Water_Rate=50.0
GasGrav=0.65
WaterGrav=1.07                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
ID=2.44
Angle=90.0
FBHT=150.0
FWHP=150.0
FWHT=100.0
Depth=5000.0
nSteps=20
API=30.0
Tgrad= (FBHT-FWHT)/ Depth

Pressure=3000.0
Thickness=75.0
GasGrav=0.65
GOR= 375.0
Temp=150.0
rw=0.328
re=1053.0
k=75.0
s=-1.5
visc= 1.55
Psat = FluidProps.Pbub(Temp,75.0,100.0,GasGrav, API, GOR)
Wcut= Water_Rate/(Oil_Rate+Water_Rate)

def Press_Function(rateArray, pressArray, x):
        return np.interp(x,rateArray, pressArray)

IPR=Vogel.Vogel_DarcyIPR(Pressure,k,Thickness, visc,re,rw,s, 1.21, Temp, Psat, 20)

BB_Pwf=[]
BB_Rate=[]
for i in IPR[0]:
    if i==0:
        rate=0.1
    else:
        rate=i  
    wr=rate*Wcut/(1-Wcut)   
    #print 'req rate '   ,rate
    Pf=HB.Pwf_q(FWHP,FWHT,rate, wr,GOR,GasGrav,API,WaterGrav,ID,Angle, Depth, FBHT)
    BB_Pwf.append(Pf)
    BB_Rate.append(rate)
    #print Pf
#print "PWF " + str(Pf)

x= np.asarray(BB_Rate)
f= np.asarray(IPR[1])
g= np.asarray(BB_Pwf)

q=x

QNew =np.linspace(0.1, BB_Rate[-1],10000)
P_IPR=np.interp(QNew, x, f)
P_VLP=np.interp(QNew, x, g)

SolRate= [fsolve(lambda x:Press_Function(QNew,P_IPR,x)-Press_Function(QNew,P_VLP,x), 0.0)]
print0(SolRate)

plt.plot(QNew, P_IPR )
plt.plot(QNew, P_VLP )

#def findIntersection(fun1, fun2, x0):
    #return fsolve(lambda x:fun1[x]-fun2[x], x0)

#t= findIntersection(P_IPR,P_VLP,0.0)

#print t

idx = np.argwhere(np.diff(np.sign(P_IPR - P_VLP)) != 0.0).reshape(-1) + 0
plt.plot(IPR[0], IPR[1])
plt.plot(BB_Rate, BB_Pwf )

plt.plot(QNew[idx], P_IPR[idx], 'ro', ms=10)

plt.ylabel('Pwf')
plt.xlabel('Rate')
#plt.gca().invert_yaxis()
plt.show()

