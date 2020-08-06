import matplotlib.pyplot as plt
import BeggsandBrill as BB
import Hagendornandbrown   as HB
import FluidProps

DPs = []          ## Start as the empty list
Temps=[]
Press=[]
PressList=[]
DepthList=[]
#list.append('a')   ## Use append() to add elements
#list.append('b')

Oil_Rate=682.0
Water_Rate=76.0
GOR=84.0
GasGrav=0.7
WaterGrav=1.05                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
ID=2.995
Angle=90
FBHT=180.0
FWHP=200.0
FWHT=80.0
Depth=9700.0
nSteps=20
API=40.0
Tgrad= (FBHT-FWHT)/ Depth

i=0
PressList.append(FWHP)
Temps.append(FWHT)
DepthList.append(0)
DPs.append(0)

i=1
while (i<= nSteps):
   
    DeltaD= Depth/nSteps*i
    DepthList.append(DeltaD) 
   
    T=Temps[i-1]+Tgrad*DeltaD
    Temps.append(T)
    p=PressList[i-1]
    dp= HB.Pgrad(p,T,Oil_Rate,Water_Rate,GOR,GasGrav,API, WaterGrav, ID, Angle)
    DPs.append(dp)
    #print dp
    p=PressList[i-1]+DPs[i]*(DepthList[i]-DepthList[i-1])     
   
    
    PressList.append(p)
    #print p  ,  DeltaD
    
    i=i+1



plt.plot(PressList, DepthList)
plt.ylabel('Depth')
plt.xlabel('Pressure')
plt.gca().invert_yaxis()
plt.show()

