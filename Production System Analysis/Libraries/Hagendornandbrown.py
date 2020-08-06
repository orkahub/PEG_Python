# coding=utf-8
import FluidProps
import math


def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """Function to Calculate the Flowing Pressure Gradient by the Method of Beggs and Brill"""
    #P          pressure, psia
    #T          temperature, °F
    #oil_rate   oil flowrate, stb/d
    #wtr_rate   water flowrate, stb/d
    #Gor        producing gas-oil ratio, scf/stb
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    #wtr_grav   water specific gravity
    #d          pipe I.D., in.
    #angle      angle of pipe inclination in degrees
    #               90° = vertical
    #               0°  = horizontal

    print( 'P', P, ' T ', T, ' WOR ', oil_rate/wtr_rate, ' GOR ', Gor)
    
    #Set constants
    pi = math.pi   #4 * math.atan(1)                                               #Define pi
    Psep = 114.7                                                        #Separator pressure, psia
    Tsep = 50                                                           #Separator temperature, °F
   
    #Convert pipe angle from degrees to radians
    angle = angle * pi / 180
    
    #Calculate fluid properties
    Z = FluidProps.zfact((T + 460) / FluidProps.Tc(gas_grav), P / FluidProps.Pc(gas_grav))               #Gas compressibility factor
    if T==145.0:
        print ('error')
       
    Wor = wtr_rate / oil_rate                                           #Water-oil ratio, stb/stb
    TDS = FluidProps.salinity(wtr_grav)                                            #Water salinity, wt% total dissolved solids
    Pb = FluidProps.Pbub(T, Tsep, Psep, gas_grav, oil_grav, Gor)                   #Bubble point pressure, psia
    Rso = FluidProps.sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav)             #Solution gas-oil ratio, scf/stb
    Rsw = FluidProps.sol_gwr(P, T, TDS)                                            #Solution gas_water ratio, scf/stb
    Bo = FluidProps.oil_fvf(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)         #Oil formation volume factor, rb/stb
    Bw = FluidProps.wtr_fvf(P, T, TDS)                                             #Water formation volume factor, rb/stb
    Bg = FluidProps.gas_fvf(P, T, gas_grav)                                        #Gas formation volume factor, ft_/scf
    muo = FluidProps.oil_visc(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)       #Oil viscosity, cp
    muw = FluidProps.wtr_visc(P, T, TDS)                                           #Water viscosity, cp
    mug = FluidProps.gvisc(P, (T + 460), Z, gas_grav)                                #Gas viscosity, cp
    rhoo = FluidProps.oil_dens(T, P, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)  #Oil density, lb/ft_
    rhow = 62.368 * wtr_grav / Bw                                                  #Water density, lb/ft_
    rhog = 2.699 * gas_grav * P / (T + 460) / Z                                    #Gas density, lb/ft_
    sigo = FluidProps.oil_tens(P, T, oil_grav)                                     #Gas-oil interfacial tension, dynes/cm
    sigw = FluidProps.wtr_tens(P, T)                                               #Gas-water interfacial tension, dynes/cm
   
    #Volume fraction weighted liquid properties
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)              #Liquid density
    mul = (Bw * Wor * rhow) / (Bw * Wor * rhow + Bo * rhoo) * muw + (Bo * rhoo) / (Bw * Wor * rhow + Bo * rhoo) * muo             #Liquid viscosity
    sigl = (Bw * Wor * rhow) / (Bw * Wor * rhow + Bo * rhoo) * sigw + (Bo * rhoo) / (Bw * Wor * rhow + Bo * rhoo) * sigo           #Gas-liquid interfacial tension
    
    #Calculate downhole fluid flowrates in ft_/s
    qo = Bo * oil_rate / 15387                                          #Oil flowrate
    qw = Bw * Wor * oil_rate / 15387                                    #Water flowrate
    ql = qo + qw                                                        #Liquid flowrate
    if ((Gor - Rso) < 0):                                        #If gas flowrate is negative, set to zero
        qg = 0
    else:
        qg = Bg * (Gor - Rso - Rsw * Wor) * oil_rate / 86400
    
        
    #Calculate fluid superficial velocities in ft/s
    Axs = pi / 4 * (d / 12) ** 2                                         #X-sectional area of pipe, ft_
    usl = ql / Axs                                                      #Liquid superficial velocity
    usg = qg / Axs                                                      #Gas superficial velocity
    um = usl + usg                                                      #Mixture superficial velocity
    

    #Determine flow regime
    A= 1.071 -((0.2218*um**2)/(d))

    if A<0.13:
        A=0.13
    B= usg/(um)

    if B-A>=0:
        ## use Griffing Lquid holp up Correlation 
        us=0.8*0.3048
        x=(1+um/us)**2-4*usg/us
        HL=1-0.5*(1+um/us-math.sqrt(x))

    else:
        ## Use H&B holdup correlation
        ## Calculate liquid viscosity number and coefficient
        g=32.174 
        NL=0.15726*mul*(1/(rhol*sigl**3))**0.25
        CNL=0.061*NL**3-0.0929*NL**2+0.0505*NL+0.0019
        #print ' CNL ', CNL  
        ##CNL=(0.0019+0.0322*NL-0.6642*NL**2+4.9951*NL**3)/(1-10.0147*NL+33.8696*NL**2+277.2817*NL**3)

        #if NL<== 0.002:
        #    CNL=0.0019
        #if NL>= 0.4:
        #    CNL=0.0115

        #Calculate liquid, gas velocity number and pipe diametre number
        # 0
        NLv=1.938*usl*(rhol/(sigl))**0.25  
        NGv=1.938*usg*(rhol/(sigl))**0.25  
        ND=120.872*d/12*math.sqrt(rhol/(sigl))

        if NGv==0:
            H=0
        else:   
            H=NLv/(NGv**0.575)*(P/14.7)**0.1*CNL/ND
    
        H_Phi=math.sqrt((0.047+(1123.32)*H+729489.64*H**2)/(1+1097.1556*H+722153.97*H**2))

        B=NGv*(NLv**0.38)/(ND**2.14)

        if B<=0.025:
            PHI=27170*B**3-317.52*B**2+0.5472*B+0.9999
        if (B>0.025 and B<=0.055):
            PHI=-5333.33*B**2+58.524*B+0.1171
        if B>0.055:
            PHI=2.5714*B+1.5962

        HL=H_Phi*PHI

     #Liquid velocity number
    laml = usl / um                                                     #Input liquid fraction
    lamg = 1 - laml                                                     #Input gas fraction
    L1 = 316 * laml ** 0.302                                             #Dimensionless constants
    L2 = 0.0009252 * laml ** -2.4684
    L3 = 0.1 * laml ** -1.4516
    L4 = 0.5 * laml ** -6.738

    yl=HL
    yg=1-HL

    #Calculate fluid mixture properties
    rhom = rhol * laml + rhog * lamg                                      #Input fraction weighted density, lb/ft_
    #mum = mul * laml + mug * lamg                                        #Input fraction weighted viscosity, cp
    mum=mul**yl*mug**(yg)
    rhobar = rhol * yl + rhog * yg                                       #In-situ average density, lb/ft_
    
    #Calculate friction factor
    Nre = 1488 * rhom * um * (d / 12) / mum                              #Reynolds number
    fn = Fric(Nre, 0.0006)   
    #print 'friction ' , fn , ' Nre ', Nre     ,' HL ' ,HL , ' Oil ', oil_rate ,' water ', wtr_rate , ' Pressure ', P                                      #No-slip friction factor
    x = laml / HL ** 2
    if ((x > 1) and (x < 1.2)):
        s = math.log(2.2 * x - 1.2)
    else:
        s = math.log(x) / (-0.0523 + 3.182 * math.log(x) - 0.8725 * (math.log(x)) ** 2 + 0.01853 * (math.log(x)) ** 4)
    
    ftp = fn * math.exp(s)                                                    #Two-phase friction factor
    
    #Calculate gradients
    Pgrad_pe = rhobar * math.sin(angle) / 144                                 #Potential energy pressure gradient, psi/ft
    Pgrad_f = 2 * ftp * rhom * um ** 2 / 32.17 / (d / 12) / 144           #Frictional pressure gradient, psi/ft
    Ek = um * usg * rhobar / 32.17 / P / 144                              #Kinetic energy factor
    return (Pgrad_pe + Pgrad_f) / (1 - Ek)                               #Overall pressure gradient, psi/ft



def Fric(Nre, eps):
    """Calculate Fanning Friction Factor using the Chen Equation """
    try:
        math.log
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * math.log10((eps ** 1.1098 / 2.8257) + (7.149 / Nre) ** 0.8981) ) 
    except Exception as inst:
         print(type(inst))    # the exception instance
         print(inst.args)     # arguments stored in .args
         print(inst)       
    
    return (1 / Temp) ** 2




def Pwf_q(FWHP, FWHT,Oil_Rate,Water_Rate,GOR,GasGrav,API, WaterGrav, ID, Angle, Depth, FBHT):
    """Function to calculate the Pwf as function of rate"""
    DPs = []          ## Start as the empty list
    Temps=[]
    Press=[]
    PressList=[]
    DepthList=[]
    i=0
    PressList.append(FWHP)
    Temps.append(FWHT)
    DepthList.append(0)
    DPs.append(0)
    nSteps=60
    i=1
    Tgrad= (FBHT-FWHT)/ Depth

    while (i<= nSteps):
    
        DeltaD= Depth/nSteps*i
        DepthList.append(DeltaD) 
    
        T=FWHT+Tgrad*DeltaD
        Temps.append(T)
        p=PressList[i-1]+DPs[i-1]*(DepthList[i]-DepthList[i-1])     
        dp= Pgrad(p,T,Oil_Rate,Water_Rate,GOR,GasGrav,API, WaterGrav, ID, Angle)
        DPs.append(dp)
        
        PressList.append(p)
        i=i+1

    return PressList[-1]







#print(Pgrad(150,101,100,50,300,0.65,35,1.07,2.44,90))