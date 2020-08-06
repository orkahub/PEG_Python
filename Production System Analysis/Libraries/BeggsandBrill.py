# coding=utf-8
#import sys
#sys.path.append('C:\Users\frank\PEG Commmunity\VLP\Libraries')
import Libraries.FluidProps as FluidProps
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
    
    #Set constants
    pi = math.pi   #4 * math.atan(1)                                               #Define pi
    Psep = 114.7                                                        #Separator pressure, psia
    Tsep = 50                                                           #Separator temperature, °F
   
    #Convert pipe angle from degrees to radians
    angle = angle * pi / 180
    
    #Calculate fluid properties
    Z = FluidProps.zfact((T + 460) / FluidProps.Tc(gas_grav), P / FluidProps.Pc(gas_grav))               #Gas compressibility factor
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
    mug = FluidProps.gvisc(P, T + 460, Z, gas_grav)                                #Gas viscosity, cp
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
    Nfr = um ** 2 / (d / 12) / 32.174                                    #Froude number
    Nvl = 1.938 * usl * (rhol / sigl) ** 0.25                            #Liquid velocity number
    laml = usl / um                                                     #Input liquid fraction
    lamg = 1 - laml                                                     #Input gas fraction
    L1 = 316 * laml ** 0.302                                             #Dimensionless constants
    L2 = 0.0009252 * laml ** -2.4684
    L3 = 0.1 * laml ** -1.4516
    L4 = 0.5 * laml ** -6.738

    regime = Flow_regime(Nfr, laml, L1, L2, L3, L4)
                 
    #Calculate holdups
    if (regime == 2):
        a = (L3 - Nfr) / (L3 - L2)
        yl_seg = Liq_holdup(Nfr, Nvl, laml, angle, 1)
        yl_int = Liq_holdup(Nfr, Nvl, laml, angle, 3)
        yl = a * yl_seg + (1 - a) * yl_int
    else:
        yl = Liq_holdup(Nfr, Nvl, laml, angle, regime)
    
    yg = 1 - yl
    
    #Calculate fluid mixture properties
    rhom = rhol * laml + rhog * lamg                                     #Input fraction weighted density, lb/ft_
    mum = mul * laml + mug * lamg                                        #Input fraction weighted viscosity, cp
    rhobar = rhol * yl + rhog * yg                                       #In-situ average density, lb/ft_
    
    #Calculate friction factor
    Nre = 1488 * rhom * um * (d / 12) / mum                              #Reynolds number
    fn = Fric(Nre, 0.0006)                                                    #No-slip friction factor
    x = laml / yl ** 2
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

def Flow_regime(Nfr, laml, L1, L2, L3, L4):
    """Function to Determine the Flow Regime by the Method of Beggs and Brill"""
    #The function returns a number indicating the flow regime
    #   1 = Segregated flow
    #   2 = Transition flow
    #   3 = Intermittent flow
    #   4 = Distributed flow    
    #Nfr        Froude Number
    #laml       Input liquid fraction
    #L1,2,3,4   Dimensionless constants

    #Regime 1 - Segregated flow
    if (((laml < 0.01) and (Nfr < L1)) or ((laml >= 0.01) and (Nfr < L2))):
        regime = 1
    
        
    #Regime 2 - Transition flow
    if ((laml >= 0.01) and (L2 < Nfr) and (Nfr <= L3)):
        regime = 2
    
        
    #Regime 3 - Intermittent flow
    if ((((0.01 <= laml) and (laml < 0.4)) and ((L3 < Nfr) and (Nfr < L1))) or ((laml >= 0.4) and (L3 < Nfr) and (Nfr <= L4))):
        regime = 3
    
        
    #Regime 4 - Distributed flow
    if (((laml < 0.4) and (Nfr >= L1)) or ((laml >= 0.4) and (Nfr > L4))):
        regime = 4
    
    
    return regime

def Liq_holdup(Nfr, Nvl, laml, angle, regime):
    """Function to Calculate Liquid Holdup for the Segregated, Intermittent and Distributed Regimes
    by the Method of Beggs and Brill"""
    #Nfr        Froude number
    #Nvl        Liquid velocity number
    #angle      pipe inclination in degrees
    #regime     flow regime 1 = segregated, 3 = intermittent, 4 = distributed
    
    #Define constants
    if (regime == 1):
        a = 0.98
        b = 0.4846
        c = 0.0868
        if (angle >= 0):
            d = 0.011
            e = -3.768
            f = 3.539
            g = -1.614
        else:
            d = 4.7
            e = -0.3692
            f = 0.1244
            g = -0.5056
           
    
    if (regime == 3):
        a = 0.845
        b = 0.5351
        c = 0.0173
        if (angle >= 0):
            d = 2.96
            e = 0.305
            f = -0.4473
            g = 0.0978
        else:
            d = 4.7
            e = -0.3692
            f = 0.1244
            g = -0.5056
            
    
    if (regime == 4):
        a = 1.065
        b = 0.5824
        c = 0.0609
        if (angle >= 0):
            d = 1
            e = 0
            f = 0
            g = 0
        else:
            d = 4.7
            e = -0.3692
            f = 0.1244
            g = -0.5056
          
    #Calculate holdup
    corr = (1 - laml) * math.log(d * laml ** e * Nvl ** f * Nfr ** g)
    if (corr < 0):
        corr = 0
    
    psi = 1 + corr * (math.sin(1.8 * angle) - (math.sin(1.8 * angle)) ** 3 / 3)
    ylo = a * laml ** b / Nfr ** c
    if (ylo < laml):
        ylo = laml
    
    yl = ylo * psi
    return yl

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



def  Pgrad2(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """Function to Calculate the Flowing Pressure Gradient by the no slip method"""
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
    
    #Set constants
    pi = math.pi #4 * Atn(1)                                                     #Define pi
    Psep = 114.7                                                        #Separator pressure, psia
    Tsep = 50                                                           #Separator temperature, °F
   
    #Convert pipe angle from degrees to radians
    angle = angle * pi / 180
    
    #Calculate fluid properties
    Z = FluidProps.zfact((T + 460) / Tc(gas_grav), P / Pc(gas_grav))               #Gas compressibility factor
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
    mug = FluidProps.gvisc(P, T + 460, Z, gas_grav)                                #Gas viscosity, cp
    rhoo = FluidProps.oil_dens(T, P, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)  #Oil density, lb/ft_
    rhow = 62.368 * wtr_grav / Bw                                       #Water density, lb/ft_
    rhog = 2.699 * gas_grav * P / (T + 460) / Z                         #Gas density, lb/ft_

    #Volume fraction weighted liquid properties
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)              #Liquid density
    mul = (Bw * Wor * rhow) / (Bw * Wor * rhow + Bo * rhoo) * muw + (Bo * rhoo) / (Bw * Wor * rhow + Bo * rhoo) * muo             #Liquid viscosity
    
    #Calculate downhole fluid flowrates in ft_/s
    qo = Bo * oil_rate / 15387                                          #Oil flowrate
    qw = Bw * Wor * oil_rate / 15387                                    #Water flowrate
    ql = qo + qw                                                        #Liquid flowrate
    if ((Gor - Rso) < 0):                                              #If gas flowrate is negative, set to zero
        qg = 0
    else:
        qg = Bg * (Gor - Rso - Rsw * Wor) * oil_rate / 86400
    
        
    #Calculate fluid superficial velocities in ft/s
    Axs = pi / 4 * (d / 12) ** 2                                         #X-sectional area of pipe, ft_
    usl = ql / Axs                                                      #Liquid superficial velocity
    usg = qg / Axs                                                      #Gas superficial velocity
    um = usl + usg                                                      #Mixture superficial velocity

    #Calculate holdups
    laml = usl / um                                                     #Input liquid fraction
    lamg = 1 - laml                                                     #Input gas fraction
                  
    #Calculate fluid mixture properties
    rhom = rhol * laml + rhog * lamg                                    #Input fraction weighted density, lb/ft_
    mum = mul * laml + mug * lamg                                       #Input fraction weighted viscosity, cp
    
    #Calculate friction factor
    Nre = 1488 * rhom * um * (d / 12) / mum                             #Reynolds number
    fn = Fric(Nre, 0.01)                                                #No-slip friction factor
    
    #Calculate gradients
    Pgrad_pe = rhom * math.sin(angle)                                        #Potential energy pressure gradient, psi/ft
    Pgrad_f = 2 * fn * rhom * um ** 2 / 32.17 / (d / 12)                 #Frictional pressure gradient, psi/ft
    return (Pgrad_pe + Pgrad_f) / 144                                 #Overall pressure gradient, psi/ft



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
    
        T=Temps[i-1]+Tgrad*DeltaD
        Temps.append(T)
        p=PressList[i-1]+DPs[i-1]*(DepthList[i]-DepthList[i-1])     
        dp= Pgrad(p,T,Oil_Rate,Water_Rate,GOR,GasGrav,API, WaterGrav, ID, Angle)
        DPs.append(dp)
        
        PressList.append(p)
        i=i+1

    return PressList[-1]







#print(Pgrad(150,101,100,50,300,0.65,35,1.07,2.44,90))