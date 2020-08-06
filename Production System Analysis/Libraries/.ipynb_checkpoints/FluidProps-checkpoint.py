# coding=utf-8
import math


def Pbub(T, Tsep, Psep, gas_grav, oil_grav, Gor):
    """ CFunction to Calculate Bubble Point Pressure in psia using Standing Correlation"""
    #T          temperature, °F
    #Tsep       separator temperature, °F
    #Psep       separator pressure, psia
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    #Gor        producing gas-oil ratio, scf/stb
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav<= 30) :
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    else:
        C1 = 0.0178
        C2 = 1.187
        C3 = 23.931
    
    Pbubl = (Gor / (C1 * gas_grav_corr * math.exp(C3 * oil_grav / (T + 460)))) **(1 / C2)
    return Pbubl

def correct(Tsep, Psep, gas_grav, oil_grav):
    """Function to Calculate Corrected Gas Gravity"""
    #Tsep       separator temperature, °F
    #Psep       separator pressure, psia
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity

    return  gas_grav * (1 + 5.912 * 10 ** -5 * oil_grav * Tsep * math.log10(Psep / 114.7) / math.log(10))

def  sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav):
    """Function to Calculate Solution Gas-Oil Ratio in scf/stb"""
    #T          temperature, °F
    #P          pressure, psia
    #Tsep       separator temperature, °F
    #Psep       separator pressure, psia
    #Pb         bubble point pressure, psia
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav <= 30):
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    else:
        C1 = 0.0178
        C2 = 1.187
        C3 = 23.931
    
    if (P <= Pb):
        Rs = C1 * gas_grav_corr * P** C2 * math.exp(C3 * oil_grav / (T + 460))
    else:
        Rs = C1 * gas_grav_corr * Pb ** C2 * math.exp(C3 * oil_grav / (T + 460))
    
    return Rs

def oil_fvf(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
    """Function to Calculate Oil Formation Volume Factor in bbl/stb"""
    #'T          temperature, °F
    #P          pressure, psia
    #Tsep       separator temperature, °F
    #Psep       separator pressure, psia
    #Pb         bubble point pressure, psia
    #Rs         solution gas-oil ratio, scf/stb
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav <= 30) :
        C1 = 0.0004677
        C2 = 1.751E-05
        C3 = -1.811E-08
    else:
        C1 = 0.000467
        C2 = 1.1E-05
        C3 = 1.337E-09
    

    if (P <= Pb):
        Bo = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr) + C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
    else:
        Bob = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr)+ C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
        co = oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav)
        Bo = Bob * math.exp(co * (Pb - P))
    
    return  Bo
def  oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav):
    """Function to Calculate Oil Isothermal Compressibility in 1/psi"""
    #'T          temperature, °F
    #'P          pressure, psia
    #'Tsep       separator temperature, °F
    #'Psep       separator pressure, psia
    #'Rs         solution gas-oil ratio, scf/stb
    #'gas_grav   gas specific gravity
    #'oil_grav   API oil gravity
    
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    oil_compr = (5 * Rs + 17.2 * T - 1180 * gas_grav_corr + 12.61 * oil_grav - 1433) / (P * 10 ** 5)
    return oil_compr

def oil_visc(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
    """Function to Calculate Oil Viscosity in cp"""
    #'T          temperature, °F
    #'P          pressure, psia
    #'Tsep       separator temperature, °F
    #'Psep       separator pressure, psia
    #'Pb         bubble point pressure, psia
    #'Rs         solution gas-oil ratio, scf/stb
    #'gas_grav   gas specific gravity
    #'oil_grav   API oil gravity
    
    a = 10.715 * (Rs + 100) ** (-0.515)
    b = 5.44 * (Rs + 150) ** (-0.338)
    Z = 3.0324 - 0.0203 * oil_grav
    Y = 10 **Z
    x = Y * T ** (-1.163)
    visc_oD = 10 ** x - 1
    if (P <= Pb):
        visc_o = a * visc_oD ** b
    else:
        M = 2.6 * P ** 1.187 * math.exp(-11.513 - 8.98E-05 * P)
        visc_ob = a * visc_oD ** b
        visc_o = visc_ob * (P / Pb) ** M

    return visc_o

def oil_dens(T, P, Tsep, Psep, Pb, Bo, Rs, gas_grav, oil_grav):
    """Function to Calculate Oil Density in lb/ft"""
    #'T          temperature, °F
    #'P          pressure, psia
    #'Tsep       separator temperature, °F
    #'Psep       separator pressure, psia
    #'Pb         bubble point pressure, psia
    #'Bo         oil formation volume factor, bbl/stb
    #'Rs         solution gas-oil ratio, scf/stb
    #'gas_grav   gas specific gravity
    #'oil_grav   API oil gravity
    oil_grav_sp = 141.5 / (oil_grav + 131.5)
    if (P <= Pb):
        rho_o = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bo)
    else:
        co = oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav)
        Bob = Bo / (math.exp(co * (P - Pb)))
        rho_ob = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bob)
        rho_o = rho_ob * Bo / Bob
    
    return rho_o

def oil_tens(P, T, oil_grav):
    """Function to Calculate Gas-Oil Interfacial Tension in dynes/cm"""
    #P          pressure, psia
    #T          temperature, °F
    #oil_grav   API oil gravity
    s68 = 39 - 0.2571 * oil_grav
    s100 = 37.5 - 0.2571 * oil_grav
    if (T <= 68):
        st = s68
    elif(T >= 100):
        st = s100
    else:
        st = s68 - (T - 68) * (s68 - s100) / 32
    
    c = 1 - 0.024 * P ** 0.45
    so = c * st
    if (so < 1):
        so = 1
    
    return so

def Tc(grav):
    """Function to Calculate Gas Critical Temperature in °R"""
    #grav       gas specific gravity
    return 169.2 + 349.5 * grav - 74 * grav ** 2

def Pc(grav):
    """Function to Calculate Gas Critical Pressure in psia"""
    #grav       gas specific gravity
    return 756.8 - 131 * grav - 3.6 * grav**2

def zfact(Tr, Pr):
    """Function to Calculate Gas Compressibility Factor"""
    #'Tr         reduced temperatue
    #'Pr         reduced pressure
    a = 1.39 * (Tr - 0.92) ** 0.5 - 0.36 * Tr - 0.101
    b = (0.62 - 0.23 * Tr) * Pr + (0.066 / (Tr - 0.86) - 0.037) * Pr ** 2 + 0.32 * Pr ** 6 / (10** (9 * (Tr - 1)))
    c = (0.132 - 0.32 * math.log10(Tr))
    d = 10 ** (0.3106 - 0.49 * Tr + 0.1824 * Tr **2)
    zfact = a + (1 - a) * math.exp(-b) + c * Pr **d
    return zfact

def  gvisc(P, T, Z, grav):
    """Function to Calculate Gas Viscosity in cp"""
    #P          pressure, psia
    #T          temperature, °R
    #Z          gas compressibility factor
    #grav       gas specific gravity
    M = 28.964 * grav
    x = 3.448 + 986.4 / T + 0.01009 * M
    Y = 2.447 - 0.2224 * x
    rho = (1.4926 / 1000) * P * M / Z / T
    if Y<0 or rho<0:
       print ('epa')
    K = (9.379 + 0.01607 * M) * T ** 1.5 / (209.2 + 19.26 * M + T)
 
    return K * math.exp(x * rho ** Y) / 10000

def  gas_fvf(P, T, grav):
    """Function to Calculate Gas Formation Volume Factor in ft_/scf"""
    #P          pressure, psia
    #T          temperature,°F
    #grav       gas specific gravity
    Tr = (T + 460) / Tc(grav)
    Pr = P / Pc(grav)
    Z = zfact(Tr, Pr)
    return 0.0283 * Z * (T + 460) / P

def wtr_fvf(P, T, TDS):
    """Function to Calculate Water Formation Volume Factor in bbl/stb"""
    #P          pressure, psia
    #T          temperature, °F
    #TDS        total dissolved solids, wt%
    Y = 10000 * TDS
    x = 5.1 * 10 ** -8 * P + (T - 60) * (5.47 * 10 ** -6 - 1.95 * 10 ** -10 * P) + (T - 60) ** 2 * (-3.23 * 10 ** -8 + 8.5 * 10 ** -13 * P)
    C1 = 0.9911 + 6.35E-05 * T + 8.5 * 10 ** -7 * T ** 2
    C2 = 1.093 * 10 ** -6 - 3.497 * 10 ** -9 * T + 4.57 * 10 ** -12 * T ** 2
    C3 = -5 * 10 ** -11 + 6.429 * 10 ** -13 * T - 1.43 * 10 ** -15 * T ** 2
    Bwp = C1 + C2 * P + C3 * P ** 2
    Bw = Bwp * (1 + 0.0001 * x * Y)
    return Bw

def  sol_gwr(P, T, TDS):
    """Function to Calculate Solution Gas-Water Ratio in scf/stb"""
    #P          pressure, psia
    #T          temperature, °F
    #TDS        total dissolved solids, wt%
    Y = 10000 * TDS
    x = 3.471 * T ** -0.837
    C1 = 2.12 + 0.00345 * T - 3.59E-05 * T ** 2
    C2 = 0.0107 - 5.26E-05 * T + 1.48 * 10 ** -11 * T ** 2
    C3 = -8.75 * 10 ** -7 + 3.9 * 10 ** -9 * T - 1.02 * 10 ** -11 * T ** 2
    Rswp = C1 + C2 * P + C3 * P ** 2
    Rsw = Rswp * (1 - 0.0001 * x * Y)
    return Rsw

def wtr_dens(P, T, Bw, TDS):
    """Function to Calculate Water Density in lb/ft"""
    #P          pressure, psia
    #T          temperature, °F
    #Bw         water formation volume factor, bbl/stb
    #TDS        total dissolved solids, wt%
    rho = (62.368 + 0.438603 * TDS + 1.60074 * 10 ** -3 * TDS ** 2) / Bw
    return rho

def wtr_visc(P, T, TDS):
    """Function to Calculate Water viscosity in cp"""
    #P          pressure, psia
    #T          temperature, °F
    #TDS        total dissolved solids, wt%
    Y = 10000 * TDS
    a = -0.04518 + 9.313 * 10 ** -7 * Y - 3.93 * 10 ** -12 * Y ** 2
    b = 70.634 + 9.576 * 10 ** -10 * Y ** 2
    muwd = a + b / T
    mu = muwd * (1 + 3.5 * 10 ** -12 * P ** 2 * (T - 40))
    return mu

def salinity(wtr_grav):
    """Function to Calculate Water Salinity at 60°F and 1 atm"""
    #wtr_grav   specific gravity of water
    rho = 62.368 * wtr_grav
    a = 0.00160074
    b = 0.438603
    c = 62.368 - rho
    s = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
    return s

def wtr_tens(P, T):
    """Function to Calculate Gas-Water Interfacial Tension in dynes/cm"""
    #P          pressure, psia
    #T          temperature, °F
    s74 = 75 - 1.108 * P ** 0.349
    s280 = 53 - 0.1048 * P ** 0.637
    if (T <= 74):
        sw = s74
    elif(T >= 280):
        sw = s280
    else:
        sw = s74 - (T - 74) * (s74 - s280) / 206
    
    if (sw < 1):
        sw = 1
    
    return sw

#print(oil_dens(150,2500,65,90,2500,1.23,300,.65,35))


