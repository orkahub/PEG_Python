import math
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
