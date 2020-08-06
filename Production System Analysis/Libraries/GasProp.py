import math
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
    b = (0.62 - 0.23 * Tr) * Pr + (0.066 / (Tr - 0.86) - 0.037) * Pr ** 2 + 0.32 * Pr ** 6 / (10 ** (9 * (Tr - 1)))
    c = (0.132 - 0.32 * math.log(Tr) / math.log(10))
    d = 10 ** (0.3106 - 0.49 * Tr + 0.1824 * Tr ** 2)
    return a + (1 - a) * math.exp(-b) + c * Pr ** d

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
