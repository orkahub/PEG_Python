import math
def Pbub(T, Tsep, Psep, gas_grav, oil_grav, GOR):
    """ CFunction to Calculate Bubble Point Pressure in psia using Standing Correlation"""
    #T          temperature, °F
    #Tsep       separator temperature, °F
    #Psep       separator pressure, psia
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    #GOR        producing gas-oil ratio, scf/stb
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav<= 30) :
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    else:
        C1 = 0.0178
        C2 = 1.187
        C3 = 23.931
    
    Pbubl = (GOR / (C1 * gas_grav_corr * math.exp(C3 * oil_grav / (T + 460)))) **(1 / C2)
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
    elif (T >= 100):
        st = s100
    else:
        st = s68 - (T - 68) * (s68 - s100) / 32
    
    c = 1 - 0.024 * P ** 0.45
    so = c * st
    if so < 1:
        so = 1
    
    return so


#print(oil_dens(150,2500,65,90,2500,1.23,300,.65,35))


