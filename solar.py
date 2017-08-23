import math
def solar_declination_angle(caldat):
    dec_angle = 23.45*math.cos(2*math.pi*(caldat-172.)/365.)
    return dec_angle

#Computes the solar insolation at the TOA as a function of
#latitude and calendar day expressed value between 1-365.
#Method based on Coakley, 1979, JGR, 36, 260-269. And for daily see:
#http://nit.colorado.edu/atoc5560/week13.pdf
def solar_radiation_daily(latitude,caldat):
    PHI = latitude*math.pi/180.
    DEC = solar_declination_angle(caldat)*math.pi/180.
    Hdeg = (-1.)*math.tan(PHI)*math.tan(DEC)

    if Hdeg > -1.0 and Hdeg < 1.0:
        H   = math.acos( (-1.)*math.tan(PHI)*math.tan(DEC) )
    else:
        H = math.pi

    day_length = 24.*H/math.pi #hours

    t = (2 * math.pi * caldat) / 365

    #Sun-earth distance correction
    sfac = 1.000110 + 0.034221*math.cos(t) + 0.001280*math.sin(t) + 0.000719*math.cos(2*t) + 0.000077*math.sin(2*t)

    #Solar Constant
    S0 = 1361.

    #solar insolation averaged over the full day at the specified latitude
    F = (S0/math.pi)*sfac*(H*math.sin(PHI)*math.sin(DEC)+math.sin(H)*math.cos(PHI)*math.cos(DEC))

    if F < 0.:
        F=0.

    return F
