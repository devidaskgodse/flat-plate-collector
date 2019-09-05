from math import sin, cos, pi


# Input
d = pi / 180



phi = d * float(input("Enter latitude in degrees: "))
longitude = float(input("Enter longitude in degrees: "))
lsm = float(input("Enter longitude of standard time meridian in degrees: "))

Ig = float(input("Enter total irradiration on horizontal: "))
Id = float(input("Enter diffuse irradiration on horizontal: "))
n = float(input("Enter no. of the day: "))
LST = float(input("Enter time of the day in hours: "))

gamma = d * float(input("Enter surface azimuth angle in degrees: "))
beta = d * (float(input("Enter tilt angle of surface in degrees: ")))


def IT(phi, longitude, lsm, Ig, Id, n, LST, gamma, beta):
    # Process
    Ib = Ig - Id

    delta = d * (23.45 * sin(360 * (284.0 + n) * d / 365.0))

    B = d * (360.0 * (n - 1) / 365.0)
    ET = (9.87 * sin(2.0*B) - 7.53 * cos(B) - 1.5 * sin(B)) / 60.0
    ST = LST + ET + (lsm - longitude) / 15.0
    omega = d * (15 * (ST - 12.0))

    costheta = sin(phi) * sin(delta) * cos(beta) + \
               + sin(phi) * cos(delta) * cos(gamma) * cos(omega) * sin(beta) \
               + cos(phi) * cos(delta) * cos(omega) * cos(beta)\
           - cos(phi) * sin(delta) * cos(gamma) * sin(beta)\
           + cos(delta) * sin(gamma) * sin(omega) * sin(beta)
    costhetaz = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega)

    ro = 0.2
    Rb = costheta / costhetaz
    Rd = cos(beta * d / 2.0)**2.0
    Rr = sin(beta * d / 2.0)**2.0

    # Output
    It = Ib * Rb + Id * Rd + Ig * Rr * ro

    return It

print(IT(phi,longitude,lsm,Ig,Id,n,LST,gamma,beta))

'''
Inputs-
Enter latitude in degrees: 27.96
Enter longitude in degrees: -82.54
Enter longitude of standard time meridian in degrees: -75
Enter total irradiration on horizontal: 714.5160
Enter diffuse irradiration on horizontal: 99
Enter no. of the day: 32
Enter time of the day in hours: 12
Enter surface azimuth angle in degrees: 10
Enter tilt angle of surface in degrees: 30

Output-
952.0124661790134
Value in book = 932
'''
