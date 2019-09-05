from math import sin, cos, pi
import numpy as np


Jan = [0,0,0,0,0,0,0,13,136,340,529,662,726,711,619,458,254,67,2,0,0,0,0,0]
Feb = [0,0,0,0,0,0,0,19,171,398,600,731,822,820,723,560,341,120,9,0,0,0,0,0]
Mar = [0,0,0,0,0,0,1,54,247,497,714,861,929,909,804,624,389,151,15,0,0,0,0,0]
Apr = [0,0,0,0,0,0,9,120,335,577,779,913,980,960,855,676,442,192,22,0,0,0,0,0]
May = [0,0,0,0,0,0,19,161,358,565,737,852,905,884,791,631,423,200,34,0,0,0,0,0]
Jun = [0,0,0,0,0,0,16,127,273,417,528,599,646,640,580,474,330,170,40,1,0,0,0,0]
Jul = [0,0,0,0,0,0,11,100,210,325,417,478,504,492,443,360,253,137,40,1,0,0,0,0]
Aug = [0,0,0,0,0,0,6,78,194,320,426,500,530,517,463,380,265,134,25,0,0,0,0,0]
Sep = [0,0,0,0,0,0,6,91,230,374,482,547,587,577,516,405,261,109,11,0,0,0,0,0]
Oct = [0,0,0,0,0,0,4,95,290,484,639,729,762,723,617,451,249,63,1,0,0,0,0,0]
Nov = [0,0,0,0,0,0,1,55,235,441,612,722,751,701,579,408,202,30,0,0,0,0,0,0]
Dec = [0,0,0,0,0,0,0,21,164,363,536,649,694,661,555,394,196,33,0,0,0,0,0,0]

Ig = [Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec]
print(len(Ig))
Jan = [0,0,0,0,0,0,0,13,103,188,218,225,227,223,215,202,160,62,2,0,0,0,0,0]
Feb = [0,0,0,0,0,0,0,19,120,204,227,249,239,219,208,200,183,95,9,0,0,0,0,0]
Mar = [0,0,0,0,0,0,1,52,162,230,246,243,209,220,231,223,209,119,15,0,0,0,0,0]
Apr = [0,0,0,0,0,0,9,98,197,244,266,258,230,206,194,180,172,129,22,0,0,0,0,0]
May = [0,0,0,0,0,0,19,120,218,275,314,330,320,296,281,247,204,139,34,0,0,0,0,0]
Jun = [0,0,0,0,0,0,16,108,213,315,365,391,414,391,351,300,227,135,39,1,0,0,0,0]
Jul = [0,0,0,0,0,0,11,91,189,295,361,398,426,407,365,305,231,129,40,1,0,0,0,0]
Aug = [0,0,0,0,0,0,6,75,175,287,371,414,446,422,394,331,236,123,25,0,0,0,0,0]
Sep = [0,0,0,0,0,0,6,83,190,302,379,423,465,433,391,318,208,93,11,0,0,0,0,0]
Oct = [0,0,0,0,0,0,4,77,160,210,243,271,286,274,249,209,151,56,1,0,0,0,0,0]
Nov = [0,0,0,0,0,0,1,51,147,192,203,202,200,203,188,166,125,30,0,0,0,0,0,0]
Dec = [0,0,0,0,0,0,0,21,112,187,209,216,216,212,201,178,130,33,0,0,0,0,0,0]

Id = [Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec]
print(len(Id))

d = pi / 180
print("Latitude for Mumbai: 19.076 N")
phi = d * 19.076
print("Longitude for Mumbai: 72.8777 E")
longitude = 72.8777
print("Mumbai follow Time zone of 82.5 E")
lsm = 82.5

gamma = d * float(input("Enter surface azimuth angle in degrees: "))
beta = d * (float(input("Enter tilt angle of surface in degrees: ")))

def IT(phi, longitude, lsm, Ig, Id, n, LST, gamma, beta):
    # Process
    Ib = Ig - Id

    delta = d * (23.45 * sin(360 * (284 + n) * d / 365))

    B = d * (360 * (n - 1) / 365)
    ET = (9.87 * sin(2*B) - 7.53 * cos(B) - 1.5 * sin(B)) / 60
    ST = LST + ET + (lsm - longitude) / 15
    omega = d * (15 * (ST - 12))

    costheta = sin(phi) * sin(delta) * cos(beta) + \
           + sin(phi) * cos(delta) * cos(gamma) * cos(omega) * sin(beta) \
           + cos(phi) * cos(delta) * cos(omega) * cos(beta)\
           - cos(phi) * sin(delta) * cos(gamma) * sin(beta)\
           + cos(delta) * sin(gamma) * sin(omega) * sin(beta)
    costhetaz = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega)

    ro = 0.2
    Rb = costheta / costhetaz
    Rd = cos(beta * d / 2)**2
    Rr = sin(beta * d / 2)**2

    # Output
    It = Ib * Rb + Id * Rd + Ig * Rr * ro

    return It


ITilted = []
for i in range(len(Ig)):
    ITmonth = []
    for hour in range(len(Ig[i])):
        n = 15 + i * 30
        LST = hour + 0.5
        ITmonth.append(IT(phi,longitude,lsm,Ig[i][hour],Id[i][hour],n,LST,gamma,beta))
    ITilted.append(ITmonth)
print(ITilted)

file = input("Enter name of File to save data: ")
array = np.array(ITilted)
np.save(file,array)

