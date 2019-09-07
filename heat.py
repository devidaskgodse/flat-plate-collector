from math import * # sin, cos, acos, asin, pi, exp, tanh, tan
from datetime import datetime


d = pi / 180


# Location Specification
phi = d * float(input("Enter latitude in degrees: "))
''' longitudes to be taken as +ve for east and -ve for west '''
longitude = float(input("Enter longitude in degrees: "))
lsm = float(input("Enter standard time meridian in degrees: ")) 

# Time Specification
year,month,day = input("Enter day in format YYYY-MM-DD: ").split('-')
n = int(format(datetime(int(year),int(month),int(day)),'%j'))
LST = float(input("Enter time of the day in hours: "))

# Solar radiation incident on ground in W/m^2
Ig = float(input("Enter total irradiration on horizontal in W/m^2: ")) # global radiation
Id = float(input("Enter diffuse irradiration on horizontal in W/m^2: ")) # diffuse radiation

# collector orientation
gamma = d * float(input("Enter surface azimuth angle in degrees: ")) # Surface azimuthal angle
beta = d * (float(input("Enter tilt angle of surface in degrees: "))) # Tilt of the surface


Ib = Ig - Id # beam radiation incident on ground in W/m^2

delta = d * (23.45 * sin(360 * (284.0 + n) * d / 365.0)) # declination

B = d * (360.0 * (n - 1) / 365.0)
ET = (9.87 * sin(2.0*B) - 7.53 * cos(B) - 1.5 * sin(B)) / 60.0 # Equation of time
ST = LST + ET + (lsm - longitude) / 15.0
omega = d * (15 * (ST - 12.0)) # Hour angle

# cosine of angle of incidence is given by
costheta1 = sin(phi) * sin(delta) * cos(beta) + \
            + sin(phi) * cos(delta) * cos(gamma) * cos(omega) * sin(beta) \
            + cos(phi) * cos(delta) * cos(omega) * cos(beta)\
           - cos(phi) * sin(delta) * cos(gamma) * sin(beta)\
           + cos(delta) * sin(gamma) * sin(omega) * sin(beta)
# cosine of zenith angle is given by
costhetaz = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega)

ro = 0.2 # reflection coefficient for grass
Rb = costheta1 / costhetaz
Rd = cos(beta * d / 2.0)**2.0
Rr = ro * sin(beta * d / 2.0)**2.0

## Solar radiation incident on tilted surface
It = Ib * Rb + Id * Rd + Ig * Rr 

print('Solar radiation incident on tilted surface in W/m^2 -  ', It)


## calculation of incident solar flux absorbed in absorber plate
eta = float(input("Enter refractive index of cover w. r. t. air: ")) # Refractive Index of cover
M = float(input("Enter no. of transperant covers on the surface: ")) # No. of top covers
dc = float(input('Enter thickness of one transperant cover in metres: '))
theta1 = acos(costheta1) # Angle of incidence for beam radiation
theta2 = asin(sin(theta1) / eta) # Angle of refraction for beam radiation
K = 15 # Extinction coefficient for glass in per metre

# Reflectivity of two components of polarization of beam radiation
robI = (sin(theta2 - theta1))**2 / (sin(theta2 + theta1))**2 
robII = (tan(theta2 - theta1))**2 / (tan(theta2 + theta1))**2
# Transmitivity of two components of polarization
TrbI = (1 - robI) / (1 + (2*M - 1) * robI)
TrbII = (1 - robII) / (1 + (2*M - 1) * robII)
Trb = (TrbI + TrbII) / 2 # Transmittivity considering only reflection and refraction

Tab = exp(-1 * M * K * dc / cos(theta2)) # Transmittivity considering only absorption


theta1 = d * 60 # Angle of incidence for diffuse radiation
theta2 = asin(sin(theta1) / eta) # Angle of refraction for diffuse radiation
# Reflectivity of two components of polarization of diffuse radiation
rodI = (sin(theta2 - theta1))**2 / (sin(theta2 + theta1))**2 
rodII = (tan(theta2 - theta1))**2 / (tan(theta2 + theta1))**2
# Transmitivity of two components of polarization
TrdI = (1 - rodI) / (1 + (2*M - 1) * rodI)
TrdII = (1 - rodII) / (1 + (2*M - 1) * rodII)
Trd = (TrdI + TrdII) / 2 # Transmittivity considering only reflection and refraction

Tad = exp(-1 * M * K * dc / cos(theta2)) # Transmittivity considering only absorption

rodee = Tad * (1 - Trd)

alpha = float(input("Enter plate absorptivity for absorber plate: ")) # plate absorptivity

# Transmittivity absorptivity product for beam radiation
TAb = (Trb * Tab * alpha) / (1 - ((1-alpha) * rodee))
# Transmittivity absorptivity product for diffuse radiation
TAd = (Trd * Tad * alpha) / (1 - ((1-alpha) * rodee))


## Incident solar flux absorbed in absorber plate
S = Ib * Rb * TAb + (Id * Rd + Ig * Rr) * TAd
print('Incident solar flux absorbed in absorber plate in W/m^2 -  ', S)


# convective heat transfer coefficient
wind = float(input('Enter wind speed in m/s: ')) # wind speed in m/s
hw = 8.55 + 2.56 * wind  # W/m^2-K


## Top loss coefficient calculation
L = float(input('Enter spacing between cover and plate in metres: ')) # metres
C = 204.429 * (cos(beta)**0.252) / L**0.24
Ta = float(input('Enter surrounding temperature in kelvin: ')) # Kelvin
Tp = float(input('Enter average plate temperature in kelvin: ')) # Kelvin
Ec = float(input('Enter emissivity of covers for long wavelength radiation: '))
Ep = float(input('Enter emmisivity of absorber plate for long wavelength radiation: '))
f = ((9/hw) - (30/hw**2)) * (Ta / 316.9) * (1 + 0.091 * M)
sigma = 5.67 * 10**(-8) # stephan boltzman constant in W/m^2-K^4
Ut1 = (M / ((C / Tp) * ((Tp - Ta) / (M + f))**0.252)) + (1 / hw)
Ut2n = sigma * (Tp**2 + Ta**2) * (Tp + Ta)
Ut2d = (1 / (Ep + 0.0425 * M * (1-Ep))) + ((2 * M + f - 1) / Ec) - M
Ut = Ut1**(-1) + (Ut2n / Ut2d) # W/m^2-K
print('Top loss coefficient in W/m^2-K -  ', Ut)

## Bottom loss coefficient calculation
Ki = float(input("Enter thermal conductivity of insulator in W/m-K: ")) # W/m-K
db = float(input("Enter Thickness of insulator at the bottom in m: ")) # in metres
Ub = Ki / db # W/m^2-K
print('Bottom loss coefficient in W/m^2-K- ', Ub)

## Side Loss Coefficient calculation
L1 = float(input("Enter length of absorber plate in metres: ")) # length in metres
L2 = float(input("Enter breadth of absorber plate in metres: ")) # breadth in metres
L3 = M * L + db # height in metres
ds = float(input('Enter thickness of insultor at side in metres: ')) # metres
Us = ((L1 + L2) * L3 * Ki) / (L1 * L2 * ds) # W/m^2-K
print('Side loss coefficient in W/m^2-K - ', Us)

### Overall loss coefficient
Ul = Ut + Ub + Us # W/m^2-K

## plate effectiveness calculation
N = float(input('Enter no. of tubes in collector: '))
Do = float(input('Enter outer diameter of tube in metres: ')) # metres
Di = float(input('Enter inner diameter of tube in metres: ')) # metres
Kp = float(input('Enter thermal conductivity of absorber plate in W/m-K: ')) # W/m-K
dp = float(input('Enter thickness of absorber plate in metres: ')) # metres
hf = float(input('Enter heat transfer coefficient between fluid and tube in W/m^2-K: ')) # W/m^2-K
W = L2 / N # pitch of absorber plate
m = (Ul / (Kp * dp))**(1/2.0)
phy = (tanh(m * (W - Do) / 2)) / (m * (W - Do) / 2)
print("Plate effectiveness -  ", phy)

## Collector efficiency factor calculation
# Assuming the resistance by adhesive is negligible
fdash = ((W/(Do + phy * (W - Do))) + W*Ul/(pi * Di * hf))**(-1)
print('Collector efficiency factor - ', fdash)

## collector heat removal factor
flowrate = float(input('Enter flow rate of fluid in kg/s: ')) # fluid flow rate in kg/s
Cp = float(input('Enter specific heat of fluid in J/kg-K: ')) # J/kg-K
Ap = L1 * L2 # Area of absorber plate in m^2
FR = flowrate * Cp * (1 - exp(-1 * fdash * Ul * Ap / (flowrate * Cp))) / (Ul * Ap)
print('Collector heat removal factor value - ', FR)

### Useful heat gain rate of collector
Tfi = float(input('Enter inlet temperature of fluid in kelvin: ')) # kelvin
Qu = FR * Ap * (S - Ul * (Tfi - Ta)) # Watt
print('Useful heat gain rate - ', Qu)

### Outlet fluid temperature
Tfo = (Qu / (flowrate * Cp)) + Tfi 
print('Outlet temperature of fluid in Kelvin -  ', Tfo) # Kelvin

#### Instantaneous efficiency of collector 
Ac = Ap * 1.15 # Collector gross area
efficiency = Qu / (Ac * It)
print('Instantaneous efficiency of collector - ', efficiency)




