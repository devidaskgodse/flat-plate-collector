# this file is for making a numerical method to find the value of mean-plate temperature
#Initial approximation to overall loss coefficient
import numpy as np
from datetime import datetime
from math import *

# constants
sigma = 5.67 * 10**-8 # W/m^2-K^4
Cp = 4184 #J/kg-K for water

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




''' Process'''

Ib = Ig - Id # beam radiation incident on ground in W/m^2

delta = d * (23.45 * sin(360 * (284.0 + n) * d / 365.0)) # declination

B = d * (360.0 * (n - 1) / 365.0)
ET = (9.87 * sin(2.0*B) - 7.53 * cos(B) - 1.5 * sin(B)) / 60.0 # Equation of time
ST = LST + ET + (lsm - longitude) / 15.0
omega = d * (15 * (ST - 12.0)) # Hour angle

costheta1 = sin(phi) * sin(delta) * cos(beta) + \
            + sin(phi) * cos(delta) * cos(gamma) * cos(omega) * sin(beta) \
            + cos(phi) * cos(delta) * cos(omega) * cos(beta)\
           - cos(phi) * sin(delta) * cos(gamma) * sin(beta)\
           + cos(delta) * sin(gamma) * sin(omega) * sin(beta)
costhetaz = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega)

ro = 0.2
Rb = costheta1 / costhetaz
Rd = cos(beta * d / 2.0)**2.0
Rr = ro * sin(beta * d / 2.0)**2.0

''' Solar radiation incident on tilted surface'''
It = Ib * Rb + Id * Rd + Ig * Rr 

print('Solar radiation incident on tilted surface in W/m^2 -  ', It)


''' For calculation of incident solar flux absorbed in absorber plate'''

eta = float(input("Enter refractive index of cover w. r. t. air: ")) # Refractive Index
M = float(input("Enter no. of transperant covers on the surface: ")) # No. of top covers
dc = float(input('Enter thickness of one transperant cover in metres: '))
theta1 = acos(costheta1) # Angle of incidence for beam radiation
theta2 = asin(sin(theta1) / eta) # Angle of refraction for beam radiation

K = 19 # Extinction coefficient for glass in per metre. #Generally 15 per metre.

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


# Incident solar flux absorbed in absorber plate
S = Ib * Rb * TAb + (Id * Rd + Ig * Rr) * TAd
print('Incident solar flux absorbed in absorber plate in W/m^2 -  ', S)




def Pr(T):
  Pr = 0.7418 - 0.0001373*T #Prandtl number
  return Pr

def v(T):
  v = (-17.2813 + 0.1092 * T)*10**(-6)
  return v

def k(T):
  k = 0.004245 + 0.0000742245*T # thermal conductivity
  return k

def effRaL(Tp,Tc,beta,L):
  mean = (Tp + Tc)/2
  RaL = 9.81 * (Tp - Tc)*Pr(mean)*L**3/(mean * v(mean)**2) # effective rayleigh number
  return RaL*np.cos(beta)

def NuL(effRaL):
  if effRaL <= 1708:
    NuL = 1.0
  elif 1708 < effRaL <= 5900:
    NuL = 1 + 1.446*(1-1708.0/effRaL)
  elif 5900 < effRaL <= 9.23e4:
    NuL = 0.229*(effRaL**0.252)
  elif 9.23e4 < effRaL <= 1.0e6:
    NuL = 0.157*(effRaL**0.285)
  return NuL




# Bottom loss coefficient calculation
Ki = float(input("Enter thermal conductivity of insulator in W/m-K: ")) # W/m-K
db = float(input("Enter Thickness of insulator at the bottom in m: ")) # in metres
Ub = Ki / db # W/m^2-K

# Side Loss Coefficient calculation
L1 = float(input("Enter length of absorber plate in metres: ")) # length in metres
L2 = float(input("Enter breadth of absorber plate in metres: ")) # breadth in metres
L = float(input('Enter length of spacing between covers in metres: ')) # metres
L3 = M * L + db # height in metres
Ap = L1 * L2 # Area of absorber plate in m^2
ds = float(input('Enter thickness of insultor at side in metres: ')) # metres
Us = ((L1 + L2) * L3 * Ki) / (L1 * L2 * ds) # W/m^2-K

Ec = float(input('Enter emissivity of covers for long wavelength radiation: '))
Ep = float(input('Enter emmisivity of absorber plate for long wavelength radiation: '))
Kp = float(input('Enter thermal conductivity of absorber plate in W/m-K: ')) # W/m-K
dp = float(input('Enter thickness of absorber plate in metres: ')) # metres
hf = float(input('Enter heat transfer coefficient between fluid and tube in W/m^2-K: '))
flowrate = float(input('Enter flow rate of fluid in kg/s: ')) # fluid flow rate in kg/s

Tfi = float(input('Enter inlet temperature of fluid in kelvin: '))
Ta = float(input('Enter surrounding temperature in kelvin: '))
Do = float(input('Enter outer diameter of tube in metres: ')) # metres
Di = float(input('Enter inner diameter of tube in metres: ')) # metres
W = float(input('Enter tube center-to-center distance in metres: ')) # pitch of absorber plate

wind = float(input('Enter wind speed in m/s: '))
hw = 8.55 + 2.56 * wind

Ul1 = 0 # W/m^2-K
Ul2 = 4

while np.abs(Ul1 - Ul2) > 0.05:
  m = (Ul2/(Kp*dp))**0.5
  phy = (tanh(m * (W - Do) / 2)) / (m * (W - Do) / 2)
  
  fdash = ((W/(Do + phy * (W - Do))) + W*Ul2/(pi * Di * hf))**(-1)
  
  FR = flowrate * Cp * (1 - exp(-1 * fdash * Ul2 * Ap / (flowrate * Cp))) / (Ul2 * Ap)
  
  Qu = FR * Ap * (S - Ul2 * (Tfi - Ta))
  Ql = S*Ap - Qu

  Tp = (Ql/(Ul2*Ap)) + Ta
  Tc = 307 # Kelvin
  hp = NuL(effRaL(Tp,Tc,beta,L))*k((Tp+Tc)/2)/L
  
  Qt1 = hp * (Tp - Tc) + sigma*(Tp**4 - Tc**4)/(1/Ec + 1/ Ep - 1)
  Qt2 = hw * (Tc - Ta) + sigma*Ec*(Tc**4 - (Ta - 6)**4)
    
  Qt = (Qt1 + Qt2)/2
  Ut = Qt / (Tp - Ta)  

  Ul1 = Ul2
  Ul2 = Ut + Ub + Us
  print('Overall Loss coefficient',Ul2)


# Useful heat gain rate of collector
Qu = FR * Ap * (S - Ul2 * (Tfi - Ta)) # Watt
print('Useful heat gain rate - ', Qu)

# Outlet fluid temperature
Tfo = (Qu / (flowrate * Cp)) + Tfi 
print('Outlet temperature of fluid in Kelvin -  ', Tfo) # Kelvin

# Instantaneous efficiency of collector
length = float(input('Enter length of collector in metres: '))
breadth = float(input('Enter breadth of collector in metres: ')) 
Ac = length * breadth # Collector gross area
efficiency = Qu / (Ac * It)
print('Instantaneous efficiency of collector - ', efficiency)

