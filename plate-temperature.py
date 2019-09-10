# this file is for making a numerical method to find the value of mean-plate temperature
#Initial approximation to overall loss coefficient
import numpy as np
Ul1 = 0 # W/m^2-K
Ul2 = 4

def effRaL():
  RaL = 9.81 * (Tp - Tc)*Pr(T)

def NuL(effRaL):
  if EffRaNumber <= 1708:
    NuL = 1.0
  elif 1708 < effRaL <= 5900:
    NuL = 1 + 1.446*(1-1708.0/effRaL)
  elif 5900 < effRaL <= 9.23e4:
    NuL = 0.229*(effRaL**0.252)
  elif 9.23e4 < effRaL <= 1.0e6:
    NuL = 0.157*(effRaL**0.285)
  return NuL


while np.abs(Ul1 - Ul2) > 0.05:
  m = (Ul2/(Kp*dp))**0.5
  phy = (tanh(m * (W - Do) / 2)) / (m * (W - Do) / 2)
  
  fdash = ((W/(Do + phy * (W - Do))) + W*Ul2/(pi * Di * hf))**(-1)
  
  FR = flowrate * Cp * (1 - exp(-1 * fdash * Ul * Ap / (flowrate * Cp))) / (Ul2 * Ap)
  
  Qu = FR * Ap * (S - Ul2 * (Tfi - Ta))
  Ql = S*Ap - Qu
  Tp = (Ql/(Ul2*Ap)) + Ta
