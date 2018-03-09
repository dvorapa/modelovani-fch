### Skript pro výpočet změny U, H, S, F a G pro 1 mol chlóru
from sympy import *
t = Symbol('t')

## Nadefinované podmínky
T1 = 298 # K
p1 = 0.1*10**6 # Pa

T2 = 500 # K
p2 = 10*10**6 # Pa

## Tabelované hodnoty a funkce (webbooknistgov/cgi/cbookcgi?ID=C7782505&Mask=1)
# (použitelné pouze v rozmezí 298-1000 K)
# Konstanty:
A = 33.05060
B = 12.22940
C = -12.06510
D = 4.385330
E = -0.159494
G = 259.0290

# Funkce
# St = @(t) A*log(t) + B*t + (C*t**2)/2 + (D*t**3)/3 -E/(2*t**2) + G
Cpot = A + B*t + C*t**2 + D*t**3 - E/(t**2)
# S = @(t) Sot(t/1000) # J/(mol*K)
Cpo = Cpot.subs(t, t/1000) # J/(mol*K)

# Hodnoty
So = 223.08 # J/(mol*K)
S1 = So

## Mezivýpočty
# Přepočet Cp na Cv
R = 8.314
Cvo = Cpo - R

# Výpočet V1 a V2
V1 = R*T1/p1
V2 = R*T2/p2

## Výpočet delS
funkce_ = Cpo/t
S2 = S1 + 0 + integrate(funkce_, (t, T1, T2)) - R*log(p2/p1) + 0 # 0 znamenají nulové integrály R/p* - R/p, protože počítáme plyn ve stavu IP
delS = S2 - S1

## Výpočet delU
funkce_ = Cvo
delU = 0 + integrate(funkce_, (t, T1, T2)) + 0

## Výpočet delH
funkce_ = Cpo
delH = 0 + integrate(funkce_, (t, T1, T2)) + 0

## Výpočet delG
delG = delH - T2*S2 + T1*S1

## Výpočet delF
delF = delU - T2*S2 + T1*S1

# Výsledné výsledky
print('delS = ' + str(delS.evalf()) + ' J/(mol*K)\n')
print('delU = ' + str(delU.evalf()) + ' J\n')
print('delH = ' + str(delH.evalf()) + ' J\n')
print('delG = ' + str(delG.evalf()) + ' J\n')
print('delF = ' + str(delF.evalf()) + ' J\n')
