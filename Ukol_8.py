### Skript pro výpočet řádu reakce a rychlostní konstanty
from sympy import *
import numpy as np
R = 8.314

## Nadefinované podmínky
T = 293.15 + 50 # K

# Koncentrace
t = [0, 40, 70, 220, 320, 404, 642] # s
c = [0.04, 0.037, 0.035, 0.028, 0.025, 0.023, 0.019] # mol/l

## Mezivýpočet rychlosti změny
r = []
ci = []
i = 1
while i < len(t):
    r.append(log((c[i-1]-c[i])/(t[i]-t[i-1])))
    ci.append(log((c[i]+c[i-1])/2))
    i += 1

## Výpočet lineární regrese
x = np.array(ci, dtype=float)
y = np.array(r, dtype=float)
vysl = np.polyfit(x, y, 1)

## Výsledný řád reakce a rychlostní konstanta
print('alpha = ' + str(N(vysl[0])))
print('k = ' + str(N(exp(vysl[1]))) + ' s^-1')
