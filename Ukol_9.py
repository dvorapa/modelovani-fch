### Skript pro výpočet Langmuirových adsorpčních izoterm
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
R = 8.314

## Nadefinované podmínky
V_smesi = 0.25 # l
m_adsorbentu = 3 # g

# Koncentrace
c_pred = [0.75, 0.85, 1.0, 1.1] # mol/l
c_po = [0.35, 0.42, 0.54, 0.61] # mol/l

## Mezivýpočet změny koncentrace a množství rozp. látky
a = []
c = []
i = 0
while i < len(c_pred):
    del_c = c_pred[i] - c_po[i]
    del_c *= V_smesi # mol
    an = del_c/m_adsorbentu # mol/g
    a.append(an)

    c.append(c_po[i] * V_smesi) # mol

    i += 1

## Příprava vektorů x a y
x = []
y = []
i = 0
while i < len(a):
    x.append(-a[i])
    y.append(a[i]/c[i])
    i += 1

## Výpočet lineární regrese
x = np.array(x, dtype=float)
y = np.array(y, dtype=float)
vysl = np.polyfit(x, y, 1)

## Tisk grafu
plt.plot(x, y)

yfit = []
i=0
while i < len(x):
    yfit.append(vysl[0]*x[i] + vysl[1])
    i += 1
plt.plot(x, yfit)

plt.xlabel('−a [−mol/g]')
plt.ylabel('a/c [mol^-1]')
plt.title('Graf korelace regresní křivky (červeně) a dat (modře)')
plt.savefig('Ukol_9.png')

## Výsledný řád reakce a rychlostní konstanta
print('b = ' + str(N(vysl[0])) + ' mol^-1')
print('am = ' + str(N(vysl[1]/vysl[0])) + ' mol/g')
