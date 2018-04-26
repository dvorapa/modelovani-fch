### Skript pro výpočet rovnovážného složení reakce C(s) + CO2(g) = 2 CO(g)
from sympy import *
zet = Symbol('zet')
R = 8.314


gam = 1


## Nadefinované podmínky
T = 598.15 # K
p = 30 * 10**6 # Pa

## Tabulkové hodnoty (Chase, 1998; NIST, 2018)
delHf_CO = -110.53 * 10**3 # J/mol
S_CO = -197.66 # J/(mol*K)

delHf_CO2 = -393.52 * 10**3 # J/mol
S_CO2 = -213.79 # J/(mol*K)

S_C = 5.6 # J/(mol*K) průměr 10 hodnot

S_O2 = 205.15 # J/(mol*K)

## Mezivýpočet změny G v reakci
delGf_C = 0
delGf_CO2 = delHf_CO2 - T*S_CO2 + T*(S_C + S_O2)
delGf_CO = delHf_CO - T*S_CO + T*(S_C + (S_O2/2))

delGr = 2*delGf_CO - delGf_CO2 - delGf_C

## Mezivýpočet K
K = exp(-delGr/(R*T))

## Mezivýpočet rozsahu reakce
rovn = zet/((1-zet)**2) - K*gam
koreny = solve(rovn, zet)
for vysl in koreny:
    if vysl >= 0 and vysl <= 1:
        zet = vysl
        break

## Výsledky jednotlivých mezivýpočtů a výsledek
print('delGr = ' + str(N(delGr)) + ' J/mol')
print('K = ' + str(N(K)))
print('zeta = ' + str(zet))
