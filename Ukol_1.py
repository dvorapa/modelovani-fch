### Skript pro výpočet změny S, U, H, G a F jednoho molu plynného chlóru
from sympy import *
T = Symbol('T')
R = 8.314

## Nadefinované podmínky
T1 = 298 # K
p1 = 0.1 * 10**6 # Pa

T2 = 500 # K
p2 = 10 * 10**6 # Pa

## Tabulkové hodnoty (Chase, 1998; NIST, 2018)
S1 = 223.08 # J/(mol*K)

A = 33.0506
B = 12.2294
C = -12.0561
D = 4.38533
E = -0.159494

## Mezivýpočet tepelných kapacit
Cp = A + B*T + C*(T**2) + D*(T**3) + E/(T**2)
Cp = Cp.subs(T, T/1000)
Cv = Cp - R

## Výpočet změny entropie
S2 = S1 + integrate(Cp/T, (T, T1, T2)) + R * log(p2/p1)
delS = S2 - S1

## Výpočet změny vnitřní energie
delU = integrate(Cv, (T, T1, T2))

## Výpočet změny entalpie
delH = integrate(Cp, (T, T1, T2))

## Výpočet změny Gibbsovy volné energie
delG = delH - T2*S2 + T1*S1

## Výpočet změny Helmholtzovy volné energie
delF = delU - T2*S2 + T1*S1

## Výsledné změny termodynamických veličin
print('delS = ' + str(N(delS)) + ' J/(mol*K)')
print('delU = ' + str(N(delU)) + ' J/mol')
print('delH = ' + str(N(delH)) + ' J/mol')
print('delG = ' + str(N(delG)) + ' J/mol')
print('delF = ' + str(N(delF)) + ' J/mol')
