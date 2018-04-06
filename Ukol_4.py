### Skript pro výpočet změny entalpie reakce C(s) + CO2(g) = 2 CO(g)
from sympy import *
T = Symbol('T')
A = Symbol('A')
B = Symbol('B')
C = Symbol('C')
D = Symbol('D')
E = Symbol('E')
x = Symbol('x')
Cpx = Symbol('Cpx')

## Nadefinované podmínky
T1 = 298.15 # K
T2 = 798.15 # K

## Stechiometrie
r = 2
s = 0
a = 1
b = 1

## Tabelované hodnoty (Chase, 1998; NIST, 2018)
# slučovací entalpie
del_Hf_r = -110.53 * 10**3 # J/mol
del_Hf_s = 0 # J/mol
del_Hf_a = 0 # J/mol
del_Hf_b = -393.52 * 10**3 # J/mol

# tepelné kapacity
Cp = A + B*T + C*(T**2) + D*(T**3) + E/(T**2)
Ar = 25.56759
Br = 6.096130
Cr = 4.054656
Dr = -2.671301
Er = 0.131021
Cpr = Cp.subs([(A, Ar), (B, Br), (C, Cr), (D, Dr), (E, Er), (T, T/1000)])
Cps = 0
Cpa = 8.11 # J/(mol*K) (Sheindlich et al., 1972; pro grafit v rozmezí teplot od 273 do 1000 K)
Ab = 24.99735
Bb = 55.18696
Cb = -33.69137
Db = 7.948387
Eb = -0.136638
Cpb = Cp.subs([(A, Ab), (B, Bb), (C, Cb), (D, Db), (E, Eb), (T, T/1000)])

## Mezivýpočty
# Výpočet st. reakční entalpie
del_Hr_ = r*del_Hf_r + s*del_Hf_s - a*del_Hf_a - b*del_Hf_b

# Výpočet entalpické bilance
funkce_ = x*Cpx
del_Hr = del_Hr_\
+ integrate(funkce_.subs([(x, a), (Cpx, Cpa)]), (T, T2, T1))\
+ integrate(funkce_.subs([(x, b), (Cpx, Cpb)]), (T, T2, T1))\
+ integrate(funkce_.subs([(x, r), (Cpx, Cpr)]), (T, T1, T2))\
+ integrate(funkce_.subs([(x, s), (Cpx, Cps)]), (T, T1, T2))

# Výsledná entalpie reakce
print('del_Hr = ' + str(del_Hr) + ' J/mol\n')
