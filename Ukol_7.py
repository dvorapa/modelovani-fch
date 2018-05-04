### Skript pro výpočet množství rozpuštěného argonu ve vodě
from sympy import *
R = 8.314

## Nadefinované podmínky
T = 298.15 # K
p = 101325 # Pa
M_Ar = 39.944 # g/mol
ro_H2O = 0.997 # g/cm3

## Tabulkové hodnoty
# Tabulky pro tepelnou techniku: Hašek, 1980
V_zlomek = 0.93 * 10**-2
m_zlomek = 1.286 * 10**-2
# Rettich et al., 1992
kH = 1.4 * 10**(-5) # mol/(Pa*m^3)

# Převod jednotek
kH = 1/kH # Pa*m^3/mol

## Dosazení do rovnice
pi = V_zlomek*p
c = pi/kH # mol/m^3

## Jednotková konverze
c = c/1000 # mol/l
c *= M_Ar*1000 # mg/l

## Přepočet tabulkové koncentrace
c_tab = 260.63/20.381 # µmol/kg (Hamme a Emerson, 2004)
c_tab *= 10**-6 # mol/kg
c_tab *= ro_H2O # mol/l
c_tab *= M_Ar*1000 # mg/l

## Výsledné změny termodynamických veličin
print('c = ' + str(N(c)) + ' mg/l')
print('c = ' + str(N(c_tab)) + ' mg/l (Hamme a Emerson, 2004)')
