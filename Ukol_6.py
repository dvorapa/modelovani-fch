### Skript pro výpočet pH vodného roztoku H3BO3
from sympy import *
m1, m2, m3, m4, m5, m6, m7 = symbols('m1 m2 m3 m4 m5 m6 m7', real=True)
x1, x2, x3, x4 = symbols('x1 x2 x3 x4', real=True)
zA, zK = symbols('zA zK', real=True)

## Nadefinované podmínky
T = 273.15 + 25 # K
c = 0.07 # M (mol/l)
cw = 997/18.016 # mol/l

## Tabulkové hodnoty disociačních koeficientů (Wikipedie, 2018)
pK_H2BO3_ = 9.24
pK_HBO3_2 = 12.4
pK_BO3_3 = 13.3
pK_H2O = 14

## Mezivýpočet disociačních konstant K
K1 = 10**(-pK_H2BO3_)
K2 = 10**(-pK_HBO3_2)
K3 = 10**(-pK_BO3_3)
Kw = 10**(-pK_H2O)

## Iterace
x_old = [6.32200549982070e-6, 3.95048911992951e-13, 3.10775588446771e-21, 8.68626027795196e-8]
while True:
    # Bilance
    m1 = c-x1 # H3BO3
    m2 = x1-x2 # H2BO3 -
    m3 = x2-x3 # HBO3 2-
    m4 = x3 # BO3 3-
    m5 = x1+x2+x3+x4 # H +
    m6 = x4 # OH -
    m7 = cw-x4 # H2O

    # Bilanční rovnice
    ion = (m2*(1**2) + m3*(2**2) + m4*(3**2) + m5*(1**2) + m6*(1**2))/2

    # Debye-Hückelova rovnice
    gam = exp(-(zA*zK)*0.3915*((sqrt(ion)/(1+1.2*sqrt(ion)))+((2*log(1+1.2*sqrt(ion)))/1.2)))

    # Výpočet koeficientů gamma z původních x
    gam = gam.subs([(x1, x_old[0]), (x2, x_old[1]), (x3, x_old[2]), (x4, x_old[3])])
    print('gamma = ' + str(gam))
    gam11 = gam.subs([(zA, 1), (zK, 1)])
    gam21 = gam.subs([(zA, 2), (zK, 1)])
    gam31 = gam.subs([(zA, 3), (zK, 1)])

    # Výpočet nových x
    rovn1 = m5*m2*gam11*gam11/m1 - K1
    rovn2 = m5*m3*gam21/m2 - K2
    rovn3 = m5*m4*gam11*gam31/(m3*gam21) - K3
    rovnw = m5*m6*gam11*gam11/m7 - Kw
    x_new = nsolve([rovn1, rovn2, rovn3, rovnw], [x1, x2, x3, x4], x_old)

    # Porovnání původních a nových x
    if x_old == [x for x in x_new]:
        break
    else:
        x_old = [x for x in x_new]

## Výsledná x
for i, x in enumerate(x_old, 1):
    print('x' + str(i) + ' = ' + str(x) + ' mol/l')

## Výsledná m
m_obc = [m1, m2, m3, m4, m5, m6, m7]
for i, m in enumerate(m_obc, 1):
    print('m' + str(i) + ' = ' + str(m.subs([(x1, x_old[0]), (x2, x_old[1]), (x3, x_old[2]), (x4, x_old[3])])) + ' mol/l')
