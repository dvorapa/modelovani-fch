### Skript pro výpočet aktivitních koeficientů pomocí metod NRTL a Van Laar
from sympy import *
import numpy as np

from Semestralka_vystupy import *


## Parametry a hodnoty společné všem metodám
x1, T = symbols('x1 T')
x2 = 1 - x1

chem = ('A', 'T', 'C')
mix1 = ('A', 'T')
mix2 = ('A', 'C')
mixy = (mix1, mix2)

# Známé konstanty
R = 8.314
T0 = 273.15

M = (120.1485, 131.388, 84.1595) # Zdroj databáze https://webbook.nist.gov/chemistry/

## Experimentální data dle Anily et al. (https://doi.org/10.1016/j.molliq.2014.12.026)
ro = (1019.86, 1451.32, 773.85)
ro = tuple(i / 1000 for i in ro)

x1__ = (0.0000, 0.0777, 0.1441, 0.2017, 0.2519, 0.2963, 0.3356, 0.3708, 0.4025, 0.4311, 0.5310, 0.5601, 0.5927, 0.6294, 0.6708, 0.7181, 0.7725, 0.8359, 0.9106, 1.0000)

# Experimentální hodnoty T
T__mix1 = (358.10, 361.46, 364.24, 366.12, 368.42, 369.25, 370.12, 371.10, 372.14, 374.16, 379.16, 380.25, 383.15, 386.54, 389.66, 395.16, 405.16, 420.56, 437.66, 472.30)
T__mix2 = (351.60, 354.16, 355.16, 356.16, 357.16, 358.16, 359.16, 360.16, 362.16, 363.16, 369.66, 373.16, 376.16, 380.16, 388.16, 421.55, 442.45, 445.45, 460.34, 472.30)
T__ = (T__mix1, T__mix2)


## Metoda Van Laar
metoda = 'Van Laar'
print(metoda)

# Výpočet Van Laar koeficientů A a B
AB = []
for (ro_VL, M_VL) in zip(ro, M):
  AB.append(M_VL/ro_VL)

A, B = symbols('A B')

# Příprava rovnic pro aktivitní koeficienty a Gibbsovu energii
rovn1 = exp(B / (R * (T - T0) * ((1 + ((B * x1)/(A * x2)))**2)))
rovn2 = exp(A / (R * (T - T0) * ((1 + ((A * x2)/(B * x1)))**2)))
rovnQ = (A * B * x1 * x2) / (A * x2 + B * x1)

# Výpočet aktivitních koeficientů γ1 a γ2 pro jednotlivé směsi
vysl = []
for (mix, T__mixu) in zip(mixy, T__):
  # Pro každou směs
  print(' + '.join(mix))
  vysl_mix = [[], [], [], [], []]
  for (x1_, T_) in zip(x1__, T__mixu):
    # Pro každou experimentální hodnotu x1 a T
    vgam1 = rovn1.subs(((x1, x1_), (T, T_), (A, AB[chem.index(mix[0])]), (B, AB[chem.index(mix[1])])))
    vgam2 = rovn2.subs(((x1, x1_), (T, T_), (A, AB[chem.index(mix[0])]), (B, AB[chem.index(mix[1])])))
    vQ = rovnQ.subs(((x1, x1_), (A, AB[chem.index(mix[0])]), (B, AB[chem.index(mix[1])])))
    vysl_mix[0].append(x1_)
    vysl_mix[1].append(vgam1)
    vysl_mix[2].append(vgam2)
    vysl_mix[3].append(vQ)
    vysl_mix[4].append(log(vgam1/vgam2))
  vysl.append(vysl_mix)

  # Výpis výsledných hodnot
  tabulka_vysledku(T__mixu, vysl_mix)

  # Graf aktivitních koeficientů směsi
  graf_gammy(metoda, mix, vysl_mix)
  
  # Graf logaritmu poměru aktivitních koeficientů směsi
  graf_logaritmu(metoda, mix, vysl_mix)

# Graf Gibbsových energií směsí
graf_Q(metoda, vysl)


## Metoda NRTL
metoda = 'NRTL'
print('\n' + metoda)

# Experimentální data a konstanty pro NRTL dle Anily et al. (https://doi.org/10.1016/j.molliq.2014.12.026)
# Poznámka: U směsi A + C autoři uvádí konstanty zhruba 6,4× menší, než by odpovídalo skutečnosti.
b__mix1 = (51.23, 49.88)
b__mix2 = (333.72, 390.89)
b__mix2 = tuple(b / 6.4 for b in b__mix2)
b__ = (b__mix1, b__mix2)

alpha = 0.3

# Výpočet NRTL koeficientů tau1 a tau2 a funkcí g12 a g21
b12, b21 = symbols('b12 b21')

tau12 = b12 / (R * (T - T0))
tau21 = b21 / (R * (T - T0))

g12 = exp(-alpha * tau12)
g21 = exp(-alpha * tau21)

# Příprava rovnic pro aktivitní koeficienty a Gibbsovu energii
rovn1 = exp((x2**2) * (((tau12 * g12**2) / ((x1 + x2 * g12)**2)) + ((tau21 * g21) / ((x2 + x1 * g21)**2))))
rovn2 = exp((x1**2) * (((tau21 * g21**2) / ((x2 + x1 * g21)**2)) + ((tau12 * g12) / ((x1 + x2 * g12)**2))))
rovnQ = (x1 * x2 * (((tau21 * g21) / (x2 + x1 * g21)) + ((tau12 * g12) / (x1 + x2 * g12)))) * (R * (T - T0))

# Výpočet aktivitních koeficientů γ1 a γ2 pro jednotlivé směsi
vysl = []
for (mix, T__mixu, b__mixu) in zip(mixy, T__, b__):
  # Pro každou směs
  print(' + '.join(mix))
  vysl_mix = [[], [], [], [], []]
  for (x1_, T_) in zip(x1__, T__mixu):
    # Pro každou experimentální hodnotu x1 a T
    vgam1 = rovn1.subs(((x1, x1_), (b12, b__mixu[0]), (b21, b__mixu[1]), (T, T_)))
    vgam2 = rovn2.subs(((x1, x1_), (b12, b__mixu[0]), (b21, b__mixu[1]), (T, T_)))
    vQ = rovnQ.subs(((x1, x1_), (b12, b__mixu[0]), (b21, b__mixu[1]), (T, T_)))
    vysl_mix[0].append(x1_)
    vysl_mix[1].append(vgam1)
    vysl_mix[2].append(vgam2)
    vysl_mix[3].append(vQ)
    vysl_mix[4].append(log(vgam1/vgam2))
  vysl.append(vysl_mix)

  # Výpis výsledných hodnot
  tabulka_vysledku(T__mixu, vysl_mix)

  # Graf aktivitních koeficientů směsi
  graf_gammy(metoda, mix, vysl_mix)
  
  # Graf logaritmu poměru aktivitních koeficientů směsi
  graf_logaritmu(metoda, mix, vysl_mix)

# Graf Gibbsových energií směsí
graf_Q(metoda, vysl)