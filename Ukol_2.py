### Skript pro výpočet hustoty dvojice kapalin kyselina octová - 1,2-dichloropropan 1:3
from sympy import *
Vm = Symbol('Vm')
Tk = Symbol('Tk')
pk = Symbol('pk')
Vk = Symbol('Vk')
R = 8.314

## Nadefinované podmínky
T = 273.15 + 25 # K
p = 101325 # Pa

# Kritické veličiny pro kys. octovou (Prausnitz et al., 2001; NIST, 2018)
Tk_oct = 592.7 # K
pk_oct = 57.9 * 10**5 # Pa
rok_oct = 5.84 # mol/l (Vandana a Teja, 1995; NIST, 2018)
Vk_oct = 1/(rok_oct*1000)

# Kritické veličiny pro 1,2-dichloropropan (Steele et al., 1997; NIST, 2018)
Tk_dcp = 578 # K
pk_dcp = 46.5 * 10**5 # Pa
rok_dcp = 3.452 # mol/l
Vk_dcp = 1/(rok_dcp*1000)

## Redlich-Kwong
a = (0.42748*(R**2)*(Tk**2.5))/pk
b = (0.08664*R*Tk)/pk
rovnice_ = ((R*T)/(Vm-b)) - (a/(sqrt(T)*Vm*(Vm+b))) - p

koreny_oct = solve(rovnice_.subs([(Tk, Tk_oct), (pk, pk_oct)]), Vm)
koreny_oct_re = []
for k in koreny_oct:
    if im(k) < 10**21:
        koreny_oct_re.append(re(k))
Vm_oct = min(koreny_oct_re)
ro_oct_rk = 1/(Vm_oct*(10**3))

koreny_dcp = solve(rovnice_.subs([(Tk, Tk_dcp), (pk, pk_dcp)]), Vm)
koreny_dcp_re = []
for k in koreny_dcp:
    if im(k) < 10**21:
        koreny_dcp_re.append(re(k))
Vm_dcp = min(koreny_dcp_re)
ro_dcp_rk = 1/(Vm_dcp*(10**3))

# Amagat
Vm_rk = (0.25*Vm_oct + 0.75*Vm_dcp)
ro_rk_a = 1/(Vm_rk*(10**3))

# Kay
Tk_ = 0.25*Tk_oct + 0.75*Tk_dcp
pk_ = 0.25*pk_oct + 0.75*pk_dcp

koreny_ = solve(rovnice_.subs([(Tk, Tk_), (pk, pk_)]), Vm)
koreny_re_ = []
for k in koreny_:
    if im(k) < 10**21:
        koreny_re_.append(re(k))
Vm_ = min(koreny_re_)
ro_rk_k = 1/(Vm_*(10**3))

# Joffe
K1 = (((0.25*Tk_oct)/sqrt(pk_oct)) + ((0.75*Tk_dcp)/sqrt(pk_dcp)))
K2 = 0
seznam_x = (0.25, 0.75)
seznam_Tk = (Tk_oct, Tk_dcp)
seznam_pk = (pk_oct, pk_dcp)
for i in range(2):
    for j in range(2):
        K2 += seznam_x[i]*seznam_x[j]*((((seznam_Tk[i]/seznam_pk[i])**(1/3))+((seznam_Tk[j]/seznam_pk[j])**(1/3)))**3)
K2 *= 0.125

Tk_ = (K1**2)/K2
pk_ = (K1/K2)**2

koreny_ = solve(rovnice_.subs([(Tk, Tk_), (pk, pk_)]), Vm)
koreny_re_ = []
for k in koreny_:
    if im(k) < 10**21:
        koreny_re_.append(re(k))
Vm_ = min(koreny_re_)
ro_rk_j = 1/(Vm_*(10**3))

## Rackett
rovnice_ = Vk*(((pk*Vk)/(R*Tk))**(((1-T)/Tk)**(2/7)))

Vm_oct = N(rovnice_.subs([(Tk, Tk_oct), (pk, pk_oct), (Vk, Vk_oct)]))
if im(Vm_oct) < 10**21:
    Vm_oct = re(Vm_oct)
ro_oct_r = 1/(Vm_oct*(10**3))

Vm_dcp = N(rovnice_.subs([(Tk, Tk_dcp), (pk, pk_dcp), (Vk, Vk_dcp)]))
if im(Vm_dcp) < 10**21:
    Vm_dcp = re(Vm_dcp)
ro_dcp_r = 1/(Vm_dcp*(10**3))

# Amagat
Vm_r = (0.25*Vm_oct + 0.75*Vm_dcp)
ro_r_a = 1/(Vm_r*(10**3))

# Kay
Tk_ = 0.25*Tk_oct + 0.75*Tk_dcp
pk_ = 0.25*pk_oct + 0.75*pk_dcp
Vk_ = 0.25*Vk_oct + 0.75*Vk_dcp

Vm_ = N(rovnice_.subs([(Tk, Tk_), (pk, pk_), (Vk, Vk_)]))
if im(Vm_) < 10**21:
    Vm_ = re(Vm_)
ro_r_k = 1/(Vm_*(10**3))

# Joffe
K1 = (((0.25*Tk_oct)/sqrt(pk_oct)) + ((0.75*Tk_dcp)/sqrt(pk_dcp)))
K2 = 0
seznam_x = (0.25, 0.75)
seznam_Tk = (Tk_oct, Tk_dcp)
seznam_pk = (pk_oct, pk_dcp)
for i in range(2):
    for j in range(2):
        K2 += seznam_x[i]*seznam_x[j]*((((seznam_Tk[i]/seznam_pk[i])**(1/3))+((seznam_Tk[j]/seznam_pk[j])**(1/3)))**3)
K2 *= 0.125

Tk_ = (K1**2)/K2
pk_ = (K1/K2)**2

Vm_ = N(rovnice_.subs([(Tk, Tk_), (pk, pk_), (Vk, Vk_)]))
if im(Vm_) < 10**21:
    Vm_ = re(Vm_)
ro_r_j = 1/(Vm_*(10**3))

# Výsledné hustoty látek a směsí
print('ro_oct_rk = ' + str(ro_oct_rk) + ' mol/l\n')
print('ro_dcp_rk = ' + str(ro_dcp_rk) + ' mol/l\n')
print('ro_oct_r = ' + str(ro_oct_r) + ' mol/l\n')
print('ro_dcp_r = ' + str(ro_dcp_r) + ' mol/l\n')
print('\n\n')
print('ro_smes_rk_a = ' + str(ro_rk_a) + ' mol/l\n')
print('ro_smes_rk_k = ' + str(ro_rk_k) + ' mol/l\n')
print('ro_smes_rk_j = ' + str(ro_rk_j) + ' mol/l\n')
print('\n\n')
print('ro_smes_r_a = ' + str(ro_r_a) + ' mol/l\n')
print('ro_smes_r_k = ' + str(ro_r_k) + ' mol/l\n')
print('ro_smes_r_j = ' + str(ro_r_j) + ' mol/l\n')
