### Skript pro výpočet hustoty dvojice kapalin kyselina octová - 1,2-dichloropropan
from sympy import *
Vm = Symbol('Vm')
R = 8.314

## Nadefinované podmínky
T = 273.15 + 25 # K
p = 101325 # Pa

# Kritické veličiny pro kys. octovou (Prausnitz, 2001; NIST, 2018)
Tk_oct = 592.7 # K
pk_oct = 57.9 * 10**5 # Pa
Mw_oct = 60.0520 # g/mol

# Kritické veličiny pro 1,2-dichloropropan (Steele, et al., 1997; NIST, 2018)
Tk_dcp = 578 # K
pk_dcp = 46.5 * 10**5 # Pa
Mw_dcp = 112.986 # g/mol

## Redlich-Kwong
a = 0.42748*(R**2)*(Tk_oct**2.5)/pk_oct
b = 0.08664*R*Tk_oct/pk_oct
rovnice_ = R*T/(Vm-b) - a/(sqrt(T)*Vm*(Vm+b)) - p
koreny_oct = solve(rovnice_, Vm)
koreny_oct_re = []
for k in koreny_oct:
    koreny_oct_re.append(re(k))
Vm_oct = min(koreny_oct_re)

a = 0.42748*(R**2)*(Tk_dcp**2.5)/pk_dcp
b = 0.08664*R*Tk_dcp/pk_dcp
rovnice_ = R*T/(Vm-b) - a/(sqrt(T)*Vm*(Vm+b)) - p
koreny_dcp = solve(rovnice_, Vm)
koreny_dcp_re = []
for k in koreny_dcp:
    koreny_dcp_re.append(re(k))
Vm_dcp = min(koreny_dcp_re)

# Amagat
Vm_rk = (0.25*Vm_oct + 0.75*Vm_dcp)
Mw_rk = (0.25*Mw_oct + 0.75*Mw_dcp)
ro_rk_a = Mw_rk/(Vm_rk*(10**6))

## Rackett
Vk_oct = R*Tk_oct/pk_oct
Vm_oct = Vk_oct*((pk_oct*Vk_oct/(R*Tk_oct))**(((1-T)/Tk_oct)**(2/7)))

Vk_dcp = R*Tk_dcp/pk_dcp
Vm_dcp = Vk_dcp*((pk_dcp*Vk_dcp/(R*Tk_dcp))**(((1-T)/Tk_dcp)**(2/7)))

# Amagat
Vm_r = (0.25*Vm_oct + 0.75*Vm_dcp)
Mw_r = (0.25*Mw_oct + 0.75*Mw_dcp)
ro_r_a = Mw_r/(Vm_r*(10**6))

# Výsledné výsledky
print('ro_rk_a = ' + str(re(ro_rk_a)) + ' g/cm^3\n')
print('ro_r_a = ' + str(re(ro_r_a)) + ' g/cm^3\n')
