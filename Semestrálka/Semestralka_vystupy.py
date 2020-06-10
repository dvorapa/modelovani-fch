## Funkce společné všem metodám
import matplotlib.pyplot as plt

format_souboru = 'Semestralka %s.png'

mix1 = ('A', 'T')
mix2 = ('A', 'C')
mixy = (mix1, mix2)

x1__ = (0.0000, 0.0777, 0.1441, 0.2017, 0.2519, 0.2963, 0.3356, 0.3708, 0.4025, 0.4311, 0.5310, 0.5601, 0.5927, 0.6294, 0.6708, 0.7181, 0.7725, 0.8359, 0.9106, 1.0000)

def tabulka_vysledku(teploty, vysledky):
  # Tabulka výsledků dle https://stackoverflow.com/a/9536084
  # Hlavička
  titulky = ('x1 [hm. %]', 'Texp [K]', 'γ1', 'γ2', 'ln(γ1/γ2)', 'Q [J/(K⋅mol)]')
  format_popisku = '{:>15}' * len(titulky)
  print(format_popisku.format(*titulky))

  # Řádky
  presnost = (4, 2, 4, 4, 6, 2)
  format_radku = ''.join(tuple('{{:>15.{}f}}'.format(p) for p in presnost))
  for T, x1, gam1, gam2, Q, ln in zip(teploty, *vysledky):
    print(format_radku.format(x1, T, float(gam1), float(gam2), float(ln), float(Q)))

def graf_gammy(metoda, mix, vysledky):
  # Graf aktivitních koeficientů
  plt.figure()
  plt.plot(vysledky[0], vysledky[1], 'b--', vysledky[0], vysledky[2], 'g--', ms=5)
  if (mix == ('A', 'T')):
    gam_exp1 = (1.1782, 1.1467, 1.1232, 1.1043, 1.0901, 1.0783, 1.0687, 1.0608, 1.0543, 1.0446, 1.0318, 1.0277, 1.0234, 1.0181, 1.0148, 1.0106, 1.0066, 1.0033, 1.0009, 1.0000)
    gam_exp2 = (1.0000, 1.0008, 1.0029, 1.0057, 1.0080, 1.0120, 1.0154, 1.0186, 1.0218, 1.0248, 1.0367, 1.0407, 1.0450, 1.0502, 1.0564, 1.0633, 1.0713, 1.0799, 1.0906, 1.1004)
  elif (mix == ('A', 'C')):
    gam_exp1 = (1.1816, 1.1492, 1.1256, 1.1071, 1.0924, 1.0805, 1.0708, 1.0684, 1.0628, 1.0550, 1.0418, 1.0320, 1.0283, 1.0239, 1.0194, 1.0149, 1.0105, 1.0066, 1.0032, 1.0000)
    gam_exp2 = (1.0000, 1.0010, 1.0036, 1.0070, 1.0109, 1.0149, 1.0190, 1.0231, 1.0258, 1.0270, 1.0455, 1.0500, 1.0554, 1.0617, 1.0684, 1.0764, 1.0857, 1.0964, 1.1080, 1.122)
  plt.plot(x1__, gam_exp1, 'bD', x1__, gam_exp2, 'g^', ms=5)
  plt.legend(('γ2 calc.', 'γ1 calc.', 'γ2 exp.', 'γ1 exp.'))
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('γ')
  jmeno = ('Závislost aktivitních koeficientů na poměru složek směsi\n%s metodou %s'
           % (' + '.join(mix), metoda))
  plt.title(jmeno)
  soubor = 'γ %s %s' % ('+'.join(mix), metoda)
  plt.savefig(format_souboru % soubor)

def graf_logaritmu(metoda, mix, vysledky):
  # Graf logaritmu poměru aktivitních koeficientů
  plt.figure(figsize = (9, 4))
  plt.plot(vysledky[0], vysledky[4], 'g--', ms=5)
  if (mix == ('A', 'T')):
    ln_exp = (0.163945, 0.136089, 0.113286, 0.093528, 0.078301, 0.063457, 0.0512, 0.040604, 0.031341, 0.019166, -0.004757, -0.012532, -0.020886, -0.031042, -0.040175, -0.050833, -0.062295, -0.073574, -0.085828, -0.095674)
  elif (mix == ('A', 'C')):
    ln_exp = (0.166869, 0.138067, 0.114723, 0.094768, 0.077556, 0.062594, 0.049573, 0.043335, 0.035406, 0.026928, -0.003545, -0.017291, -0.026032, -0.036281, -0.046928, -0.058852, -0.07178, -0.085417, -0.099362, -0.115113)
  plt.plot(x1__, ln_exp, 'gD', ms=5)
  plt.legend(('calc.', 'exp.'))
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('ln(γ1/γ2)')
  jmeno = ('Závislost logaritmu poměru aktivitních koeficientů na poměru složek směsi\n%s metodou %s'
           % (' + '.join(mix), metoda))
  plt.title(jmeno)
  soubor = 'ln %s %s' % ('+'.join(mix), metoda)
  plt.savefig(format_souboru % soubor)
  

def graf_Q(metoda, vysledky):
  # Graf Gibbsovy energie
  plt.figure(figsize = (9, 4))
  plt.plot(vysledky[0][0], vysledky[0][3], 'b--', vysledky[1][0], vysledky[1][3], 'r--', ms=5)
  Q_mix1_exp = (0.00, 6.99, 12.12, 15.90, 18.71, 20.79, 22.32, 23.44, 24.24, 24.80, 25.25, 25.26, 24.84, 24.09, 22.92, 21.13, 18.46, 14.52, 8.69, 0.00)
  Q_mix2_exp = (0.00, 8.30, 14.22, 18.48, 21.55, 23.75, 25.33, 26.43, 27.18, 27.28, 27.85, 27.51, 26.88, 25.91, 24.45, 22.33, 19.31, 15.00, 8.85, 0.00)
  plt.plot(x1__, Q_mix1_exp, 'b^', x1__, Q_mix2_exp, 'rD', ms=5)
  plt.legend(['směs {} + {} calc.'.format(mix[0], mix[1]) for mix in mixy] + ['směs {} + {} exp.'.format(mix[0], mix[1]) for mix in mixy])
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('Q [J/(K⋅mol)]')
  jmeno = 'Závislost Gibbsovy energie na poměru složek směsí metodou ' + metoda
  plt.title(jmeno)
  soubor = 'Gibbs ' + metoda
  plt.savefig(format_souboru % soubor)