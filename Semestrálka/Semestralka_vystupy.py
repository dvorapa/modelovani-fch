## Funkce společné všem metodám
import matplotlib.pyplot as plt

format_souboru = 'Semestralka %s.png'

mix1 = ('A', 'T')
mix2 = ('A', 'C')
mixy = (mix1, mix2)

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
  plt.plot(vysledky[0], vysledky[1], 'bD--', vysledky[0], vysledky[2], 'g^--', ms=5)
  plt.legend(('γ2', 'γ1'))
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('γ')
  jmeno = ('Závislost aktivitních koeficientů na poměru složek směsi\n%s metodou %s'
           % (' + '.join(mix), metoda))
  plt.title(jmeno)
  soubor = 'γ %s %s' % ('+'.join(mix), metoda)
  plt.savefig(format_souboru % soubor)

def graf_logaritmu(metoda, mix, vysledky):
  # Graf logaritmu poměru aktivitních koeficientů
  plt.figure()
  plt.plot(vysledky[0], vysledky[4], 'gD--', ms=5)
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('ln(γ1/γ2)')
  jmeno = ('Závislost logaritmu poměru aktivitních koeficientů na poměru složek směsi\n%s metodou %s'
           % (' + '.join(mix), metoda))
  plt.title(jmeno)
  soubor = 'ln %s %s' % ('+'.join(mix), metoda)
  plt.savefig(format_souboru % soubor)
  

def graf_Q(metoda, vysledky):
  # Graf Gibbsovy energie
  plt.figure(figsize = (9, 3))
  plt.plot(vysledky[0][0], vysledky[0][3], 'b^--', vysledky[1][0], vysledky[1][3], 'rD--', ms=5)
  plt.legend(tuple('směs ' + ' + '.join(mix) for mix in mixy))
  plt.xlabel('x1 [hm. %]')
  plt.ylabel('Q [J/(K⋅mol)]')
  jmeno = 'Závislost Gibbsovy energie na poměru složek směsí metodou ' + metoda
  plt.title(jmeno)
  soubor = 'Gibbs ' + metoda
  plt.savefig(format_souboru % soubor)
