# MKWS 2026 – Optymalizacja termochemiczna silnika rakietowego N₂O/Ethanol

[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![Cantera](https://img.shields.io/badge/Cantera-3.2%2B-orange)](https://cantera.org)

Analiza równowagi termochemicznej i optymalizacja parametrów pracy silnika rakietowego na mieszankę **N₂O (utleniacz) / C₂H₅OH (paliwo)** przy użyciu biblioteki **Cantera** z mechanizmem **GRI-Mech 3.0**.

Projekt wykonany w ramach kursu *"Metody komputerowe w spalaniu" (MKWS 2026)*.

---

## Spis treści

- [Opis projektu](#opis-projektu)
- [Wyniki](#wyniki)
- [Wymagania](#wymagania)
- [Uruchomienie](#uruchomienie)
- [Struktura plików](#struktura-plików)
- [Figury](#figury)
- [Referencje](#referencje)

---

## Opis projektu

Celem projektu jest znalezienie optymalnego stosunku masowego utleniacza do paliwa (O/F) dla hipotetycznego silnika rakietowego na ciekłe N₂O i etanol, przy założeniach:

| Parametr             | Wartość  |
|----------------------|----------|
| Ciśnienie w komorze  | 45 bar   |
| Ciśnienie otoczenia  | 1 atm    |
| Temperatura wlotowa  | 300 K    |

### Metodologia

1. **Generowanie mechanizmu** – Redukcja GRI-Mech 3.0 do 11 gatunków + dodanie etanolu (C₂H₅OH) ze współczynnikami NASA7 z bazy Burcata.
2. **Obliczenia równowagowe** – Dla 281 punktów (O/F = 1.0–15.0, krok 0.05) wyznaczono skład równowagowy metodą HP (stała entalpia, stałe ciśnienie) w Cantera.
3. **Parametry silnika** – Obliczono charakterystyczną prędkość `c*`, współczynnik ciągu `Cf`, ciśnienie na wyjściu z dyszy, liczbę Macha na wyjściu `Me`, stopień rozszerzenia dyszy `ε` oraz impuls właściwy `Isp` (funkcja celu) według teorii dyszy Vandenkerckhove'a.
4. **Optymalizacja** – Maksymalizacja `Isp` względem O/F.

---

## Wyniki

### Optymalny punkt pracy

| Parametr                      | Wartość     | Jednostka   |
|-------------------------------|-------------|-------------|
| **Optymalny O/F**             | **4.10**    | –           |
| **Maksymalny Isp**            | **221.03**  | s           |
| Prędkość charakterystyczna c* | 1422.7      | m/s         |
| Temperatura w komorze Tc      | 3194.4      | K           |
| Średnia masa molowa M_mean    | 24.78       | kg/kmol     |
| Wykładnik izentropy γ         | 1.234       | –           |
| Liczba Macha na wyjściu Me    | 3.000       | –           |
| Stopień rozszerzenia ε        | 6.092       | –           |
| Współczynnik ciągu Cf         | 1.524       | –           |

### Skład produktów (przy O/F = 4.10)

| Gatunek | Udział molowy |
|---------|--------------:|
| N₂      |        45.01% |
| H₂O     |        24.49% |
| CO      |        13.70% |
| CO₂     |         7.39% |
| H₂      |         5.84% |
| OH      |         1.61% |
| H       |         1.00% |
| NO      |         0.51% |
| O₂      |         0.27% |
| O       |         0.17% |

---

## Wymagania

- Python 3.10+
- [Cantera](https://cantera.org) 3.2+
- NumPy
- Matplotlib

Instalacja zależności:

```bash
pip install cantera numpy matplotlib
```

---

## Uruchomienie

```bash
# Krok 1 – wygeneruj plik mechanizmu (tylko raz)
python generate_yaml.py

# Krok 2 – uruchom optymalizację (generuje figury i wydruk w konsoli)
python lre_optimization.py

# Krok 3 – eksportuj tabele CSV (opcjonalnie, do raportu)
python export_tables.py
```

Skrypty należy uruchamiać z głównego katalogu projektu.

---

## Struktura plików

```
.
├── generate_yaml.py          # Generator mechanizmu Cantera
├── lre_optimization.py       # Główny skrypt optymalizacyjny
├── export_tables.py          # Eksport wyników do CSV
├── README.md                 # Ten plik
├── Raport_MKWS2026_ML.pdf    # Pełny raport (PDF, jęz. polski)
├── mechanisms/
│   └── n2o_ethanol.yaml      # Mechanizm N₂O/etanol dla Cantery
├── figures/
│   ├── lre_performance.png   # Isp, c*, Tc, M_mean vs O/F
│   ├── lre_composition.png   # Skład produktów vs O/F
│   └── lre_isp_vs_epsilon.png# Isp vs stopień rozszerzenia dyszy
└── tables/
    ├── performance_sweep.csv # Pełny przegląd O/F (281 punktów)
    ├── optimal_point.csv     # Metryki w punkcie optymalnym
    └── composition_optimal.csv# Skład przy optymalnym O/F
```

---

## Figury

![Performance](figures/lre_performance.png)
*Parametry silnika w funkcji stosunku O/F: Isp, c*, Tc, średnia masa molowa.*

![Composition](figures/lre_composition.png)
*Skład molowy produktów spalania w funkcji stosunku O/F.*

![Isp vs epsilon](figures/lre_isp_vs_epsilon.png)
*Zależność impulsu właściwego od stopnia rozszerzenia dyszy (przy optymalnym O/F).*

---

## Referencje

1. **Sutton, G. P., Biblarz, O.** – *Rocket Propulsion Elements*, 9th ed., Wiley, 2017.
2. **Goodwin, D. G., et al.** – *Cantera: An Object-oriented Software Toolkit for Chemical Kinetics, Thermodynamics, and Transport Processes*, https://cantera.org.
3. **Smith, G. P., et al.** – *GRI-Mech 3.0*, http://combustion.berkeley.edu/gri-mech/.
4. **Burcat, A., Ruscic, B.** – *Third Millennium Ideal Gas and Condensed Phase Thermochemical Database for Combustion*.
5. **Vandenkerckhove, J. A.** – *Isentropic Nozzle Theory* (referenced in Sutton & Biblarz).
