import numpy as np
import Stirling_harm_osc as sho
import Stirling_spin_pol as ss2
import Stirling_harm_osc_bezreg as shob
import Stirling_spin_pol_bezreg as ss2b

# Obliczenia dla kolejnych warunkow silnika Stirlinga
# Oscylator harmoniczny z regeneratorem
def StirlingHarmonicOscillator(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    # Wartosci omega1 i kappa, dla ktorych obliczono sprawnosc
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    # Stworzenie dwuwymiarowej sieci dla ulatwienia obliczen
    momega1, mkappa = np.meshgrid(omega1, kappa)
    # Wartosc eta dla kazdej komorki sieci
    eta_grid = sho.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    # Przygotowanie tablicy do zapisu
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    # Zapis do pliku
    np.savetxt('Stirling_harm_osc_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

# Uklad o spinie 1/2 z regeneratorem
def StirlingSpin2(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    # Wartosci omega1 i kappa, dla ktorych obliczono sprawnosc
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    # Stworzenie dwuwymiarowej sieci dla ulatwienia obliczen
    momega1, mkappa = np.meshgrid(omega1, kappa)
    # Wartosc eta dla kazdej komorki sieci
    eta_grid = ss2.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    # Przygotowanie tablicy do zapisu
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    # Zapis do pliku
    np.savetxt('Stirling_spin_pol_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

# Oscylator harmoniczny bez regeneratora
def StirlingHarmonicOscillatorBezReg(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    # Wartosci omega1 i kappa, dla ktorych obliczono sprawnosc
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    # Stworzenie dwuwymiarowej sieci dla ulatwienia obliczen
    momega1, mkappa = np.meshgrid(omega1, kappa)
    # Przygotowanie tablicy do zapisu
    eta_grid = shob.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    # Przygotowanie tablicy do zapisu
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    # Zapis do pliku
    np.savetxt('Stirling_harm_osc_bezreg_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

# Uklad o spinie 1/2 bez regeneratora
def StirlingSpin2BezReg(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    # Wartosci omega1 i kappa, dla ktorych obliczono sprawnosc
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    # Stworzenie dwuwymiarowej sieci dla ulatwienia obliczen
    momega1, mkappa = np.meshgrid(omega1, kappa)
    # Przygotowanie tablicy do zapisu
    eta_grid = ss2b.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    # Zapis do pliku
    np.savetxt('Stirling_spin_pol_bezreg_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

def main():
    # Zadanie podstawowych parametrow
    grid_size = 51    # Ilosc badanych punktow
    omega1_min = 0.1  # Zakres wartosci omega1
    omega1_max = 5
    kappa_min = 1.1   # zakres wartosci kappa
    kappa_max = 5
    Tc = 1            # Temperatura Th wyrazona w wielokrotnosci Tc
    Th = 3
    
    # Uruchomienie obliczen
    StirlingHarmonicOscillator(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)
    StirlingSpin2(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)
    StirlingHarmonicOscillatorBezReg(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)
    StirlingSpin2BezReg(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)

if __name__ == "__main__":
    main()

