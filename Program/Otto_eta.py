import numpy as np
from scipy import optimize
import Otto_harm_osc as oho
import Otto_spin_pol as os2

# Funkcje liczace wartosci temperatur T1, T2, T3, T4
def fT1(kappa, alphac, alphah, tauc, tauh, Tc, Th):
    return (Tc * np.exp(alphah * tauh) * \
        (np.exp(alphac * tauc) - 1) + \
        kappa * Th * (np.exp(alphah * tauh) - 1)) / \
        (np.exp(alphah * tauh + alphac * tauc) - 1)

def fT2(kappa, alphac, alphah, tauc, tauh, Tc, Th):
    return (Tc * np.exp(alphah * tauh) * \
        (np.exp(alphac * tauc) - 1) + \
        kappa * Th * (np.exp(alphah * tauh) - 1)) / \
        (kappa * (np.exp(alphah * tauh + alphac * tauc) - 1))

def fT3(kappa, alphac, alphah, tauc, tauh, Tc, Th):
    return (Tc * (np.exp(alphac * tauc) - 1) + \
        kappa * Th * np.exp(alphac * tauc) * \
        (np.exp(alphah * tauh) - 1)) / \
        (kappa * (np.exp(alphah * tauh + alphac * tauc) - 1))
            
def fT4(kappa, alphac, alphah, tauc, tauh, Tc, Th):
    return (Tc * (np.exp(alphac * tauc) - 1) + \
        kappa * Th * np.exp(alphac * tauc) * \
        (np.exp(alphah * tauh) - 1)) / \
        (np.exp(alphah * tauh + alphac * tauc) - 1)

# Maksymalizacja P zostala zaimplementowana jako minimalizacja -P
#def OttoHOFuncToMin(x, params):
def OttoHOFuncToMin(kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta):
    # Wartosci temperatur i mocy dla danych parametrow
    T1 = fT1(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T2 = fT2(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T3 = fT3(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T4 = fT4(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    P = oho.P(omega2, kappa, T1, T2, T3, T4, alphac, alphah, \
        tauc, tauh, zeta)
    return -P

#def OttoS2FuncToMin(x, params):
def OttoS2FuncToMin(kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta):
    # Wartosci temperatur i mocy dla danych parametrow
    T1 = fT1(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T2 = fT2(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T3 = fT3(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    T4 = fT4(kappa, alphac, alphah, tauc, tauh, Tc, Th)
    P = os2.P(omega2, kappa, T1, T2, T3, T4, alphac, alphah, \
        tauc, tauh, zeta)
    return -P

# Obliczenia dla oscylatora harmonicznego
def OttoHarmonicOscillator(alphac, alphah, Th, zeta, grid_size):
    # Wartosci, dla ktorych obliczono sprawnosc
    Tc_grid = np.linspace(0.01, 1, grid_size)
    # Macierz zapisana do pliku
    eta = np.empty((grid_size, 2))
    
    # Obliczenia dla omega2/Tc = 0.1
    for i in range(grid_size):
        omega2 = 0.1 * Tc_grid[i]
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [0.5, 1, 1]
        # Optymalizacja
        # x[0] - kappa
        # x[1] - tauc
        # x[2] - tauh
        fun = lambda x: OttoHOFuncToMin(x[0], x[1], x[2], omega2, \
            alphac, alphah, Tc_grid[i], Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, 1), (0.01, None), (0.01, None)))
        eta[i, 0] = Tc_grid[i]       # Tc
        eta[i, 1] = 1 - result.x[0]  # eta
    # Zapis do pliku
    np.savetxt('Otto_harm_osc_PT_01_data.dat', eta, fmt = '%.6f', \
        delimiter=' ', newline='\n')

    # Obliczenia dla omega2/Tc = 10
    for i in range(grid_size):
        omega2 = 10 * Tc_grid[i]
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [0.5, 2, 2]
        # Optymalizacja
        # x[0] - kappa
        # x[1] - tauc
        # x[2] - tauh
        fun = lambda x: OttoHOFuncToMin(x[0], x[1], x[2], omega2, \
            alphac, alphah, Tc_grid[i], Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, 1), (0.01, None), (0.01, None)))
        eta[i, 0] = Tc_grid[i]       # Tc
        eta[i, 1] = 1 - result.x[0]  # eta
    # Zapis do pliku
    np.savetxt('Otto_harm_osc_PT_10_data.dat', eta, fmt = '%.6f', \
        delimiter=' ', newline='\n')

# Obliczenia dla ukladu o spinie 1/2
def OttoSpin2(alphac, alphah, Th, zeta, grid_size):
    # Wartosci, dla ktorych obliczono sprawnosc
    Tc_grid = np.linspace(0.01, 1, grid_size)
    # Macierz zapisana do pliku
    eta = np.empty((grid_size, 2))

    # Obliczenia dla omega2/Tc = 0.1
    for i in range(grid_size):
        omega2 = 0.1 * Tc_grid[i]
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [0.5, 1, 1]
        # Optymalizacja
        # x[0] - kappa
        # x[1] - tauc
        # x[2] - tauh
        fun = lambda x: OttoS2FuncToMin(x[0], x[1], x[2], omega2, \
            alphac, alphah, Tc_grid[i], Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, 1), (0.01, None), (0.01, None)))
        eta[i, 0] = Tc_grid[i]       # Tc
        eta[i, 1] = 1 - result.x[0]  # eta
    # Zapis do pliku
    np.savetxt('Otto_spin_pol_PT_01_data.dat', eta, fmt = '%.6f', \
        delimiter=' ', newline='\n')

    # Obliczenia dla omega2/Tc = 10
    for i in range(grid_size):
        omega2 = 10 * Tc_grid[i]
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [0.5, 1, 1]
        # Optymalizacja
        # x[0] - kappa
        # x[1] - tauc
        # x[2] - tauh
        fun = lambda x: OttoS2FuncToMin(x[0], x[1], x[2], omega2, \
            alphac, alphah, Tc_grid[i], Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, 1), (0.01, None), (0.01, None)))
        eta[i, 0] = Tc_grid[i]       # Tc
        eta[i, 1] = 1 - result.x[0]  # eta
    # Zapis do pliku
    np.savetxt('Otto_spin_pol_PT_10_data.dat', eta, fmt = '%.6f', \
        delimiter=' ', newline='\n')

def OttoHarmonicOscillatorkappa(alphac, alphah, Th, zeta, grid_size):
    # Wartosci, dla ktorych obliczono sprawnosc
    kappa_grid = np.linspace(1, 0.01, grid_size)
    # Macierz zapisana do pliku
    P = np.empty((grid_size, 2))
    
    # Obliczenia dla omega2/Tc = 0.1
    for i in range(grid_size):
        #params = (omega2, alphac, alphah, Tc_grid[i], Th, zeta)
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [1, 1]
        Tc = Th / 3
        omega2 = 0.1 * Tc
        # Optymalizacja
        #fun = lambda x: OttoHOFuncToMin(x, params)
        # kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta
        fun = lambda x: OttoHOFuncToMin(kappa_grid[i], x[0], x[1], omega2, \
            alphac, alphah, Tc, Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, None), (0.01, None)))
        P[i, 0] = 1 - kappa_grid[i]     # eta
        P[i, 1] = - result.fun if result.fun<0 else 0 # P
    # Zapis do pliku
    np.savetxt('Otto_harm_osc_Peta_01_data.dat', P, fmt = '%.6f', \
        delimiter=' ', newline='\n')

    # Obliczenia dla omega2/Tc = 10
    for i in range(grid_size):
        #params = (omega2, alphac, alphah, Tc_grid[i], Th, zeta)
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [1, 1]
        Tc = Th / 3
        omega2 = 10 * Tc
        # Optymalizacja
        #fun = lambda x: OttoHOFuncToMin(x, params)
        # kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta
        fun = lambda x: OttoHOFuncToMin(kappa_grid[i], x[0], x[1], omega2, \
            alphac, alphah, Tc, Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, None), (0.01, None)))
        P[i, 0] = 1 - kappa_grid[i]     # eta
        P[i, 1] = - result.fun if result.fun<0 else 0 # P
    # Zapis do pliku
    np.savetxt('Otto_harm_osc_Peta_10_data.dat', P, fmt = '%.6f', \
        delimiter=' ', newline='\n')

def OttoSpin2kappa(alphac, alphah, Th, zeta, grid_size):
    # Wartosci, dla ktorych obliczono sprawnosc
    kappa_grid = np.linspace(1, 0.01, grid_size)
    # Macierz zapisana do pliku
    P = np.empty((grid_size, 2))
    
    # Obliczenia dla omega2/Tc = 0.1
    for i in range(grid_size):
        #params = (omega2, alphac, alphah, Tc_grid[i], Th, zeta)
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [1, 1]
        Tc = Th / 3
        omega2 = 0.1 * Tc
        # Optymalizacja
        #fun = lambda x: OttoHOFuncToMin(x, params)
        # kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta
        fun = lambda x: OttoS2FuncToMin(kappa_grid[i], x[0], x[1], omega2, \
            alphac, alphah, Tc, Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, None), (0.01, None)))
        P[i, 0] = 1 - kappa_grid[i]  # eta
        P[i, 1] = - result.fun if result.fun<0 else 0 # P
    # Zapis do pliku
    np.savetxt('Otto_spin_pol_Peta_01_data.dat', P, fmt = '%.6f', \
        delimiter=' ', newline='\n')

    # Obliczenia dla omega2/Tc = 10
    for i in range(grid_size):
        #params = (omega2, alphac, alphah, Tc_grid[i], Th, zeta)
        # Wartosci poczatkowe podczas optymalizacji
        x_init = [1, 1]
        Tc = Th / 3
        omega2 = 10 * Tc
        # Optymalizacja
        #fun = lambda x: OttoHOFuncToMin(x, params)
        # kappa, tauc, tauh, omega2, alphac, alphah, Tc, Th, zeta
        fun = lambda x: OttoS2FuncToMin(kappa_grid[i], x[0], x[1], omega2, \
            alphac, alphah, Tc, Th, zeta)
        result = optimize.minimize(fun, x_init, \
            bounds=((0.01, None), (0.01, None)))
        P[i, 0] = 1 - kappa_grid[i]  # eta
        P[i, 1] = - result.fun if result.fun<0 else 0  # P
    # Zapis do pliku
    np.savetxt('Otto_spin_pol_Peta_10_data.dat', P, fmt = '%.6f', \
        delimiter=' ', newline='\n')

def main():
    # Zadanie podstawowych parametrow
    grid_size = 201  # Ilosc badanych punktow
    alphac = 1       # Wspolczynnik przenikalnosci. Wartosci nie
    alphah = 1       # zmieniaja parametru kappa, a jedynie tau i P
    Th = 1           # Temperatura Tc jest wyrazona w wielokrotnosci Th
    zeta = 1         # Parametr zmienia jedynie wartosc P
    
    # Uruchomienie obu optymalizacji
    OttoHarmonicOscillator(alphac, alphah, Th, zeta, grid_size)
    OttoSpin2(alphac, alphah, Th, zeta, grid_size)
    OttoHarmonicOscillatorkappa(alphac, alphah, Th, zeta, grid_size)
    OttoSpin2kappa(alphac, alphah, Th, zeta, grid_size)

if __name__ == "__main__":
    main()

