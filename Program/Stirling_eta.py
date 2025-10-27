import numpy as np
import Stirling_harm_osc as sho
import Stirling_spin_pol as ss2

def StirlingHarmonicOscillator(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    momega1, mkappa = np.meshgrid(omega1, kappa)
    eta_grid = sho.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    np.savetxt('Stirling_harm_osc_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

def StirlingSpin2(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size):
    omega1 = np.linspace(omega1_min, omega1_max, grid_size)
    kappa = np.linspace(kappa_min, kappa_max, grid_size)
    momega1, mkappa = np.meshgrid(omega1, kappa)
    eta_grid = ss2.eta_grid(momega1, mkappa, Tc, Th, grid_size)
    eta = np.empty((grid_size**2,3))
    for i in range(grid_size):
        eta[grid_size * i:grid_size * (i+1), 0] = omega1[i]
        for j in range(grid_size):
            eta[grid_size * i + j, 1] = kappa[j]
            eta[grid_size * i + j, 2] = eta_grid[j,i]
    np.savetxt('Stirling_spin_pol_data.dat', eta, fmt='%.6f', \
        delimiter=' ', newline='\n')

def main():
    grid_size = 51
    omega1_min = 0.1
    omega1_max = 5
    kappa_min = 1.1
    kappa_max = 5
    Tc = 1
    Th = 3

    StirlingHarmonicOscillator(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)
    StirlingSpin2(omega1_min, omega1_max, \
        kappa_min, kappa_max, Tc, Th, grid_size)

if __name__ == "__main__":
    main()

