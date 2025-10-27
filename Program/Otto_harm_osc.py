import numpy as np

# Funkcje realizuja wyprowadzone w pracy zaleznosci, 
# kappa = omega1 / omega2
def fU(T, omega):
    return 0.5 * omega * 1/np.tanh(0.5 * omega / T)

def fS(T, omega):
    return 0.5 * omega / T * 1/np.tanh(0.5 * omega / T) - \
        np.log(2 * np.sinh(0.5 * omega / T))

def fQ23(omega2, kappa, T2, T3):
    return fU(T3, omega2) - fU(T2, omega2)

def fQ41(omega2, kappa, T1, T4):
    return fU(T1, kappa * omega2) - fU(T4, kappa * omega2)

def fW(Q23, Q41):
    return Q23 + Q41

def fP(W, tau):
    return W/tau

# Funkcja zwraca moc P dla podanych parametrow
def P(omega2, kappa, T1, T2, T3, T4, alphac, alphah, \
        tauc, tauh, zeta):
    Q23 = fQ23(omega2, kappa, T2, T3)
    Q41 = fQ41(omega2, kappa, T1, T4)
    W = fW(Q23, Q41)
    tau = zeta * (tauh + tauc)
    P = fP(W, tau)
    return P

