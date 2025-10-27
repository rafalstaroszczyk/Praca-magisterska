import numpy as np

def fU(T, omega):
    return 0.5 * omega * 1/np.tanh(0.5 * omega / T)

def fS(T, omega):
    return 0.5 * omega / T * 1/np.tanh(0.5 * omega / T) - \
        np.log(2 * np.sinh(0.5 * omega / T))

def fQ12(omega1, kappa, Tc, Th):
    return Tc * fS(Tc, kappa * omega1) - Tc * fS(Tc, omega1)

def fQ23(omega1, kappa, Tc, Th):
    return fU(Th, kappa * omega1) - fU(Tc, kappa * omega1)

def fQ34(omega1, kappa, Tc, Th):
    return Th * fS(Th, omega1) - Th * fS(Th, kappa * omega1)

def fQ41(omega1, kappa, Tc, Th):
    return fU(Tc, omega1) - fU(Th, omega1)

def fDeltaQ(Q23, Q41):
    return Q23 + Q41

def fW(Q12, Q23, Q34, Q41):
    return Q12 + Q23 + Q34 + Q41

def fQd(Q34, DeltaQ):
    dDeltaQ = DeltaQ
    dDeltaQ[dDeltaQ < 0] = 0
    return Q34 + dDeltaQ

def feta(W, Qd):
    return W/Qd

def eta_grid(momega1, mkappa, Tc, Th, grid_size):
    Q12 = fQ12(momega1, mkappa, Tc, Th)
    Q23 = fQ23(momega1, mkappa, Tc, Th)
    Q34 = fQ34(momega1, mkappa, Tc, Th)
    Q41 = fQ41(momega1, mkappa, Tc, Th)
    DeltaQ = fDeltaQ(Q23, Q41)
    W = fW(Q12, Q23, Q34, Q41)
    Qd = fQd(Q34, DeltaQ)
    eta = feta(W, Qd)
    return eta

