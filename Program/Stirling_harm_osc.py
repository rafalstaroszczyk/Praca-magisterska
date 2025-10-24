import numpy as np

def fQ12(omega1, kappa, Tc, Th):
    return Tc * np.log(2 * np.sinh(0.5 * omega1 / Tc)) - \
        Tc * np.log(2 * np.sinh(kappa * 0.5 * omega1 / Tc)) + \
        kappa * 0.5 * omega1 * 1/np.tanh(kappa * 0.5 * omega1 / Tc) - \
        0.5 * omega1 * 1/np.tanh(0.5 * omega1 / Tc)

def fQ23(omega1, kappa, Tc, Th):
    return kappa * 0.5 * omega1 * 1/np.tanh(kappa * 0.5 * omega1 / Th) - \
        kappa * 0.5 * omega1 * 1/np.tanh(kappa * 0.5 * omega1 / Tc)

def fQ34(omega1, kappa, Tc, Th):
    return Th * np.log(2 * np.sinh(kappa * 0.5 * omega1 / Th)) - \
        Th * np.log(2 * np.sinh(0.5 * omega1 / Th)) + \
        0.5 * omega1 * 1/np.tanh(0.5 * omega1 / Th) - \
        kappa * 0.5 * omega1 * 1/np.tanh(kappa * 0.5 * omega1 / Th)

def fQ41(omega1, kappa, Tc, Th):
    return 0.5 * omega1 * 1/np.tanh(0.5 * omega1 / Tc) - \
        0.5 * omega1 * 1/np.tanh(0.5 * omega1 / Th)

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

