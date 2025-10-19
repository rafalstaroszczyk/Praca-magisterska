import numpy as np


def fQ12(omega1, kappa, Tc, Th):
    return Tc * np.log(2 * np.cosh(kappa * 0.5 * omega1 / Tc)) - \
            Tc * np.log(2 * np.cosh(0.5 * omega1 / Tc)) + \
            0.5 * omega1 * np.tanh(0.5 * omega1 / Tc) - \
            kappa * 0.5 * omega1 * np.tanh(kappa * 0.5 * omega1 / Tc)


def fQ23(omega1, kappa, Tc, Th):
    return kappa * 0.5 * omega1 * np.tanh(kappa * 0.5 * omega1 / Tc) - \
            kappa * 0.5 * omega1 * np.tanh(kappa * 0.5 * omega1 / Th)


def fQ34(omega1, kappa, Tc, Th):
    return Th * np.log(2 * np.cosh(0.5 * omega1 / Th)) - \
            Th * np.log(2 * np.cosh(kappa * 0.5 * omega1 / Th)) + \
            kappa * 0.5 * omega1 * np.tanh(kappa * 0.5 * omega1 / Th) - \
            0.5 * omega1 * np.tanh(0.5 * omega1 / Th)


def fQ41(omega1, kappa, Tc, Th):
    return 0.5 * omega1 * np.tanh(0.5 * omega1 / Th) - \
            0.5 * omega1 * np.tanh(0.5 * omega1 / Tc)


def fDeltaQ(Q23, Q41):
    return Q23 + Q41


def fW(Q12, Q23, Q34, Q41):
    return Q12 + Q23 + Q34 + Q41


def fQd(Q34, DeltaQ):
    dDeltaQ = DeltaQ
    dDeltaQ[dDeltaQ > 0] = 0
    return Q34 + dDeltaQ


def feta(W, Qd):
    return W/Qd


def main():
    samples = 101
    Tc = 1
    Th = 3
    omega1 = np.linspace(0, 5, samples)
    kappa = np.linspace(1, 5, samples)
    momega1, mkappa = np.meshgrid(omega1, kappa)
    Q12 = fQ12(momega1, mkappa, Tc, Th)
    Q23 = fQ23(momega1, mkappa, Tc, Th)
    Q34 = fQ34(momega1, mkappa, Tc, Th)
    Q41 = fQ41(momega1, mkappa, Tc, Th)
    DeltaQ = fDeltaQ(Q23, Q41)
    W = fW(Q12, Q23, Q34, Q41)
    Qd = fQd(Q34, DeltaQ)
    eta = feta(W, Qd)
    #print(np.sign(Q12))
    #print(np.sign(Q23))
    #print(np.sign(Q34))
    #print(np.sign(Q41))
    #print(np.sign(DeltaQ))
    #print(np.sign(W))
    print(Qd)
    print(eta)

    


if __name__ == "__main__":
    main()

