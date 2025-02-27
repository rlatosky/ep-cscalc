from ep_cscalc import *
import sys
import matplotlib
import matplotlib.pyplot as plt

def graph_photo_pol_cs():
    alpha = 1 / 137
    proton_mass = .938272089  # Units - GeV/c2
    energy = np.array([
        1.148,
        1.148,
        1.882,
        2.235,
        2.235,
        2.235,
        2.235,
        2.235,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
    ])
    Q_squared = np.array([
        0.6200,
        0.8172,
        0.8995,
        0.6182,
        1.1117,
        1.6348,
        2.2466,
        2.7802,
        0.4241,
        0.6633,
        0.9312,
        2.0354,
        2.6205,
        3.1685,
        3.7561,
        4.2330
    ])
    theta = np.array([
        47.97,
        59.99,
        33.95,
        21.97,
        31.95,
        42.97,
        58.97,
        79.97,
        12.47,
        15.97,
        19.46,
        32.97,
        40.97,
        49.97,
        64.97,
        77.97,
    ])
    tau = np.zeros(len(energy))
    photon_pol = np.zeros(len(energy))
    energy_prime = np.zeros(len(energy))
    ep_cs = np.zeros(len(energy))

    for x in range(len(energy)):
        tau[x] = Q_squared[x] / (4 * (proton_mass**2))
        photon_pol[x] = (1 + 2*(1+tau[x])*(math.tan(math.radians(theta[x]/2)))**2)**(-1)
        energy_prime[x] = energy[x] - (Q_squared[x]**2) / (2 * proton_mass * x)
        ep_cs[x] = calculate(energy=energy[x], theta=theta[x], Q_squared=Q_squared[x])
    matplotlib.use("TkAgg")
    plt.plot(photon_pol, ep_cs)
    plt.xlabel("photon_pol")
    plt.ylabel("sigma_ep")
    plt.show()

def graph_G_q_squared():
    qs = np.linspace(0,6.5,num=1000)
    g = np.zeros(1000)
    nucleon_mm = 3.152e-17

    for x in range(1000):
        G_E, G_M = calculate_G_proton(qs[x])
        g[x] = (nucleon_mm * G_E) / G_M
    matplotlib.use("TkAgg")
    plt.plot(qs, g)
    plt.xlabel("Q_squared")
    plt.ylabel("mu*G_E/G_M")
    plt.show()

def graph_E_prime_cs():
    alpha = 1 / 137
    proton_mass = .938272089  # Units - GeV/c2
    energy = np.array([
        1.148,
        1.148,
        1.882,
        2.235,
        2.235,
        2.235,
        2.235,
        2.235,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
        3.114,
    ])
    Q_squared = np.array([
        0.6200,
        0.8172,
        0.8995,
        0.6182,
        1.1117,
        1.6348,
        2.2466,
        2.7802,
        0.4241,
        0.6633,
        0.9312,
        2.0354,
        2.6205,
        3.1685,
        3.7561,
        4.2330
    ])
    theta = np.array([
        47.97,
        59.99,
        33.95,
        21.97,
        31.95,
        42.97,
        58.97,
        79.97,
        12.47,
        15.97,
        19.46,
        32.97,
        40.97,
        49.97,
        64.97,
        77.97,
    ])
    energy_prime = np.zeros(len(energy))
    ep_cs = np.zeros(len(energy))
    for x in range(len(energy)):
        energy_prime[x] = energy[x] - (Q_squared[x]**2) / (2 * proton_mass * x)
        ep_cs[x] = calculate(energy=energy[x], theta=theta[x], Q_squared=Q_squared[x])
    matplotlib.use("TkAgg")
    plt.plot(energy_prime, ep_cs)
    plt.xlabel("E'")
    plt.ylabel("sigma_ep")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        graph_E_prime_cs()
        graph_photo_pol_cs()
        graph_G_q_squared()
