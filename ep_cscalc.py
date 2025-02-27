import math
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(
    description="Calculate EP cross section using different known variables."
)
parser.add_argument(
    "-e", "--energy", type=float, help="Energy of the particle in units of MeV."
)
parser.add_argument(
    "-t",
    "--theta",
    type=float,
    help="Cross section scattering angle in units of radians.",
)
parser.add_argument(
    "-ep",
    "--eprime",
    type=float,
    help="Energy prime of the particle interaction in units of MeV.",
)
parser.add_argument(
    "-q", "--qsquared", type=float, help="Q-squared parameter in particle interaction."
)
parser.add_argument("-x", "--x", type=float, help="x input parameter for form factors.")

args = parser.parse_args()

# GOAL: event rate of protons = luminosity * cross section

def calculate(
    energy=args.energy,
    theta=args.theta,
    energy_prime=args.eprime,
    x=args.x,
    Q_squared=args.qsquared,
):
    alpha = 1 / 137
    proton_mass = .938272089  # Units - GeV/c2

    if energy_prime is not None and theta is not None:
        print("Using energy, energy_prime, and theta with the following parameters:")
        print("Energy (MeV):", energy)
        print("Energy Prime (MeV):", energy_prime)
        print("Theta (degrees):", theta)

        Q_squared = 2 * energy * energy_prime * (1 - math.cos(math.radians(theta)))
    elif Q_squared is not None and x is not None:
        print("Using x, Q-squared, and energy with the following parameters:")
        print("X:", x)
        print("Q-squared:", Q_squared)
        print("Energy (MeV):", energy)

        energy_prime = energy - (Q_squared**2) / (2 * proton_mass * x)
        theta = math.acos(((-Q_squared) / (2 * energy * energy_prime)) + 1)
    elif Q_squared is not None and theta is not None:
        print("Using energy, theta, and Q-squared with the following parameters:")
        print("Energy (MeV):", energy)
        print("Theta (degrees):", theta)
        print("Q-squared:", Q_squared)

        energy_prime = (Q_squared) / ((2 * energy) * (1 - math.cos(math.radians(theta))))

        # Additional Equations
        # x = Q_squared/(2*proton_mass*(energy - energy_prime))
        # x = ((2*energy*energy_prime)*(1-math.cos(theta)))/((2*proton_mass)*(energy-energy_prime))
    mott_cs = (
        (alpha**2)
        * (energy_prime / energy)
        #* ((math.cos(math.radians(theta)/2))**2)
        / (4 * (energy**2) * ((math.sin(math.radians(theta/2))) ** 4))
    )
    print("===","Mott Cross-section:", mott_cs)
    tau = Q_squared / (4 * (proton_mass**2))
    print("Tau:", tau)

    photon_pol = (1 + 2*(1+tau)*(math.tan(math.radians(theta/2)))**2)**(-1)
    print("photon_pol:", photon_pol)

    G_E_p, G_M_p = calculate_G_proton(Q_squared)
    G_E_n, G_M_n = calculate_G_neutron(Q_squared)

    G_D = 1 / (1+(Q_squared**2/0.71))**2

    reduced_cs = photon_pol*(G_E_p**2)+tau*(G_M_p**2)
    print("reduced Cross-section:", reduced_cs)

    return ((mott_cs) * ((reduced_cs))) / (photon_pol*(1+tau))
#    return (mott_cs) * (
#        (((G_E_p**2) + (tau) * (G_M_p**2)) / (1 + tau)) * ((math.cos(math.radians(theta) / 2)) ** 2)
#        + 2 * tau * (G_M_p**2) * ((math.sin(math.radians(theta) / 2)) ** 2)
#    )


def calculate_G_proton(Q_squared):
    # Taken from Supplemental materials, available at: https://doi.org/10.1016/j.physletb.2017.11.023
    a_G_E_p = np.array(
        [
            0.239163298067,
            -1.109858574410,
            1.444380813060,
            0.479569465603,
            -2.286894741870,
            1.126632984980,
            1.250619843540,
            -3.631020471590,
            4.082217023790,
            0.504097346499,
            -5.085120460510,
            3.967742543950,
            -0.981529071103,
        ]
    )  # a fit parameter
    a_G_M_p = np.array(
        [
            0.264142994136,
            -1.095306122120,
            1.218553781780,
            0.661136493537,
            -1.405678925030,
            -1.356418438880,
            1.447029155340,
            4.235669735900,
            -5.334045653410,
            -2.916300520960,
            8.707403067570,
            -5.706999943750,
            1.280814375890,
        ]
    )

    G_E_p = np.zeros(12)  # values for G - make sure to sum over
    G_M_p = np.zeros(12)  # values for G - make sure to sum over

    pion_mass = 0.13957  # GeV
    t_0 = -0.7  # GeV2
    t_cut = 4 * (pion_mass**2)


    z = (np.sqrt(t_cut + Q_squared) - np.sqrt(t_cut - t_0)) / (
        np.sqrt(t_cut + Q_squared) + np.sqrt(t_cut - t_0)
    )

    for x in range(12):
        G_E_p[x] = a_G_E_p[x] * (z**x)
        G_M_p[x] = a_G_M_p[x] * (z**x)

    sum_G_E_p = np.sum(G_E_p)
    sum_G_M_p = np.sum(G_M_p)
    print("Form factor G_E (proton):", sum_G_E_p)
    print("Form factor G_M (proton):", sum_G_M_p)
    return sum_G_E_p, sum_G_M_p


def calculate_G_neutron(Q_squared):

    # Taken from Supplemental materials, available at: https://doi.org/10.1016/j.physletb.2017.11.023
    a_G_E_n = np.array(
        [
            0.048919981379,
            -0.064525053912,
            -0.240825897382,
            0.392108744873,
            0.300445258602,
            -0.661888687179,
            -0.175639769687,
            0.624691724461,
            -0.077684299367,
            -0.236003975259,
            0.090401973470,
        ]
    )  # a fit parameter
    a_G_M_n = np.array(
        [
            0.257758326959,
            -1.079540642058,
            1.182183812195,
            0.711015085833,
            -1.348080936796,
            -1.662444025208,
            2.624354426029,
            1.751234494568,
            -4.922300878888,
            3.197892727312,
            -0.712072389946,
        ]
    )

    G_E_n = np.zeros(12)  # values for G - make sure to sum over
    G_M_n = np.zeros(12)  # values for G - make sure to sum over

    pion_mass = .13957  # GeV/c^2
    t_0 = -0.7  # GeV2
    t_cut = 4 * (pion_mass**2)
    z = (np.sqrt(t_cut + Q_squared) - np.sqrt(t_cut - t_0)) / (
        np.sqrt(t_cut + Q_squared) + np.sqrt(t_cut - t_0)
    )

    for x in range(10):
        G_E_n[x] = a_G_E_n[x] * (z**x)
        G_M_n[x] = a_G_M_n[x] * (z**x)

    sum_G_E_n = np.sum(G_E_n)
    sum_G_M_n = np.sum(G_M_n)
    print("Form factor G_E (neutron):", sum_G_E_n)
    print("Form factor G_M (neutron):", sum_G_M_n)
    return sum_G_E_n, sum_G_M_n



if __name__ == "__main__":
    if len(sys.argv) > 3:
        print(
            "===","E_P scattering cs (lab frame):",
            calculate(),
            "barns/steradians"
        )
    else:
        print("Utilize -e for energy (in GeV)")
        print("Utilize -t for scattering angle (in degrees)")
        print("Utilize -ep energy prime of particle interaction (in GeV)")
        print("Utilize -q for the Q-squared parameter")
        print("Utilize -x for form factor")
        print("Example: python .py -e [value] -ep [value] -t [value]")
