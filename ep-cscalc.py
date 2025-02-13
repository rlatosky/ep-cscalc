import math
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description='Calculate EP cross section using different known variables.')
parser.add_argument("-e", "--energy", type=float, help="Energy of the particle in units of MeV.")
parser.add_argument("-t", "--theta", type=float, help="Cross section scattering angle in units of radians.")
parser.add_argument("-ep", "--eprime", type=float, help="Energy prime of the particle interaction in units of MeV.")
parser.add_argument("-q", "--qsquared", type=float, help="Q-squared parameter in particle interaction.")
parser.add_argument("-x", "--x", type=float, help="x input parameter for form factors.")

args = parser.parse_args()

def calculate(energy=args.energy, theta=args.theta, energy_prime=args.eprime, x=args.x, Q_squared=args.qsquared):
    alpha = 1/137
    proton_mass = 938.272089 # Units - MeV/c2

    if energy_prime is not None and theta is not None:
        print("Using energy, energy_prime, and theta with the following parameters:")
        print("Energy (MeV):", energy)
        print("Energy Prime (MeV):", energy_prime)
        print("Theta (radians):", theta)

        Q_squared = 2*energy*energy_prime*(1-math.cos(theta))
    elif Q_squared is not None and x is not None:
        print("Using x, Q-squared, and energy with the following parameters:")
        print("X:", x)
        print("Q-squared:", Q_squared)
        print("Energy (MeV):", energy)

        energy_prime = (energy - (Q_squared**2)/(2*proton_mass*x))
        theta = math.acos((-Q_squared)/(2*energy*energy_prime)+1)
   # elif x is not None and energy_prime is not None:
   #     print("Using energy, energy_prime, and x with the following parameters:")
   #     print("Energy (MeV):", energy)
   #     print("Energy Prime (MeV):", energy_prime)
   #     print("X:", x)

   #     Q_squared = 2*energy*energy_prime*(1-math.cos(theta))
   #     theta = math.acos(((-x*proton_mass)*(energy-energy_prime)/(energy*energy_prime))+1)
    elif Q_squared is not None and theta is not None:
        print("Using energy, theta, and Q-squared with the following parameters:")
        print("Energy (MeV):", energy)
        print("Theta (radians):", theta)
        print("Q-squared:", Q_squared)

        energy_prime = (Q_squared)/((2*energy)*(1-math.cos(theta)))

        #Additional Equations
        #x = Q_squared/(2*proton_mass*(energy - energy_prime))
        #x = ((2*energy*energy_prime)*(1-math.cos(theta)))/((2*proton_mass)*(energy-energy_prime))

    mott_cs = (alpha**2)/(4*(energy**2)*(math.sin(theta/2))**4)*(energy_prime/energy)
    print(f"===\nMott Cross-section:", mott_cs)
    tau = Q_squared/(4*(proton_mass)**2)
    print("Tau:", tau)
    G_E_p, G_M_p = calculate_G(Q_squared)

    return (mott_cs)*( (((G_E_p**2)+(tau)*(G_M_p**2))/(1+tau))*((math.cos(theta/2))**2) + 2*tau*(G_M_p**2)*((math.sin(theta/2))**2) )

def calculate_G(Q_squared):
    # Taken from Supplemental materials, available at: https://doi.org/10.1016/j.physletb.2017.11.023
    a_G_E_p = np.array([0.239163298067,
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
                        -0.981529071103
                    ]) # a fit parameter
    a_G_M_p = np.array([0.264142994136,
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
                        1.280814375890
                    ])

    G_E_p = np.zeros(shape=(12, 1)) # values for G - make sure to sum over
    G_M_p = np.zeros(shape=(12, 1)) # values for G - make sure to sum over

    pion_mass = 139.57 # MeV
    t_0 = -700 # MeV2
    t_cut = 4*(pion_mass**2)
    z = (np.sqrt(t_cut+Q_squared)-np.sqrt(t_cut-t_0))/(np.sqrt(t_cut+Q_squared)+np.sqrt(t_cut-t_0))

    for x in range(12):
        G_E_p[x] = [a_G_E_p[x]*(z**x)]
        G_M_p[x] = [a_G_M_p[x]*(z**x)]

    sum_G_E_p = np.sum(G_E_p)
    sum_G_M_p = np.sum(G_M_p)
    print("Form factor G_E (proton):", sum_G_E_p)
    print("Form factor G_M (proton):", sum_G_M_p)
    return sum_G_E_p, sum_G_M_p

if __name__=="__main__":
    if len(sys.argv) > 1:
        print(f"===\nElectron-Proton scattering cross section calculation (lab frame):", calculate())
    else:
        print("Utilize -e for energy (in MeV)")
        print("Utilize -t for scattering angle (in radians)")
        print("Utilize -ep energy prime of particle interaction (in MeV)")
        print("Utilize -q for the Q-squared parameter")
        print("Utilize -x for form factor")
        print(f"===\nExample: python .py -e [value] -ep [value] -t [value]")
