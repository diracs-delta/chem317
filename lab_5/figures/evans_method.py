#!/usr/bin/python

import sys
import pandas as pd
import numpy as np

# assume ligand to be acetylacetonate (in g/mol)
ligand_mass = 99.11
# assume solvent to be chloroform (in g/mL)
solvent_density = 1.489
# assume NMR freq (in MHz)
nmr_freq = 500

metal_masses = {"Sc" : 44.96,
                "Ti" : 47.87,
                "V"  : 50.94,
                "Cr" : 52.00,
                "Mn" : 54.94,
                "Fe" : 55.85,
                "Co" : 58.93,
                "Ni" : 58.69,
                "Cu" : 63.55,
                "Zn" : 65.38}

# number of d-electrons per metal
metal_e = {"Sc" : 3,
           "Ti" : 4,
           "V"  : 5,
           "Cr" : 6,
           "Mn" : 7,
           "Fe" : 8,
           "Co" : 9,
           "Ni" : 10,
           "Cu" : 11,
           "Zn" : 12}

# magmnetic moments per number of unpaired e
mag_moments = {1 : 1.73,
               2 : 2.83,
               3 : 3.87,
               4 : 4.90,
               5 : 5.92}

def calc_assignments(data):
    data['solvent_vol'] = data['solvent_mass'] / solvent_density
    data['delta_hz'] = data['delta_ppm'] * nmr_freq

    # data frame of calculated magnetic moments based on identity
    calculations = pd.DataFrame(data['sample'])

    for element in metal_masses:
        MW = metal_masses[element] + ligand_mass * data['ligand_no']
        conc = (data['solid_mass'] / MW) / data['solvent_vol']
        mag_susc = (3 * data['delta_hz']) / (4 * np.pi * (nmr_freq * 1E6) * conc)
        mag_moment = np.sqrt(8 * data['temp'] * mag_susc)

        unpaired_e = - np.abs(metal_e[element] - data['ligand_no'] - 5) + 5
        calculations[element] = np.abs(mag_moment - unpaired_e.apply(lambda n: mag_moments[n] if 0 < n < 6 else np.inf))


    print(calculations.to_markdown())

def main():
    col_names = ['sample', 'color', 'solid_mass', 'solvent_mass', 'temp', 'ligand_no', 'delta_ppm']
    data = pd.read_csv(sys.argv[1], skiprows=[0], names=col_names)

    print(data)
    calc_assignments(data)

if __name__ == "__main__":
    if(len(sys.argv) != 2):
        print("USAGE: evans_method.py [csv_file]")
    else:
        main()
