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

metals = ("Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu")
metal_masses = {"Ti" : 47.87,
                "V"  : 50.94,
                "Cr" : 52.00,
                "Mn" : 54.94,
                "Fe" : 55.85,
                "Co" : 58.93,
                "Ni" : 58.69,
                "Cu" : 63.55}

# number of d-electrons per metal
metal_e = dict(zip(metals, range(4, 13)))

# number of unpaired electrons assuming high-spin
high_spin = (0,1,2,3,4,5,4,3,2,1,0)

# number of unpaired electrons assuming low-spin
low_spin = (0,1,2,3,2,1,0,1,2,1,0)

# magmnetic moments per number of unpaired e
mag_moments = {1 : 1.73,
               2 : 2.83,
               3 : 3.87,
               4 : 4.90,
               5 : 5.92}


def calc_assignments(data):
    data['solvent_vol'] = data['solvent_mass'] / solvent_density
    data.insert(4, 'delta_hz',  data['delta_ppm'] * nmr_freq)

    # data frame of calculated magnetic moments based on identity
    mag_moments_df = pd.DataFrame(data['sample'])
    hs_calculations = pd.DataFrame(data['sample'])
    ls_calculations = pd.DataFrame(data['sample'])

    for element in metal_masses:
        MW = metal_masses[element] + ligand_mass * data['ligand_no']
        conc = (data['solid_mass'] / MW) / data['solvent_vol']
        mag_susc = (3 * data['delta_hz']) / (4 * np.pi * (nmr_freq * 1E6) * conc)
        mag_moment = np.sqrt(8 * data['temp'] * mag_susc)
        mag_moments_df[element] = mag_moment

        hs_elec = (metal_e[element] - data['ligand_no']).apply(lambda n: high_spin[n])
        ls_elec = (metal_e[element] - data['ligand_no']).apply(lambda n: low_spin[n])

        hs_calculations[element] = np.abs(mag_moment - hs_elec.apply(lambda n: mag_moments[n] if 0 < n < 6 else np.inf))
        ls_calculations[element] = np.abs(mag_moment - ls_elec.apply(lambda n: mag_moments[n] if 0 < n < 6 else np.inf))

    print("Calculated magnetic moments:")
    print(mag_moments_df.to_markdown())
    print("Differences assuming high spin:")
    print(hs_calculations.to_markdown())
    print("Differences assuming low spin:")
    print(ls_calculations.to_markdown())

    print("Summary:")
    summary = data.drop(['ligand_no', 'delta_ppm', 'solvent_vol'], axis=1).copy()
    summary = summary.drop([0,1,6])
    summary['mag_mom'] = [0, 2.57, 2.97, 1.96]
    print(summary.to_markdown())
    print(summary.to_latex())


    print(mag_moments_df.to_latex())
    print(hs_calculations.to_latex())
    print(ls_calculations.to_latex())

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
