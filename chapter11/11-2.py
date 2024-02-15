"""基底関数系とBSSEの関係を調べる
"""

import numpy as np
import pandas as pd
import psi4


def get_intE(monomer: psi4.core.Molecule,
             dimer: psi4.core.Molecule,
             theory: str) -> list[float, float, float]:
    """
    相互作用エネルギーを見積もる
    Args:
        monomer: 単量体
        dimer: 二量体
        theory: 計算レベル

    Returns:
        float: 超分子法により推定した相互作用エネルギー
        float: CP法により補正した相互作用エネルギー
        float: BSSE

    """
    au2kcal = psi4.constants.hartree2kcalmol
    energy_monomer = psi4.energy(theory, molecule=monomer)
    energy_dimer = psi4.energy(theory, molecule=dimer)
    energy_cp = psi4.energy(theory, bsse_type='cp', molecule=dimer)

    int_supermol = au2kcal * (energy_dimer - 2 * energy_monomer)
    int_cp = au2kcal * energy_cp
    bsse = int_cp - int_supermol

    return [int_supermol, int_cp, bsse]


psi4.set_num_threads(14)
psi4.set_memory('28GB')

hf = pd.DataFrame()
mp2 = pd.DataFrame()
basis_sets = ['sto-3g', 'cc-pvdz', 'cc-pvtz',
              'aug-cc-pvtz', 'aug-cc-pvqz', 'aug-cc-pv5z']

psi4.set_output_file('basis-set-bsse.log')

water_dimer = psi4.geometry('''
0 1
 O    1.412128481591   -0.092045390825   -0.029906791666
 H    1.719638124602    0.616275423632   -0.598773926023
 H    1.722516776196    0.155785589195    0.843086003713
 --
0 1
 H   -0.525178885277   -0.028828516343   -0.009223430301
 O   -1.480679643058    0.092866715086    0.030210284186
 H   -1.829019336855   -0.756267508525   -0.239905293928
''')

water_monomer = psi4.geometry('''
0 1
 O    0.063345196812    0.019471021033    0.000000000000
 H   -0.723797357344    0.564895145188    0.000000000000
 H   -0.281536884383   -0.873914378699    0.000000000000
''')

for basis_set in basis_sets:
    hf_theory = f'hf/{basis_set}'
    mp2_theory = f'mp2/{basis_set}'

    hf[hf_theory] = get_intE(monomer=water_monomer,
                             dimer=water_dimer,
                             theory=hf_theory)
    mp2[mp2_theory] = get_intE(monomer=water_monomer,
                               dimer=water_dimer,
                               theory=mp2_theory)


print('#### HF ####')
print(hf.round(2))

print('#### MP2 ####')
print(mp2.round(2))
