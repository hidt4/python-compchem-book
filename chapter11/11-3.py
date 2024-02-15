"""三フッ化ホウ素とアンモニアの相互作用
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import psi4


def optfreq(mol: psi4.core.Molecule,
            theory: str) -> list[float, psi4.core.Wavefunction]:
    """
    分子の構造最適化と振動数計算を行う

    Args:
        mol (psi4.core.Molecule): 計算対象の分子
        theory: 計算レベル

    Returns:
        float: エネルギー
        psi4.core.Wavefunction: 計算のWavefunction

    """
    _, wfn = psi4.optimize(theory,
                           molecule=mol,
                           return_wfn=True)
    energy, wfn = psi4.frequency(theory,
                                 molecule=mol,
                                 ref_gradient=wfn.gradient(),
                                 return_wfn=True)

    return [energy, wfn]


def get_freqs(wavefunction: psi4.core.Wavefunction) -> np.ndarray:
    """
    Wavefunctionオブジェクトから振動数を取り出す
    Args:
        wavefunction: 対象のWavefunctionオブジェクト

    Returns:
        振動数のarrray

    """
    return wavefunction.frequencies().to_array()


psi4.set_num_threads(14)
psi4.set_memory('28GB')
au2kcal = psi4.constants.hartree2kcalmol
theory = 'mp2/cc-pvtz'

bf3 = psi4.geometry('''
0 1
 B                 -0.46887968    0.02074689    0.00000000
 F                 -1.19887968   -1.24365020    0.00000000
 F                  0.99112032    0.02074689    0.00000000
 F                 -1.19887968    1.28514398    0.00000000
 ''')

nh3 = psi4.geometry('''
0 1
 N                  0.10373444   -0.31950207    0.00000000
 H                  0.43705633   -1.26231516    0.00000000
 H                  0.43707354    0.15189811    0.81649673
 H                  0.43707354    0.15189811   -0.81649673
 ''')

bf3_nh3 = psi4.geometry('''
0 1
 B                 -0.46887968    0.02074689    0.00000000
 F                 -1.15713027   -1.17133812   -0.48666667
 F                  0.90762153    0.02074689   -0.48666666
 F                 -1.15713027    1.21283190   -0.48666667
--
0 1
 N                 -0.46887968    0.02074689    1.58000000
 H                  0.00252447    0.83724368    1.91333334
 H                 -1.41168872    0.02074647    1.91333333
 H                  0.00252520   -0.79574948    1.91333334
''')


psi4.set_output_file('bf3.log')
e_bf3, wfn_bf3 = optfreq(bf3, theory=theory)
print(get_freqs(wfn_bf3))

psi4.set_output_file('nh3.log')
e_nh3, wfn_nh3 = optfreq(nh3, theory=theory)
print(get_freqs(wfn_nh3))

psi4.set_output_file('bf3-nh3.log')
_, wfn_complex = optfreq(bf3_nh3, theory=theory)
print(get_freqs(wfn_complex))
e_cp = psi4.energy(theory, bsse_type='cp', molecule=bf3_nh3)

bf3_complex = ''
nh3_complex = ''

for i, line in enumerate(bf3_nh3.save_string_xyz().split('\n')):
    if i < 5:
        bf3_complex += line + '\n'
    elif i == 5:
        nh3_complex += '0 1\n'
        nh3_complex += line + '\n'
    else:
        nh3_complex += line + '\n'

bf3_complex = psi4.geometry(bf3_complex)
nh3_complex = psi4.geometry(nh3_complex)

psi4.set_output_file('bf3_complex.log')
e_bf3_complex = psi4.energy(theory, molecule=bf3_complex)
psi4.set_output_file('nh3_complex.log')
e_nh3_complex = psi4.energy(theory, molecule=nh3_complex)

energy_strain = au2kcal * (e_bf3_complex + e_nh3_complex - e_bf3 - e_nh3)
energy_int = au2kcal * e_cp
print(f'strain energy: {energy_strain: .2f} kcal/mol')
print(f'interaction energy: {energy_int: .2f} kcal/mol')
print(f'stabilized energy: {energy_strain + energy_int: .2f} kcal/mol')
