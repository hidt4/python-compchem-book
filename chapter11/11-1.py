"""水二量体の相互作用エネルギー
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


psi4.set_num_threads(10)
psi4.set_memory('20GB')
au2kcal = psi4.constants.hartree2kcalmol

# 水二量体
water_dimer = psi4.geometry('''
0 1
 O                  1.41038900   -0.12125100   -0.03877600
 H                  1.57876400    0.65723900   -0.59380500
 H                  1.58039200    0.18907100    0.86523600
--
0 1
 H                 -0.54699100   -0.10878800   -0.03565200
 O                 -1.49308500    0.11629300    0.03640700
 H                 -1.95060200   -0.69785800   -0.21683100
 ''')

psi4.set_output_file('water_dimer.log')
e_dimer, wfn_dimer = optfreq(water_dimer, theory='mp2/cc-pvtz')
print(get_freqs(wfn_dimer))

# 水単量体
water_monomer = psi4.geometry('''
0 1
 O                  1.41038900   -0.12125100   -0.03877600
 H                  1.57876400    0.65723900   -0.59380500
 H                  1.58039200    0.18907100    0.86523600
''')

psi4.set_output_file('water_monomer.log')
e_monomer, wfn_monomer = optfreq(water_monomer, theory='mp2/cc-pvtz')
print(get_freqs(wfn_monomer))

# 超分子法による相互作用エネルギーの計算
e_int = e_dimer - 2 * e_monomer
print(f'interaction energy: {au2kcal * e_int: .2f} kcal/mol')


# CP法による相互作用エネルギーの計算
psi4.set_output_file('water_dimer_cp.log')
e_cp = psi4.energy('mp2/cc-pvtz',
                   bsse_type='cp',
                   molecule=water_dimer)
print(f'CP corrected interaction energy: {au2kcal * e_cp: .2f} kcal/mol')

# BSSEの算出
print(f'BSSE at MP2/cc-pVTZ level: {au2kcal * (e_cp - e_int): .2f} kcal/mol')