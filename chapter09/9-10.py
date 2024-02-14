"""アセトアルデヒドの結合次数
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


psi4.set_output_file('aldehyde_bond_indices.log')

# アセトアルデヒドの構造定義
aldehyde = psi4.geometry('''
0 1
 C                 -1.12448134   -0.00414938    0.00000000
 O                  0.10578584   -0.00415033   -0.00211348
 H                 -1.83462756   -0.84827544    0.00202710
 C                 -1.83871980    1.36020521    0.00000000
 H                 -1.93336158    1.71324347    1.00563757
 H                 -1.26795122    2.06158523   -0.57200462
 H                 -2.81110303    1.25374715   -0.43363296
 ''')

# 構造最適化・振動数計算
_, aldehyde_wfn = optfreq(aldehyde, 'b3lyp/cc-pvdz')
print(get_freqs(aldehyde_wfn))

# 結合次数の計算
_, aldehyde_wfn = psi4.properties('b3lyp/cc-pvdz',
                                  molecule=aldehyde,
                                  properties=['WIBERG_LOWDIN_INDICES'],
                                  return_wfn=True)
aldehyde_df = pd.DataFrame(aldehyde_wfn.variable('WIBERG LOWDIN INDICES').to_array())
aldehyde_df.index = [aldehyde.symbol(i) for i in range(aldehyde.natom())]
print(aldehyde_df.round(2))
