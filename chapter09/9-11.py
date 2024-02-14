"""N-メチルギ酸アミドのWiberg結合次数
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


psi4.set_output_file('amide_bond_indices.log')

# N-メチルギ酸アミドの構造定義
amide = psi4.geometry('''
0 1
 C                 -1.12448134   -0.00414938    0.00000000
 O                  0.10578766   -0.00414938    0.00000000
 H                 -1.83462934   -0.84827638    0.00004800
 C                 -0.82455772    2.25265795   -0.26299920
 H                 -1.33919615    3.19075802   -0.25905562
 H                 -0.06763533    2.25790797    0.49326857
 H                 -0.37150536    2.09840882   -1.22000012
 N                 -1.75383934    1.19806462    0.00000000
 H                 -2.46211626    1.19926398   -0.70593369
 ''')

# 構造最適化・振動数計算
_, amide_wfn = optfreq(amide, 'b3lyp/cc-pvdz')
print(get_freqs(amide_wfn))

# 結合次数の計算
_, amide_wfn = psi4.properties('b3lyp/cc-pvdz',
                               molecule=amide,
                               properties=['WIBERG_LOWDIN_INDICES'],
                               return_wfn=True)
amide_df = pd.DataFrame(amide_wfn.variable('WIBERG LOWDIN INDICES').to_array())
amide_df.index = [amide.symbol(i) for i in range(amide.natom())]
print(amide_df.round(2))
