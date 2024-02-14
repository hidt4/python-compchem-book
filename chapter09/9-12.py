"""スチレンの結合次数
"""


from __future__ import annotations
import numpy as np
import pandas as pd
import psi4


def optfreq(mol: psi4.core.Molecule,
            theory: str = 'b3lyp-d3bj/6-31g(d)') -> list[float, psi4.core.Wavefunction]:
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


psi4.set_output_file('stylene_bond_indices.log')

# スチレンの構造を定義
stylene = psi4.geometry('''
0 1
 C                 -0.90871371    0.51037344    0.00000000
 C                  0.48644629    0.51037344    0.00000000
 C                  1.18398429    1.71812444    0.00000000
 C                  0.48633029    2.92663344   -0.00119900
 C                 -0.90849471    2.92655544   -0.00167800
 C                 -1.60609571    1.71834944   -0.00068200
 H                 -1.45847271   -0.44194356    0.00045000
 H                  1.03595429   -0.44213956    0.00131500
 H                  1.03653029    3.87877644   -0.00125800
 H                 -1.45861671    3.87883644   -0.00263100
 H                 -2.70569971    1.71853244   -0.00086200
 C                  2.72398403    1.71823647    0.00088786
 C                  3.44068836    0.60271559    0.00203654
 H                  3.18050619    2.71712045    0.00051453
 H                  4.53895470    0.60279550    0.00268872
 H                  2.98416621   -0.39616839    0.00240586
''')

# 構造最適化と振動数計算
_, stylene_wfn = optfreq(stylene, 'b3lyp/cc-pvdz')
print(get_freqs(stylene_wfn))

# 結合次数の計算
_, stylene_wfn = psi4.properties('b3lyp/cc-pvdz',
                                 molecule=stylene,
                                 properties=['WIBERG_LOWDIN_INDICES'],
                                 return_wfn=True)
stylene_df = pd.DataFrame(stylene_wfn.variable('WIBERG LOWDIN INDICES').to_array())
stylene_df.index = [stylene.symbol(i) for i in range(stylene.natom())]
print(stylene_df.round(2))