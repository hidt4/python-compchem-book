"""エチレンのMulliken電荷を求める
"""

from __future__ import annotations
import numpy as np
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


psi4.set_output_file('ethylene_Mulliken.log')

# エチレンの構造定義
ethylene = psi4.geometry('''
0 1
 C                 -0.97510375    1.11618256    0.00000000
 C                  0.35081225    1.11618256    0.00000000
 H                 -1.56868875    2.04022056    0.00000000
 H                 -1.56871975    0.19216856   -0.00002200
 H                  0.94439725    0.19214456   -0.00001900
 H                  0.94442825    2.04019656    0.00002600
 ''')

# 構造最適化・振動数計算
_, ethylene_wfn = optfreq(ethylene, 'hf/sto-3g')
print(get_freqs(ethylene_wfn))

# Mulliken電荷の計算
psi4.oeprop(ethylene_wfn, 'MULLIKEN_CHARGES')
atoms_ethylene = [ethylene.symbol(i) for i in range(ethylene.natom())]
charges_ethylene = ethylene_wfn.variable('MULLIKEN CHARGES')

for atom, charge in zip(atoms_ethylene, charges_ethylene):
    print(f'{atom}\t{charge:.3f}')