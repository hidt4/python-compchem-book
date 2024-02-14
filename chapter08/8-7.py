"""エチレンの分子軌道を基底関数系を変えて見てみる
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


def get_fchk(molecule: psi4.core.Molecule, basis: str):
    """
    分子のエネルギー計算を行いfchkファイルを保存する
    Args:
        molecule: 計算対象の分子
        basis: 使用する基底関数系

    """
    theory = f'hf/{basis}'
    _, wfn = psi4.energy(theory, molecule=molecule, return_wfn=True)
    psi4.fchk(wfn, f'hf_{basis}.fchk')


psi4.set_num_threads(14)
psi4.set_memory('32GB')
psi4.set_output_file('ethylene_basisset.log')
basis_sets = ['sto-3g', '3-21g', '6-31+g(d,p)', '6-311+g(2d,p)']

# エチレンの構造を定義
ethylene = psi4.geometry('''
0 1
 C                 -0.32780083    1.53112031    0.00000000
 C                  0.99811517    1.53112031    0.00000000
 H                 -0.92138583    2.45515831    0.00000000
 H                 -0.92141683    0.60710631   -0.00002200
 H                  1.59170017    0.60708231   -0.00001900
 H                  1.59173117    2.45513431    0.00002600
''')

# 構造最適化と振動数計算
_, wfn_ethylene = optfreq(ethylene, theory='hf/sto-3g')
print(get_freqs(wfn_ethylene))

# 各レベルでエネルギー計算
for basis in basis_sets:
    psi4.set_output_file(f'ethylene_SP_{basis}.log')
    get_fchk(molecule=ethylene, basis=basis)

