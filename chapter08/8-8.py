"""軌道エネルギーを比較する
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


psi4.set_memory('32GB')
psi4.set_num_threads(14)

# アクリロニトリル
psi4.set_output_file('acrylonitrile.log')
acrylonitrile = psi4.geometry('''
0 1
 C                 -0.54356847    0.29460580    0.00000000
 C                  0.78234753    0.29460580    0.00000000
 H                 -1.13715347    1.21864380    0.00000000
 H                 -1.13718447   -0.62940820   -0.00002200
 H                  1.37593253   -0.62943220   -0.00001900
 C                  1.61472425    1.59027127    0.00003646
 N                  2.24170932    2.56622707    0.00006392
 ''')

_, wfn_acrylonitrile = optfreq(acrylonitrile, 'b3lyp/6-31g(d)')
print(get_freqs(wfn_acrylonitrile))

_, wfn_acrylonitrile = psi4.energy('hf/6-311+g(2d,p)',
                                   molecule=acrylonitrile,
                                   return_wfn=True)
psi4.fchk(wfn_acrylonitrile, 'acrylonitriele.fchk')


# 1,1-ジシアノエチレン
psi4.set_output_file('dicyano.log')
dicyano = psi4.geometry('''
0 1
 C                 -0.56016598   -0.13692946    0.00000000
 C                  0.76575002   -0.13692946    0.00000000
 C                  1.59812674    1.15873600    0.00003646
 N                  2.22511181    2.13469180    0.00006392
 C                  1.59808067   -1.43262452   -0.00002664
 N                  2.22503103   -2.40860261   -0.00004671
 H                 -1.09516598    0.78971772    0.00002206
 H                 -1.09516598   -1.06357664   -0.00002206
''')

_, wfn_di = optfreq(dicyano, 'b3lyp/6-31g(d)')
print(get_freqs(wfn_di))

_, wfn_di = psi4.energy('hf/6-311+g(2d,p)',
                        molecule=dicyano,
                        return_wfn=True)
psi4.fchk(wfn_di, 'dicyano.fchk')


# トリシアノエチレン
psi4.set_output_file('tricyano.log')
tricyano = psi4.geometry('''
0 1
 C                 -0.56016598   -0.13692946    0.00000000
 C                  0.76575002   -0.13692946    0.00000000
 C                  1.59812674    1.15873600    0.00003646
 N                  2.22511181    2.13469180    0.00006392
 C                 -1.39254271   -1.43259492   -0.00003085
 N                 -2.01952778   -2.40855072   -0.00005409
 C                  1.59808067   -1.43262452   -0.00002664
 N                  2.22503103   -2.40860261   -0.00004671
 H                 -1.07303856    0.80214547    0.00002236
 ''')

_, wfn_tri = optfreq(tricyano, 'b3lyp/6-31g(d)')
print(get_freqs(wfn_tri))

_, wfn_tri = psi4.energy('hf/6-311+g(2d,p)',
                         molecule=tricyano,
                         return_wfn=True)
psi4.fchk(wfn_tri, 'tricyano.fchk')


# テトラシアノエテン
psi4.set_output_file('tetracyano.log')
tetracyano = psi4.geometry('''
0 1
 C                 -0.56016598   -0.13692946    0.00000000
 C                  0.76575002   -0.13692946    0.00000000
 C                  1.59812674    1.15873600    0.00003646
 N                  2.22511181    2.13469180    0.00006392
 C                 -1.39249664    1.15876560    0.00000000
 N                 -2.01944700    2.13474370   -0.00000000
 C                 -1.39254271   -1.43259492   -0.00003085
 N                 -2.01952778   -2.40855072   -0.00005409
 C                  1.59808067   -1.43262452   -0.00002664
 N                  2.22503103   -2.40860261   -0.00004671
 ''')

_, wfn_tetra = optfreq(tetracyano, 'b3lyp/6-31g(d)')
print(get_freqs(wfn_tetra))

_, wfn_tetra = psi4.energy('hf/6-311+g(2d,p)',
                           molecule=tetracyano,
                           return_wfn=True)
psi4.fchk(wfn_tetra, 'tetracyano.fchk')
