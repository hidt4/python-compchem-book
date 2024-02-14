"""ピラジンのESPを表示する
"""


from __future__ import annotations
import numpy as np
import psi4


def optfreq(mol: psi4.core.Molecule,
            theory: str )-> list[float, psi4.core.Wavefunction]:
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
psi4.set_memory('20GB')

psi4.set_output_file('pyrazine_esp.log')

# ピラジンの構造定義
pyrazine = psi4.geometry('''
0 1
 C                 -1.03319504    0.52697095    0.00000000
 C                  0.38637296    0.52697095    0.00000000
 N                  1.09861396    1.67230695    0.00000000
 C                  0.38631896    2.81756395   -0.00004800
 C                 -1.03327004    2.81757695   -0.00006100
 N                 -1.74546104    1.67226695   -0.00006200
 H                 -1.60724804   -0.41446405    0.00005200
 H                  0.96039596   -0.41446405    0.00003600
 H                  0.96034996    3.75900795   -0.00006900
 H                 -1.60729804    3.75903195   -0.00011100
 ''')

# 構造最適化と振動数計算
_, pyrazine_wfn = optfreq(pyrazine, 'b3lyp/cc-pvdz')
print(get_freqs(pyrazine_wfn))

# 最適化構造でfchkファイル生成のWavefunctionを取得
_, pyrazine_wfn= psi4.energy('b3lyp/cc-pvdz',
                             molecule=pyrazine,
                             return_wfn=True)
psi4.fchk(pyrazine_wfn, 'pyrazine.fchk')