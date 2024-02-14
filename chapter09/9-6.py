"""ベンゼンの双極子モーメント
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


psi4.set_output_file('benzene_dipole.log')
au2debye = psi4.constants.dipmom_au2debye

# ベンゼンの構造定義
benzene = psi4.geometry('''
0 1
 C                 -0.93360997    0.57676348    0.00000000
 C                  0.46155003    0.57676348    0.00000000
 C                  1.15908803    1.78451448    0.00000000
 C                  0.46143403    2.99302348   -0.00119900
 C                 -0.93339097    2.99294548   -0.00167800
 C                 -1.63099197    1.78473948   -0.00068200
 H                 -1.48336897   -0.37555352    0.00045000
 H                  1.01105803   -0.37574952    0.00131500
 H                  2.25876803    1.78459448    0.00063400
 H                  1.01163403    3.94516648   -0.00125800
 H                 -1.48351297    3.94522648   -0.00263100
 H                 -2.73059597    1.78492248   -0.00086200
''')

# 構造最適化・振動数計算
_, benzene_wfn = optfreq(benzene, 'mp2/cc-pvdz')
print(get_freqs(benzene_wfn))

# 双極子モーメントの取得
benzene_dipole = benzene_wfn.variable('SCF DIPOLE')
benzene_dipole_au = np.sqrt(np.sum(benzene_dipole ** 2))
print(f'benzene dipole moment: {au2debye * benzene_dipole_au: .2f} D')
