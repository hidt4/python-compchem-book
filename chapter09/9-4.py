"""ホルムアルデヒドのMulliken電荷を計算する
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


psi4.set_output_file('formaldehyde_Mulliken.log')

# ホルムアルデヒドの構造定義
HCHO = psi4.geometry('''
0 1
 C                 -0.48547719   -0.11203319    0.00000000
 O                  0.74183981   -0.11203319    0.00000000
 H                 -1.07762219    0.82737081    0.00000000
 H                 -1.07762219   -1.05143719    0.00003900
 ''')

# 構造最適化・振動数計算
_, HCHO_wfn = optfreq(HCHO, 'hf/sto-3g')
print(get_freqs(HCHO_wfn))

# Mulliken電荷の計算
psi4.oeprop(HCHO_wfn, 'MULLIKEN_CHARGES')
atoms_HCHO = [HCHO.symbol(i) for i in range(HCHO.natom())]
charges_HCHO = HCHO_wfn.variable('MULLIKEN CHARGES')

for atom, charge in zip(atoms_HCHO, charges_HCHO):
    print(f'{atom}\t{charge:.3f}')