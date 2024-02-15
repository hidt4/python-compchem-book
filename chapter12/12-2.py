"""ホルムアルデヒドのTD-DFT計算
"""

from __future__ import annotations
import numpy as np
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations


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


psi4.set_num_threads(14)
psi4.set_memory('20GB')
psi4.set_options({'SAVE_JK': True})
psi4.set_output_file('TD_formaldehyde2.log')

HCHO = psi4.geometry('''
0 1
 C                 -0.02074689    2.31120329    0.00000000
 O                  1.20657011    2.31120329    0.00000000
 H                 -0.61289189    3.25060729    0.00000000
 H                 -0.61289189    1.37179929    0.00003900
''')

# 構造最適化と振動数計算
_, wfn_hcho = optfreq(HCHO, theory='cam-b3lyp/cc-pvdz')
print(get_freqs(wfn_hcho))
# エネルギー計算
_, wfn_td = psi4.energy('CAM-B3LYP/aug-cc-pvdz',
                        molecule=HCHO,
                        return_wfn=True)

# 励起状態の計算
excited_states = tdscf_excitations(wfn_td, states=4)

# 第1励起状態の出力
state1 = excited_states[0]
for key in state1.keys():
    print(key)