"""1,2-ジクロロエタンの遷移状態
"""

from __future__ import annotations
import numpy as np
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


psi4.set_num_threads(14)
psi4.set_memory('28GB')
psi4.set_output_file('dce_ts.log')

# 二面角0度の初期構造
dce = psi4.geometry('''
0 1
 C                 -2.52697099    1.46473027    0.00000000
 H                 -2.17029815    1.96912846   -0.87365150
 H                 -3.59697099    1.46474345    0.00000000
 C                 -2.01365528    0.01279812    0.00000000
 H                 -2.84588587   -0.65972864    0.00000086
 H                 -1.41921341   -0.15534323   -0.87365174
 Cl                -1.94029417    2.29439458    1.43703425
 Cl                -1.03587983   -0.26377103    1.43703329
 ''')

# 遷移状態計算に設定
psi4.set_options({'opt_type': 'TS',
                  'normal_modes_write': True})

# 構造最適化・振動数計算
_, wfn_dce = optfreq(dce, 'hf/sto-3g')

print(get_freqs(wfn_dce))