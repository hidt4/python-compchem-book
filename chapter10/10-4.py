"""簡単な遷移状態を鋳型として複雑な遷移状態を作成
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

ts_2 = psi4.geometry('''
-1 1
C        0.2649260080     -0.0067094977      0.0000000000                 
H       -0.2717752670     -0.0312506114     -0.9301314632                 
C        1.3381167522      0.0418605925      0.0000000000                 
C       -0.2717752670     -0.0312506114      0.9301314632                 
Cl       0.3727604714     -2.3688191841      0.0000000000                 
Cl       0.1571287303      2.3544767471      0.0000000000                 
H        1.6944025173      0.0582577786     -1.0088070112                 
H        1.7339248047     -0.8149079124      0.5041663560                 
H        1.6549282236      0.9306064646      0.5046397392                 
H       -1.3226464201     -0.0787576909      0.7343946575                 
H       -0.0531535919      0.8530376826      1.4915036894                 
H        0.0258452199     -0.8924765946      1.4910341746                 
''')

# 制約付き最適化
psi4.set_options({
    'GEOM_MAXITER': 200,
    'FROZEN_CARTESIAN': '1 XYZ 5 XYZ 6 XYZ'
})
psi4.set_output_file('SN2_OPT2.log')
psi4.optimize('b3lyp/6-31g(d)', molecule=ts_2)

# 遷移状態の最適化
psi4.core.clean_options()
psi4.set_options({
    'OPT_TYPE': 'TS',
    'GEOM_MAXITER': 100,
    'OPT_COORDINATES': 'CARTESIAN',
    'INTRAFRAG_STEP_LIMIT_MAX': 0.5,
    'FULL_HESS_EVERY': 0,
    'NORMAL_MODES_WRITE': True,
})
psi4.set_output_file('SN2_OPT_TS2.log')
_, wfn_sn2_2 = optfreq(ts_2, 'b3lyp/6-31g(d)')
print(get_freqs(wfn_sn2_2))
print(ts_2.save_string_xyz())
ts_2.save_xyz_file('sn2_ts2.xyz', True)