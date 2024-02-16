"""反応座標軸に沿った構造のスキャン
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
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
psi4.set_output_file('sn2_scan.log')

# Z座標を変数として扱えるように構造を定義
sn2_geom = '''
 -1 1
 C                 -0.37149818    0.23751522    0.00000000
 H                  0.13279358   -0.63613612    0.35682369
 H                 -1.38041486    0.23751522    0.35635256
 H                  0.13279358    1.11116656    0.35682369
 Cl                -0.37090351    0.23751522   -1.90999991
 Cl                -0.37090351    0.23751522   {}
'''

energies = []
geometries = []
distance = np.arange(1.7, 2.8, 0.1)
psi4.set_options({'GEOM_MAXITER': 200,
                  'FROZEN_DISTANCE': '1 6',
                  'INTRAFRAG_STEP_LIMIT_MAX': 0.5,
                  'OPT_COORDINATES': 'CARTESIAN'})

# 距離を変えながら構造最適化
for dist in distance:
    ts = psi4.geometry(sn2_geom.format(dist))
    energy = psi4.optimize('b3lyp/6-31g(d)', molecule=ts)
    energies.append(energy)
    geometries.append(ts.save_string_xyz())

# エネルギー最大値を与える距離を出力
idx = np.argmax(energies)
print(f'maximum energy at {distance[idx]: .1f} angstrom')

# 遷移状態の最適化
psi4.core.clean_options()
psi4.set_options({'OPT_TYPE': 'TS',
                  'GEOM_MAXITER': 100,
                  'OPT_COORDINATES': 'CARTESIAN',
                  'INTRAFRAG_STEP_LIMIT_MAX': 0.5,
                  'FULL_HESS_EVERY': 0,
                  'normal_modes_write': True})

ts = psi4.geometry(geometries[idx])
psi4.set_output_file('SN2_OPT_TS.log')

# 構造最適化・振動数計算
_, wfn_sn2 = optfreq(ts, 'b3lyp/6-31g(d)')
# 振動数の出力
print(get_freqs(wfn_sn2))
# 最適化構造の保存
ts.save_xyz_file('sn2_ts.xyz', True)

# エネルギー図の出力
plt.plot(distance, energies)
plt.show()
