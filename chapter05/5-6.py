"""アルゴン二量体のエネルギー計算
"""

from __future__ import annotations
from itertools import product
import numpy as np
import pandas as pd
import psi4


def calc_energies(level: str, geom: str, range: np.ndarray) -> list[float]:
    """
    指定範囲内で分子を定義し，指定レベルでエネルギー計算を行う
    Args:
        level: 計算レベル
        geom: 空白入りの分子座標
        range: 計算を実行する範囲

    Returns:
        list[float]: エネルギー値のリスト

    """
    energies = []
    for distance in range:
        mol = psi4.geometry(geom.format(distance))
        energy = psi4.energy(level, molecule=mol)
        energies.append(energy)

    return energies


Ar2_geom = '''
0 1
Ar  0   0   0
Ar  0   0   {}
'''

methods = ['mp2', 'b3lyp', 'b3lyp-d3bj', 'wb97x-d']
df_ar2 = pd.DataFrame()

for method in methods:
    level = f'{method}/cc-pvtz'
    df_ar2[method] = calc_energies(level=level,
                                   geom=Ar2_geom,
                                   range=np.arange(1.5, 5.5, 0.1))


print(df_ar2.round(2))