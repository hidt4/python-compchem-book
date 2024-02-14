"""塩化ナトリウムのエネルギー計算
"""

from __future__ import annotations
from itertools import product
import numpy as np
import pandas as pd
import psi4


def calc_energies(level: str,
                  geom: str,
                  range: np.ndarray) -> list[float]:
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


NaCl_geom = '''
0 1
Na  0   0   0
Cl  0   0   {}
'''

methods = ['hf', 'mp2', 'b3lyp']
basis_sets = ['sto-3g', 'cc-pvdz', 'cc-pvtz', 'aug-cc-pvtz']
df_NaCl = pd.DataFrame()

for (method, basis_set) in product(methods, basis_sets):
    level = method + '/' + basis_set
    df_NaCl[level] = calc_energies(level=level,
                                   geom=NaCl_geom,
                                   range=np.arange(1.0, 4.0, 0.1))


print(df_NaCl.round(2))