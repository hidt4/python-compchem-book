"""電荷の基底関数系依存性の調査
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import psi4


basis_sets = ['sto-3g', '3-21g', 'cc-pvdz', 'aug-cc-pvdz', 'cc-pvtz',
              'aug-cc-pvtz', 'cc-pvqz', 'aug-cc-pvqz', 'aug-cc-pv5z']
properties = ['MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'MBIS_CHARGES']

psi4.set_output_file('basisset_atomic_charges.log')

# メタノールの構造定義
methanol = psi4.geometry('''
0 1
 C                 -1.39464072    0.45676004    0.00000000
 H                 -0.94297547   -0.51123296   -0.06235416
 H                 -1.13979098    1.02593296   -0.86947904
 H                 -2.45782841    0.35118152    0.05818167
 O                 -0.91796581    1.13086230    1.16759033
 H                 -1.32295402    1.99942959    1.22394109
 ''')

# 電荷を格納するデータフレームの準備
Mulliken = pd.DataFrame()
Lowdin = pd.DataFrame()
MBIS = pd.DataFrame()
indices = [methanol.symbol(i) for i in range(methanol.natom())]
Mulliken.index = indices
Lowdin.index = indices
MBIS.index = indices

# 各レベルで電荷の計算
for basis_set in basis_sets:
    _, wfn = psi4.properties(f'b3lyp/{basis_set}',
                             molecule=methanol,
                             properties=properties,
                             return_wfn=True)
    Mulliken[basis_set] = wfn.variable('MULLIKEN CHARGES')
    Lowdin[basis_set] = wfn.variable('LOWDIN CHARGES')
    MBIS[basis_set] = wfn.variable('MBIS CHARGES').flatten()

# 結果の出力
print(Mulliken.round(3))
print(Lowdin.round(3))
print(MBIS.round(3))
