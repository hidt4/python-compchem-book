"""エチレンの構造最適化・振動数計算
"""

import psi4

psi4.set_output_file('ethylene_optfreq_hf-3-21g.log')

# エチレンの構造定義
ethylene = psi4.geometry('''
0 1
 C                 -0.31950208    0.90041492    0.00000000
 C                  1.00641392    0.90041492    0.00000000
 H                 -0.91308708    1.82445292    0.00000000
 H                 -0.91311808   -0.02359908   -0.00002200
 H                  1.59999892   -0.02362308   -0.00001900
 H                  1.60002992    1.82442892    0.00002600
''')

# 構造最適化
_, wfn_c2h4 = psi4.optimize('hf/3-21g',
                            molecule=ethylene,
                            return_wfn=True)

# 振動数計算
psi4.set_options({'normal_modes_write': True})
_, wfn_c2h4 = psi4.frequency('hf/3-21g',
                             molecule=ethylene,
                             return_wfn=True,
                             ref_gradient=wfn_c2h4.gradient())
print(wfn_c2h4.frequencies().to_array())