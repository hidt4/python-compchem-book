"""アセチレンの構造最適化・振動数計算
"""

import psi4

psi4.set_output_file('acetylene_optfreq_hf-3-21g.log')

# アセチレンの構造定義
acetylene = psi4.geometry('''
0 1
 C                 -0.51037345    0.81742737    0.00000000
 C                  0.68462655    0.81742737    0.00000000
 H                 -1.57137345    0.81742737    0.00000000
 H                  1.74562655    0.81742737    0.00000000
''')

# 構造最適化
_, wfn_c2h2 = psi4.optimize('hf/3-21g',
                            molecule=acetylene,
                            return_wfn=True)

# 振動数計算
psi4.set_options({'normal_modes_write': True})
_, wfn_c2h2 = psi4.frequency('hf/3-21g',
                             molecule=acetylene,
                             ref_gradient=wfn_c2h2.gradient(),
                             return_wfn=True)

# 振動数を出力
print(wfn_c2h2.frequencies().to_array())