"""エタンの構造最適化・振動数計算
"""

import psi4

psi4.set_output_file('ethane_optfreq_hf-3-21g.log')
psi4.set_options({'normal_modes_write': True})

# エタン構造の定義
ethane = psi4.geometry('''
0 1
 C                 -1.25726143    1.25726139    0.00000000
 H                 -0.90060700    0.24845139    0.00000000
 H                 -0.90058859    1.76165958   -0.87365150
 H                 -2.32726143    1.25727458    0.00000000
 C                 -0.74391921    1.98321767    1.25740497
 H                  0.32608077    1.98303509    1.25750243
 H                 -1.10041427    2.99208399    1.25730739
 H                 -1.10075147    1.47893186    2.13105625
 ''')

# 構造最適化
_, wfn_c2h6 = psi4.optimize('hf/3-21g',
                            molecule=ethane,
                            return_wfn=True)
print(ethane.save_string_xyz_file())

# 振動数計算
_, wfn_c2h6 = psi4.frequency('hf/3-21g',
                             molecule=ethane,
                             return_wfn=True,
                             ref_gradient=wfn_c2h6.gradient())
print(wfn_c2h6.frequencies().to_array())