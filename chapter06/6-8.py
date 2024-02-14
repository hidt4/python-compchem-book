"""frozen_dihedralオプション
"""
import psi4

psi4.set_num_threads(14)
psi4.set_memory('20GB')

# ログファイルの設定
psi4.set_output_file('frozen_dihedral.log')

# 二面角を90度にして初期構造を定義
dce_constrained = psi4.geometry('''
0 1
 C                 -2.25311207    1.17427384    0.00000000
 H                 -1.69945761    1.74831286   -0.71333438
 H                 -3.29070391    1.18574460   -0.26109690
 C                 -1.73979635   -0.27765831    0.00000000
 H                 -2.29345264   -0.85169798   -0.71333244
 H                 -1.86707706   -0.70095780    0.97443171
 Cl                -2.04374927    1.87054265    1.60280285
 Cl                -0.03310419   -0.29652555   -0.42947162
''')

# 最適化
print('### Constrained optimization')
psi4.set_options({'FROZEN_DIHEDRAL': '7 1 4 8'})
print(psi4.p4util.prepare_options_for_set_options())
psi4.optimize('hf/sto-3g', molecule=dce_constrained)
print(dce_constrained.save_string_xyz())

print('### Clear constraints')
psi4.core.clean_options()
print(psi4.p4util.prepare_options_for_set_options())
psi4.optimize('hf/sto-3g', molecule=dce_constrained)
print(dce_constrained.save_string_xyz())