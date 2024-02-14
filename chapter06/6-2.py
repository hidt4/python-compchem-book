"""異なる終了判定基準で構造最適化
"""

import psi4

# ログファイルを設定
psi4.set_output_file('h2o_opt_hf-sto3g_2.log')

h2o = psi4.geometry('''
0 1
 O   -0.000000000000    0.000000000000   -0.071153222149
 H    0.758023512319    0.000000000000    0.564626634974
 H   -0.758023512319   -0.000000000000    0.564626634974
''')

# 構造最適化のオプション設定
psi4.set_options({'G_CONVERGENCE': 'GAU_VERYTIGHT'})

# 構造最適化
psi4.optimize('hf/sto-3g', molecule=h2o)

# 最適化構造を出力
print(h2o.save_string_xyz())