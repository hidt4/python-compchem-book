"""水分子の構造最適化
"""

import psi4

# 計算のログファイルを設定
psi4.set_output_file('h2o_opt_hf-sto3g.log')

# 水分子の構造を定義
h2o = psi4.geometry('''
0 1
O       -0.7520847362      0.3573098626      0.0156794168
H        0.2372416200      0.3907947487     -0.0110698848
H       -1.0414501312      1.0972341870     -0.5754069255
''')

# 構造最適化計算の実行
psi4.optimize('hf/sto-3g', molecule=h2o)

# 最適化構造を出力
print(h2o.save_string_xyz())