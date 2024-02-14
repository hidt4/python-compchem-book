"""水分子のエネルギー計算
"""

import sys
import numpy as np
import psi4

# 計算環境の設定
psi4.set_num_threads(6)
psi4.set_memory('40GB')
# ログファイルの設定
psi4.set_output_file('h2o_optfreq_sto3g.log')

# ライブラリのバージョン
print(f'python version:\n{sys.version}')
print(f'numpy version: {np.__version__}')
print(f'Psi4 version: {psi4.__version__}')


# 水分子の構造定義
h2o = psi4.geometry('''
0 1
O       -0.7520847362      0.3573098626      0.0156794168
H        0.2372416200      0.3907947487     -0.0110698848
H       -1.0414501312      1.0972341870     -0.5754069255
''')

# 構造最適化計算
psi4.optimize('hf/sto-3g', molecule=h2o)

# 振動数計算
_, wfn = psi4.frequency('hf/sto-3g',
                        molecule=h2o,
                        return_wfn=True)
print(wfn.frequencies().to_array())