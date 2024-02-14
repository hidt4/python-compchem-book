"""水分子の振動数計算
"""

import numpy as np
import psi4


# 水分子の構造の設定
h2o = psi4.geometry('''
0 1
O       -0.7520847362      0.3573098626      0.0156794168
H        0.2372416200      0.3907947487     -0.0110698848
H       -1.0414501312      1.0972341870     -0.5754069255
''')

# 振動数オプションの設定
psi4.set_output_file('h2o_optfreq_sto3g_2.log')
psi4.set_options({'normal_modes_write': True})

# 構造最適化計算
psi4.optimize('hf/sto-3g', molecule=h2o)

# 振動数計算
_, wfn = psi4.frequency('hf/sto-3g',
                        molecule=h2o,
                        return_wfn=True)
print(wfn.frequencies().to_array())