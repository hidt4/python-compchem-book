"""IRC計算
"""

import shutil
import psi4

psi4.set_num_threads(14)
psi4.set_memory('28GB')

ts = psi4.geometry('''
-1 1
 C    0.264926007988   -0.006709497680    0.000000000000
 H   -0.271775267024   -0.031250611390   -0.930131463207
 H    1.338116752186    0.041860592522    0.000000000000
 H   -0.271775267024   -0.031250611390    0.930131463207
CL    0.372760471437   -2.368819184104    0.000000000000
CL    0.157128730252    2.354476747149    0.000000000000
''')

# IRC計算用に分子をコピー
forward = ts.clone()

# forward方向のIRC計算
psi4.set_options({
    'OPT_TYPE': 'IRC',
    'IRC_DIRECTION': 'FORWARD',
    'IRC_POINTS': 20,
    'IRC_STEP_SIZE': 0.2,
    'FULL_HESS_EVERY': 0,
    'OPT_COORDINATES': 'CARTESIAN',
    'GEOM_MAXITER': 200
                  })

psi4.set_output_file('sn2_irc_forward.log')
psi4.optimize('b3lyp/6-31g(d)', molecule=forward)

# 計算結果の保存
forward.save_xyz_file('forward.xyz', True)
shutil.copy('ircprogress.log', 'forward.log')