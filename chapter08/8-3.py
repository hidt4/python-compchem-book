"""水分子のエネルギー計算
HF/STO-3Gレベルでcubeファイルを出力する
"""

import pathlib
import psi4

h2o = psi4.geometry('''
O       -0.1176269719      0.7387773605      0.0000000000
H        0.8523730281      0.7387773605      0.0000000000
H       -0.4409567836      1.4770262439     -0.5397651517
''')

psi4.set_output_file('h2o_sto3g_cubeprop.log')
_, wfn = psi4.energy('hf/sto-3g', molecule=h2o, return_wfn=True)

# cubeファイルの保存先を作成
target_dir = pathlib.Path('h2o_cube')
target_dir.mkdir(exist_ok=True)

# cubeファイルの作成
psi4.set_options({'CUBEPROP_TASKS': ['ORBITALS'],
                 'CUBEPROP_FILEPATH': target_dir.name})
psi4.cubeprop(wfn)