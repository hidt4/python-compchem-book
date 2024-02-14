"""水分子のエネルギー計算
異なる基底関数系を用いて計算を行い，
得られる分子軌道の数が異なることを確認する
"""

import psi4

h2o = psi4.geometry('''
O       -0.1176269719      0.7387773605      0.0000000000
H        0.8523730281      0.7387773605      0.0000000000
H       -0.4409567836      1.4770262439     -0.5397651517
''')

basis_sets = ['sto-3g', '6-31g(d)']

for basis_set in basis_sets:
    psi4.set_output_file(f'h2o_SP_{basis_set}.log')
    psi4.energy(f'hf/{basis_set}', molecule=h2o)