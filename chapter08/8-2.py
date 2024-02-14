"""水分子のエネルギー計算
"""

import psi4

h2o = psi4.geometry('''
O       -0.1176269719      0.7387773605      0.0000000000
H        0.8523730281      0.7387773605      0.0000000000
H       -0.4409567836      1.4770262439     -0.5397651517
''')

psi4.set_output_file('h2o_sto3g_fchk.log')
_, wfn = psi4.energy('hf/sto-3g', molecule=h2o, return_wfn=True)
psi4.fchk(wfn, 'h2o_sto3g.fchk')