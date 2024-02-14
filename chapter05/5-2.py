"""色々な手法で水分子のエネルギー計算
"""

import psi4

h2o = psi4.geometry('''
0 1
O       -0.1176269719      0.7387773605      0.0000000000
H        0.8523730281      0.7387773605      0.0000000000
H       -0.4409567836      1.4770262439     -0.5397651517
''')

psi4.set_memory('16GB')
psi4.set_num_threads(12)

# HF/STO-3G
psi4.set_output_file('HF_STO-3G.log')
psi4.energy('hf/sto-3g', molecule=h2o)

# MP2/STO-3G
psi4.set_output_file('MP2_STO-3G.log')
psi4.energy('mp2/sto-3g', molecule=h2o)

# B3LYP/STO-3G
psi4.set_output_file('B3LYP_STO-3G.log')
psi4.energy('b3lyp/sto-3g', molecule=h2o)

# HF/6-31G
psi4.set_output_file('HF_6-31G.log')
psi4.energy('hf/6-31g', molecule=h2o)

# HF/6-31G(d)
psi4.set_output_file('HF_6-31Gd.log')
psi4.energy('hf/6-31g(d)', molecule=h2o)

# HF/6-31G(d,p)
psi4.set_output_file('HF_6-31Gdp.log')
psi4.energy('hf/6-31g(d,p)', molecule=h2o)

# HF/6-311G
psi4.set_output_file('HF_6-311G.log')
psi4.energy('hf/6-311g', molecule=h2o)
