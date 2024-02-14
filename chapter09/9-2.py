"""psi4.propertiesとpsi4.oepropの比較
"""

import time
import psi4

psi4.set_output_file('properties_oeprop.log')

# ベンゼンの構造定義
benzene = psi4.geometry('''
0 1
 C                 -0.93360997    0.57676348    0.00000000
 C                  0.46155003    0.57676348    0.00000000
 C                  1.15908803    1.78451448    0.00000000
 C                  0.46143403    2.99302348   -0.00119900
 C                 -0.93339097    2.99294548   -0.00167800
 C                 -1.63099197    1.78473948   -0.00068200
 H                 -1.48336897   -0.37555352    0.00045000
 H                  1.01105803   -0.37574952    0.00131500
 H                  2.25876803    1.78459448    0.00063400
 H                  1.01163403    3.94516648   -0.00125800
 H                 -1.48351297    3.94522648   -0.00263100
 H                 -2.73059597    1.78492248   -0.00086200
''')

# psi4.properties
start = time.time()
_, wfn = psi4.properties('hf/6-311+g(2d,p)',
                         molecule=benzene,
                         properties=['DIPOLE', 'MULLIKEN_CHARGES'],
                         return_wfn=True)
end = time.time()
print(f'psi4.properties: {end - start: .2f} sec')

# psi4.oeprop
_, wfn = psi4.energy('hf/6-311+g(2d,p)', molecule=benzene, return_wfn=True)
start = time.time()
psi4.oeprop(wfn, 'DIPOLE', 'MULLIKEN_CHARGES')
end = time.time()
print(f'psi4.oeprop: {end - start: .2f} sec')