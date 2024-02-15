"""ベンゼンのOPAスペクトル
"""

from __future__ import annotations
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from psi4.driver.p4util import spectrum


def optfreq(mol: psi4.core.Molecule,
            theory: str) -> list[float, psi4.core.Wavefunction]:
    """
    分子の構造最適化と振動数計算を行う

    Args:
        mol (psi4.core.Molecule): 計算対象の分子
        theory: 計算レベル

    Returns:
        float: エネルギー
        psi4.core.Wavefunction: 計算のWavefunction

    """
    _, wfn = psi4.optimize(theory,
                           molecule=mol,
                           return_wfn=True)
    energy, wfn = psi4.frequency(theory,
                                 molecule=mol,
                                 ref_gradient=wfn.gradient(),
                                 return_wfn=True)

    return [energy, wfn]


def get_freqs(wavefunction: psi4.core.Wavefunction) -> np.ndarray:
    """
    Wavefunctionオブジェクトから振動数を取り出す
    Args:
        wavefunction: 対象のWavefunctionオブジェクト

    Returns:
        振動数のarrray

    """
    return wavefunction.frequencies().to_array()


mpl.style.use('seaborn-v0_8-poster')
mpl.style.use('seaborn-v0_8-whitegrid')

psi4.set_num_threads(14)
psi4.set_memory('20GB')
psi4.set_options({'SAVE_JK': True})
psi4.set_output_file('OPA_benzene.log')

benzene = psi4.geometry('''
0 1
 C                 -0.45228216    0.14522821    0.00000000
 C                  0.94287784    0.14522821    0.00000000
 C                  1.64041584    1.35297921    0.00000000
 C                  0.94276184    2.56148821   -0.00119900
 C                 -0.45206316    2.56141021   -0.00167800
 C                 -1.14966416    1.35320421   -0.00068200
 H                 -1.00204116   -0.80708879    0.00045000
 H                  1.49238584   -0.80728479    0.00131500
 H                  2.74009584    1.35305921    0.00063400
 H                  1.49296184    3.51363121   -0.00125800
 H                 -1.00218516    3.51369121   -0.00263100
 H                 -2.24926816    1.35338721   -0.00086200
''')

# 構造最適化と振動数計算
_, wfn_c6h6 = optfreq(benzene, theory='cam-b3lyp/cc-pvdz')
print(benzene.save_string_xyz())
print(get_freqs(wfn_c6h6))

# TD-DFT計算
_, wfn_td = psi4.energy('cam-b3lyp/aug-cc-pvdz',
                        molecule=benzene,
                        return_wfn=True)
excited_states = tdscf_excitations(wfn_td, states=4)

# OPAスペクトル
poles = [state['EXCITATION ENERGY'] for state in excited_states]
opa_residues = [np.linalg.norm(state['ELECTRIC DIPOLE TRANSITION MOMENT (LEN)']) ** 2
                for state in excited_states]
opa_spectrum = spectrum(poles=poles,
                        residues=opa_residues,
                        gamma=0.01,
                        out_units='nm')

# 描画
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(opa_spectrum['convolution']['x'],
        opa_spectrum['convolution']['y'],
        color='black')
ax.bar(opa_spectrum['sticks']['poles'],
       opa_spectrum['sticks']['residues'],
       color='black')
ax.set_xlabel('Wavelength (nm)', size=20)
ax.set_ylabel(r'$ \rm{\epsilon}$', size=20)
plt.show()
