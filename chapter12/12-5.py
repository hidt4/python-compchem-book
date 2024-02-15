"""計算する励起状態の数を変えてECDを比較
"""

from __future__ import annotations
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations


def get_ECD_spectrum(excited_states: list) -> list[np.ndarray]:
    """
    ECDスペクトル描画用データを作成する
    Args:
        excited_state: TD計算で得られた励起状態のデータ

    Returns:
        np.ndarray: ECDスペクトル用のX軸
        np.ndarray: ECDスペクトル用のY軸
        np.ndarray: 吸収波長
        np.ndarray: 強度

    """
    from psi4.driver.p4util import spectrum

    poles = [excited_state['EXCITATION ENERGY'] for excited_state in excited_states]
    residues = [excited_state['ROTATORY STRENGTH (LEN)'] for excited_state in excited_states]
    ecd_spectrum = spectrum(poles=poles,
                            residues=residues,
                            kind='ECD',
                            gamma=0.01,
                            out_units='nm')

    x_broaden = ecd_spectrum['convolution']['x']
    y_broaden = ecd_spectrum['convolution']['y']
    x = ecd_spectrum['sticks']['poles']
    y = ecd_spectrum['sticks']['residues']

    return [x_broaden, y_broaden, x, y]


mpl.style.use('seaborn-v0_8-poster')
mpl.style.use('seaborn-v0_8-whitegrid')

psi4.set_num_threads(14)
psi4.set_memory('20GB')
psi4.set_options({'SAVE_JK': True})
psi4.set_output_file('ECD_comparison.log')

s_form = psi4.geometry('''
0 1
 C   -1.914289575658   -1.037107371791    0.213196572292
 C   -0.521700028002   -1.056355678596    0.249244926680
 C    0.202359881702    0.128438076284    0.122685235042
 C   -0.489832386012    1.332134716868   -0.029009849201
 C   -1.880108537841    1.352172053526   -0.066401486859
 C   -2.598238855598    0.164184334968    0.053162898795
 H   -2.469854534141   -1.971605397783    0.315094831987
 H    0.021037262963   -1.990658201289    0.384959539471
 H    0.066075968845    2.269765118799   -0.112999689290
 H   -2.406504489448    2.301395070124   -0.183709491535
 H   -3.689400725785    0.177289946206    0.027304360256
 C    1.716525479583    0.121666923384    0.124676609167
 H    2.058136911594    0.949865748309    0.776784243561
 C    2.282222905013    0.339204791582   -1.276611019725
 H    3.384780697299    0.355004184935   -1.252337347666
 H    1.957049412976   -0.474288800419   -1.942035484724
 H    1.942420215434    1.294831711290   -1.701746476663
 O    2.152479713756   -1.119290971161    0.653722133977
 H    3.113144187560   -1.154439846669    0.565567660301
''')

num_states = [4, 8, 20, 40]

fig = plt.figure()
for i, num in enumerate(num_states, 1):
    _, wfn_td = psi4.energy('cam-b3lyp/aug-cc-pvdz',
                            molecule=s_form,
                            return_wfn=True)
    excited_states = tdscf_excitations(wfn_td, states=num)
    x_broaden, y_broaden, x, y  = get_ECD_spectrum(excited_states)
    ax = fig.add_subplot(2, 2, i)
    ax.plot(x_broaden,
            y_broaden,
            color='black',
            label=f'states: {num}')
    ax.bar(x, y, color='black', width=0.4)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel(r'$ \rm{\Delta \epsilon \ [L \cdot mol^{-1} \cdot cm^{-1}]} $')
    ax.legend()
plt.show()