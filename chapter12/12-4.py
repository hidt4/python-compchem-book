"""ECDスペクトルの描画
"""

from __future__ import annotations
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations


def optfreq(mol: psi4.core.Molecule,
            theory: str = 'b3lyp-d3bj/6-31g(d)') -> list[float, psi4.core.Wavefunction]:
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

    poles = [excited_state['EXCITATION ENERGY']
             for excited_state in excited_states]
    residues = [excited_state['ROTATORY STRENGTH (LEN)']
                for excited_state in excited_states]
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
psi4.set_memory('28GB')
psi4.set_options({'SAVE_JK': True})
psi4.set_output_file('ECD_sample.log')

s_form = psi4.geometry('''
0 1
 C                 -0.45228216    0.14522821    0.00000000
 C                  0.94287784    0.14522821    0.00000000
 C                  1.64041584    1.35297921    0.00000000
 C                  0.94276184    2.56148821   -0.00119900
 C                 -0.45206316    2.56141021   -0.00167800
 C                 -1.14966416    1.35320421   -0.00068200
 H                 -1.00204116   -0.80708879    0.00045000
 H                  1.49238584   -0.80728479    0.00131500
 H                  1.49296184    3.51363121   -0.00125800
 H                 -1.00218516    3.51369121   -0.00263100
 H                 -2.24926816    1.35338721   -0.00086200
 C                  3.18041558    1.35309125    0.00088786
 H                  3.53672179    2.22050769    0.51618835
 C                  3.69458408    1.37096069   -1.45063253
 H                  4.76451996    1.35925545   -1.45018335
 H                  3.32868914    0.50953646   -1.96925494
 H                  3.34778992    2.25650596   -1.94098279
 O                  3.65678760    0.17734468    0.66088007
 H                  4.51713271    0.35378051    1.04853545
''')

r_form = psi4.geometry('''
0 1
 C                 -0.45228216    0.14522821    0.00000000
 C                  0.94287784    0.14522821    0.00000000
 C                  1.64041584    1.35297921    0.00000000
 C                  0.94276184    2.56148821    0.00119900
 C                 -0.45206316    2.56141021    0.00167800
 C                 -1.14966416    1.35320421    0.00068200
 H                 -1.00204116   -0.80708879   -0.00045000
 H                  1.49238584   -0.80728479   -0.00131500
 H                  1.49296184    3.51363121    0.00125800
 H                 -1.00218516    3.51369121    0.00263100
 H                 -2.24926816    1.35338721    0.00086200
 C                  3.18041558    1.35309125   -0.00088786
 H                  3.53672179    2.22050769   -0.51618835
 C                  3.69458408    1.37096069    1.45063253
 H                  4.76451996    1.35925545    1.45018335
 H                  3.32868914    0.50953646    1.96925494
 H                  3.34778992    2.25650596    1.94098279
 O                  3.65678760    0.17734468   -0.66088007
 H                  4.51713271    0.35378051   -1.04853545
''')

# 構造最適化と振動数計算
_, wfn_s = optfreq(s_form, 'cam-b3lyp/cc-pvdz')
print(s_form.save_string_xyz())
print(get_freqs(wfn_s))

_, wfn_r = optfreq(r_form, 'cam-b3lyp/cc-pvdz')
print(r_form.save_string_xyz())
print(get_freqs(wfn_r))

# TD-DFT計算
_, wfn_tds = psi4.energy('cam-b3lyp/aug-cc-pvdz',
                         molecule=s_form,
                         return_wfn=True)
excited_states_S = tdscf_excitations(wfn_tds, states=6)

_, wfn_tdr = psi4.energy('cam-b3lyp/aug-cc-pvdz',
                         molecule=r_form,
                         return_wfn=True)
excited_states_R = tdscf_excitations(wfn_tdr, states=6)

# ECDスペクトル
x_broaden_S, y_broaden_S, x_S, y_S = get_ECD_spectrum(excited_states_S)
x_broaden_R, y_broaden_R, x_R, y_R = get_ECD_spectrum(excited_states_R)

# 描画
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_broaden_S,
        y_broaden_S,
        color='black',
        label='S-form')
ax.bar(x_S, y_S, color='black')
ax.plot(x_broaden_R,
        y_broaden_R,
        '--',
        color='gray',
        label='R-form')
ax.bar(x_R, y_R, color='gray')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel(r'$ \rm{\delta \epsilon}$')
ax.legend()
plt.show()