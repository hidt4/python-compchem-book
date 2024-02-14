"""ホモデスモチック反応
シクロプロパンの歪みエネルギーを推定する
"""

import numpy as np
import psi4

au2kcal = psi4.constants.hartree2kcalmol

psi4.set_num_threads(14)
psi4.set_memory('20GB')
psi4.set_options({'normal_modes_write': True})


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


# シクロプロパン
psi4.set_output_file('cyclopropane.log')

cyclopropane = psi4.geometry('''
0 1
 C                  0.16445241    0.38365967    0.00286974
 C                  1.66556741    0.38365967    0.00286974
 C                  0.91498241    1.68363467    0.00286974
 H                 -0.37214359    0.07378467   -0.91098826
 H                 -0.37214359    0.07378467    0.91672774
 H                  2.20212841    0.07389267   -0.91101826
 H                  2.20212841    0.07389267    0.91675774
 H                  0.91475941    2.30316267    0.91681374
 H                  0.91475941    2.30316267   -0.91107426
 ''')

energy_cpropane, wfn_cpropane = optfreq(cyclopropane, 'hf/cc-pvdz')
print(f'cyclopropane freqs:\n{get_freqs(wfn_cpropane)}')
cyclopropane_H = psi4.variable('ENTHALPY')


# メタン
psi4.set_output_file('methane.log')

methane = psi4.geometry('''
0 1
 C                 -0.40248963    0.33609958    0.00000000
 H                 -0.04583521   -0.67271042    0.00000000
 H                 -0.04581679    0.84049777    0.87365150
 H                 -0.04581679    0.84049777   -0.87365150
 H                 -1.47248963    0.33611276    0.00000000
 ''')

energy_methane, wfn_methane = optfreq(methane, 'hf/cc-pvdz')
print(f'methane freqs:\n{get_freqs(wfn_methane)}')
methane_H = psi4.variable('ENTHALPY')

# エタン
psi4.set_output_file('ethane.log')

ethane = psi4.geometry('''
0 1
 C                 -0.40248963    0.33609958    0.00000000
 H                 -0.04583521   -0.67271042    0.00000000
 H                 -0.04581679    0.84049777   -0.87365150
 H                 -1.47248963    0.33611276    0.00000000
 C                  0.11085259    1.06205585    1.25740497
 H                  1.18085259    1.06202576    1.25741439
 H                 -0.24578624    2.07087137    1.25739542
 H                 -0.24583592    0.55766838    2.13105626
 ''')

energy_ethane, wfn_ethane = optfreq(ethane, 'hf/cc-pvdz')
print(f'ethane freqs:\n{get_freqs(wfn_ethane)}')
ethane_H = psi4.variable('ENTHALPY')

# アイソデスミック反応による推定
isodesmic_h = au2kcal * (3 * ethane_H - (cyclopropane_H + 3 * methane_H))
print(f'estimated strain energy using isodesmic reaction: {-isodesmic_h: .2f} kcal/mol')