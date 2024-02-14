"""アセトアミド
"""

import numpy as np
import psi4


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


# 計算レベルとログファイル名のリスト
theories = ['hf/sto-3g', 'hf/3-21g', 'mp2/6-31g', 'wb97x-d/aug-cc-pvdz']
logfiles = ['hf_sto3g', 'hf_3-21g', 'mp2_6-31g', 'wb97xd_aug-cc-pvdz']

# アセトアミドの構造定義
acetamide = psi4.geometry('''
0 1
 C                  0.11203320    0.61825725    0.00000000
 O                  1.33489796    0.61545082   -0.10440743
 C                 -0.70916294    1.92103771    0.00000000
 H                 -0.40304404    2.53560774   -0.82066734
 H                 -1.74874842    1.68735489   -0.09774644
 H                 -0.54626743    2.44532874    0.91841383
 N                 -0.67183364   -0.62530452    0.00262330
 H                 -0.42379451   -1.17736195   -0.79343523
 H                 -0.47841040   -1.13949341    0.83820527
 ''')

# 構造最適化と振動数計算の実行
co_amide = []

for theory, logfile in zip(theories, logfiles):
    filename = f'acetamide_{logfile}.log'
    psi4.set_output_file(filename)

    _, wfn = optfreq(acetamide, theory)
    co_amide.append(get_freqs(wfn)[15])

print(co_amide)