"""メチルシクロヘキサンの自由エネルギー差
"""

import numpy as np
import psi4

au2kcal = psi4.constants.hartree2kcalmol

psi4.set_num_threads(14)
psi4.set_memory('20GB')


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


# エクアトリアル置換
psi4.set_output_file('equatorial.log')

equatorial = psi4.geometry('''
0 1
 C                 -2.15838901   -1.67295376    0.00000000
 C                 -0.64328301   -1.67295376    0.00000000
 C                 -0.09135201   -0.26187576    0.00000000
 C                 -0.64101501    0.54266124    1.16066100
 C                 -2.15614001    0.54332224    1.16017200
 C                 -2.70894001   -0.86729876    1.15887600
 H                  1.02724699   -0.29586676    0.06271400
 H                 -0.27073401   -2.21858776    0.90656200
 H                 -0.26798901   -2.22281276   -0.90191000
 H                 -2.53107701   -1.23992476   -0.96538500
 H                 -2.53398601   -2.72717376    0.06350200
 H                 -0.26848001    0.10781624    2.12527200
 H                 -0.26499201    1.59679824    1.09866600
 H                 -2.53146001    1.09275024    2.06228600
 H                 -2.52799501    1.08992324    0.25384900
 H                 -2.44293601   -1.37302076    2.12419100
 H                 -0.35997001    0.24375724   -0.96454600
 C                 -4.24559913   -0.81913202    1.06966430
 H                 -4.63527674   -1.81550142    1.08696717
 H                 -4.63411695   -0.26971208    1.90158475
 H                 -4.53508198   -0.33871606    0.15845621
 ''')

_, equatorial_wfn = optfreq(equatorial, 'hf/cc-pvdz')
print(f'equatorial freqs:\n{get_freqs(equatorial_wfn)}')
equatorial_G = psi4.variable('GIBBS FREE ENERGY')
equatorial_Gcorr = psi4.variable('GIBBS FREE ENERGY CORRECTION')


# アキシアル置換
psi4.set_output_file('axial.log')

axial = psi4.geometry('''
0 1
 C                 -2.15838901   -1.67295376    0.00000000
 C                 -0.64328301   -1.67295376    0.00000000
 C                 -0.09135201   -0.26187576    0.00000000
 C                 -0.64101501    0.54266124    1.16066100
 C                 -2.15614001    0.54332224    1.16017200
 C                 -2.70894001   -0.86729876    1.15887600
 H                  1.02724699   -0.29586676    0.06271400
 H                 -0.27073401   -2.21858776    0.90656200
 H                 -0.26798901   -2.22281276   -0.90191000
 H                 -2.53107701   -1.23992476   -0.96538500
 H                 -2.53398601   -2.72717376    0.06350200
 H                 -0.26848001    0.10781624    2.12527200
 H                 -0.26499201    1.59679824    1.09866600
 H                 -2.53146001    1.09275024    2.06228600
 H                 -2.52799501    1.08992324    0.25384900
 H                 -0.35997001    0.24375724   -0.96454600
 C                 -2.34375843   -1.56157538    2.48410141
 H                 -2.74805066   -1.00326314    3.30247494
 H                 -2.75034420   -2.55128438    2.49210909
 H                 -1.27915043   -1.61256543    2.57849329
 H                 -3.77661875   -0.83383226    1.09689125
 ''')

_, axial_wfn = optfreq(axial, 'hf/cc-pvdz')
print(f'axial freqs:\n{get_freqs(equatorial_wfn)}')
axial_G = psi4.variable('GIBBS FREE ENERGY')
axial_Gcorr = psi4.variable('GIBBS FREE ENERGY CORRECTION')


# 自由エネルギー差を計算
diff = au2kcal * (axial_G - equatorial_G)
print(f'free energy difference: {diff: .2f} kcal/mol')