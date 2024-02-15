"""リチウムカチオンとベンゼンの相互作用をSAPTで計算
"""

from __future__ import annotations
import numpy as np
import psi4


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


def get_SAPT_energies() -> list[float, float, float, float, float]:
    """
    psi4.variableからSAPT計算の結果を取得してkcal/mol単位で返す

    Returns:
        float: 静電力による寄与
        float: 交換反発力による寄与
        float: 誘起力による寄与
        float: 分散力による寄与
        float: 相互作用エネルギー合計

    """
    au2kcal = psi4.constants.hartree2kcalmol
    elst_energy = psi4.variable('SAPT ELST ENERGY') * au2kcal
    exch_energy = psi4.variable('SAPT EXCH ENERGY') * au2kcal
    ind_energy = psi4.variable('SAPT IND ENERGY') * au2kcal
    disp_energy = psi4.variable('SAPT DISP ENERGY') * au2kcal
    sapt_energy = psi4.variable('SAPT TOTAL ENERGY') * au2kcal

    return [elst_energy, exch_energy, ind_energy, disp_energy, sapt_energy]


def print_SAPT_results(sapt_level: str = 'SAPT0'):
    """
    SAPT計算の結果を出力する

    """
    elst, exch, ind, disp, sapt = get_SAPT_energies()

    print(f'{sapt_level}')
    print(f'SAPT energy: {sapt: .2f}')
    print('--------------------')
    print(f'Electrostatic term: {elst: .2f}')
    print(f'Exchange term: {exch: .2f}')
    print(f'Induction term: {ind: .2f}')
    print(f'Dispersion term: {disp: .2f}')


psi4.set_num_threads(14)
psi4.set_memory('28GB')
psi4.set_output_file('SAPT0_Li-benzene.log')
theory = 'mp2/cc-pvdz'

# リチウムカチオンとベンゼンの複合体の構造定義
Li_benzene = psi4.geometry('''
0 1
 C                 -0.42021926   -1.46772227    0.00000000
 C                  0.97494074   -1.46772227    0.00000000
 C                  1.67247874   -0.25997127    0.00000000
 C                  0.97482474    0.94853773   -0.00119900
 C                 -0.42000026    0.94845973   -0.00167800
 C                 -1.11760126   -0.25974627   -0.00068200
 H                 -0.96997826   -2.42003927    0.00045000
 H                  1.52444874   -2.42023527    0.00131500
 H                  2.77215874   -0.25989127    0.00063400
 H                  1.52502474    1.90068073   -0.00125800
 H                 -0.97012226    1.90074073   -0.00263100
 H                 -2.21720526   -0.25956327   -0.00086200
 --
 1 1
 Li                 0.35180857   -0.32725870    3.23341003
''')

# 構造最適化と振動数計算
_, wfn_Li = optfreq(Li_benzene, theory)
print(get_freqs(wfn_Li))

# SAPT計算
psi4.set_options({'SCF_TYPE': 'DF',
                  'FREEZE_CORE': True,
                  'GUESS': 'SAD'})
psi4.energy('sapt0/jun-cc-pvdz', molecule=Li_benzene)
print_SAPT_results()