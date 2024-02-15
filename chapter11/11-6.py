"""水二量体のSAPT0計算
"""

import psi4


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


def record_elapsed_time(theory: str,
                        molecule: psi4.core.Molecule) -> float:
    """
    指定メソッドでエネルギー計算を行い，実行にかかった時間を返す
    Args:
        theory: 計算レベル
        molecule: 計算対象の分子

    Returns:
        float: 計算にかかった時間（秒）

    """
    import time
    start = time.time()
    psi4.energy(theory, molecule=molecule)
    end = time.time()

    return end - start


psi4.set_num_threads(14)
psi4.set_memory('20GB')

# 水二量体の構造定義（最適化済み）
water_dimer = psi4.geometry('''
0 1
 O    1.412128481591   -0.092045390825   -0.029906791666
 H    1.719638124602    0.616275423632   -0.598773926023
 H    1.722516776196    0.155785589195    0.843086003713
 --
0 1
 H   -0.525178885277   -0.028828516343   -0.009223430301
 O   -1.480679643058    0.092866715086    0.030210284186
 H   -1.829019336855   -0.756267508525   -0.239905293928
''')

# SAPT0/jun-cc-pvdzレベル
psi4.set_output_file('SAPT0_water_dimer.log')
sapt0_time = record_elapsed_time(theory='sapt0/jun-cc-pvdz',
                                 molecule=water_dimer)

print(f'SAPT0: {sapt0_time: .2f} sec')
print_SAPT_results('SAPT0')
print('######')


# sSAPT0/jun-cc-pvdzレベル
psi4.set_output_file('sSAPT_water_dimer.log')
ssapt0_time = record_elapsed_time(theory='ssapt0/jun-cc-pvdz',
                                  molecule=water_dimer)

print(f'sSAPT0: {ssapt0_time: .2f} sec')
print_SAPT_results('sSAPT0')
print('######')

# SAPT2+/aug-cc-pvdzレベル
psi4.set_output_file('SAPT2+_water_dimer.log')
sapt2_time = record_elapsed_time(theory='sapt2+/aug-cc-pvdz',
                                 molecule=water_dimer)

print(f'SAPT2+: {sapt2_time: .2f} sec')
print_SAPT_results('SAPT2+')

