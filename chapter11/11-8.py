"""三フッ化ホウ素とアンモニアの相互作用をSAPT0計算
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


au2kcal = psi4.constants.hartree2kcalmol
psi4.set_num_threads(14)
psi4.set_memory('32GB')
psi4.set_output_file('SAPT0_BF3_NH3.log')

bf3_nh3 = psi4.geometry('''
0 1
 B    0.125118224796   -0.000000000409    0.000000000536
 F    0.455918039271   -0.662711022072   -1.147849157116
 F    0.455918029138   -0.662711022072    1.147849160853
 F    0.455918039220    1.325422040779    0.000000002202
--
0 1
 N   -1.545025507582    0.000000003860   -0.000000006706
 H   -1.894344094062    0.474567048776    0.821974219832
 H   -1.894344086834    0.474567045504   -0.821974238181
 H   -1.894344094077   -0.949134080011   -0.000000006271
''')

# SAPT計算
psi4.set_options({'SCF_TYPE': 'DF',
                  'FREEZE_CORE': True,
                  'GUESS': 'SAD'})
psi4.energy('sapt0/jun-cc-pvdz', molecule=bf3_nh3)
print_SAPT_results()

# SAPT-CT計算
psi4.set_output_file('SAPT0-CT_BF3_NH3.log')
psi4.energy('sapt0-ct/jun-cc-pvdz', molecule=bf3_nh3)
ct = psi4.variable('SAPT CT ENERGY')
print(f'Charge transfer: {au2kcal * ct: .2f}')
