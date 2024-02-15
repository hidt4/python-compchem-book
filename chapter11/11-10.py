"""アルゴン二量体のSAPT0計算
"""

from __future__ import annotations
import numpy as np
import pandas as pd
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


psi4.set_num_threads(14)
psi4.set_memory('20GB')
psi4.set_output_file('SAPT0_ar_dimer.log')

elst_energies = []
exch_energies = []
ind_energies = []
disp_energies = []
sapt_energies = []

ar_geom = '''
0 1
 Ar     0   0   0 
--
0 1
 Ar     0   0   {}
'''

for distance in np.arange(1.5, 5.5, 0.1):
    ar_dimer = psi4.geometry(ar_geom.format(distance))
    psi4.energy('sapt0/jun-cc-pvdz', molecule=ar_dimer)
    elst, exch, ind, disp, sapt = get_SAPT_energies()
    elst_energies.append(elst)
    exch_energies.append(exch)
    ind_energies.append(ind)
    disp_energies.append(disp)
    sapt_energies.append(sapt)

