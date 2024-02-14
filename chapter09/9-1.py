"""手法ごとのWavefunctionオブジェクトの違い
"""

import psi4


def calc_energy(mol: psi4.core.Molecule,
                theory: str = 'hf/sto-3g') -> psi4.core.Wavefunction:

    """
    エネルギー計算を行い，Wavefunctionオブジェクトを返す
    Args:
        mol: 計算を行う分子
        theory: 計算レベル

    Returns: エネルギー計算で得られるWavefunctionオブジェクト

    """
    _, wfn = psi4.energy(theory, molecule=mol, return_wfn=True)

    return wfn


methods = ['hf', 'mp2', 'b3lyp']

# 水分子の構造定義
h2o = psi4.geometry('''
0 1
 O                 -0.16182573   -0.32780082    0.00000000
 H                  0.79817427   -0.32780082    0.00000000
 H                 -0.48228031    0.57713501    0.00000000
 ''')

# 各手法でエネルギー計算
for method in methods:
    psi4.set_output_file(f'h2o_{method}.log')
    theory = f'{method}/sto-3g'
    wfn = calc_energy(h2o, theory=theory)

    # 各Wavefunctionオブジェクトのvariablesを出力
    print(f'================   {method.upper()}   ================')
    for key in wfn.variables():
        if type(wfn.variable(key)) == float:
            print(f'{key}\t{wfn.variable(key): .4f}')
        else:
            print(key, '\t', wfn.variable(key))
    print()