"""FULL_HESS_EVERYオプション
"""

import psi4


def record_optimization(mol: psi4.core.Molecule, theory: str) -> float:
    """
    分子を指定計算レベルで最適化する時間を計測する
    Args:
        mol: 目的とする分子
        theory: 計算レベル

    Returns:
        float: 最適化にかかった時間

    """
    import time
    start = time.time()
    psi4.optimize(theory, molecule=mol)
    end = time.time()

    return end - start


psi4.set_num_threads(14)
psi4.set_memory('20GB')

# ログファイルの設定
psi4.set_output_file('full_hess_every.log')

# 異なる値用に同じ初期構造で分子を3つ準備
hcho = psi4.geometry('''
0 1
 C                 -0.70124482    0.39419087    0.00000000
 O                  0.52607218    0.39419087    0.00000000
 H                 -1.29338982    1.33359487    0.00000000
 H                 -1.29338982   -0.54521313    0.00003900
 ''')

hcho1 = hcho.clone()
hcho2 = hcho.clone()

# 最適化計算
for (i, mol) in zip([-1, 0, 1], [hcho, hcho1, hcho2]):
    psi4.core.clean_options()
    psi4.set_options({'FULL_HESS_EVERY': i})
    time = record_optimization(mol=mol, theory='hf/3-21g')
    print(f'FULL_HESS_EVERY = {i}: {time: .2f} sec')
