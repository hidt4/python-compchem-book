"""エタンを高精度で構造最適化する
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


# 2つの方法用に同じ初期構造で分子を2つ準備
hcho = psi4.geometry('''
0 1
 C                 -0.70124482    0.39419087    0.00000000
 O                  0.52607218    0.39419087    0.00000000
 H                 -1.29338982    1.33359487    0.00000000
 H                 -1.29338982   -0.54521313    0.00003900
 ''')

hcho2 = hcho.clone()

psi4.set_num_threads(14)
psi4.set_memory('20GB')

# ログファイルの設定
psi4.set_output_file('hcho-low-to-high.log')

# 計算レベルの設定
level1 = 'mp2/cc-pvdz'
level2 = 'mp2/aug-cc-pvqz'

# 低レベルでまず最適化してから高精度で最適化
time_1_1 = record_optimization(hcho, level1)
time_1_2 = record_optimization(hcho, level2)
print(hcho.save_string_xyz())
print(f'level 1: {time_1_1: .2f} sec\tlevel 2: {time_1_2: .2f} sec')
print(f'level 1 + 2: {time_1_1 + time_1_2: .2f} sec')
print('####')

# 直接高精度で最適化
time_2 = record_optimization(hcho2, level2)
print(hcho2.save_string_xyz())
print(f'level 2: {time_2: .2f} sec')
print('####')

# 2つの方法で得た構造の差分を計算
for i in range(hcho.natom()):
    vec1 = hcho.xyz(i)
    vec2 = hcho2.xyz(i)
    dist = vec1.distance(vec2)
    symbol = hcho.symbol(i)
    print(f'{symbol}\t{dist: .3f}')