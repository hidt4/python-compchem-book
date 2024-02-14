"""ゴーシュ配座の1,2-ジクロロエタンの構造最適化
"""

import psi4

psi4.set_num_threads(14)
psi4.set_memory('20GB')

# ログファイルの設定
psi4.set_output_file('DCE_gauche.log')

# gauche配座の分子を定義
dce_gauche = psi4.geometry('''
0 1
C       -4.7587744478     -0.1541943645     -0.0559642548
C       -3.2346420969     -0.1669054169      0.0013710608
H       -5.1443980180     -1.1552590965      0.2305585920
H       -5.0986095260      0.0746426192     -1.0881859644
Cl      -5.4379603530      1.0500010949      1.0719697870
H       -2.8490185267     -0.9530856920     -0.6813781716
H       -2.8948070187     -0.3957485023      1.0335914176
Cl      -2.5554561917      1.4011538588     -0.5119631531
''')

# 構造最適化
psi4.optimize('hf/sto-3g', molecule=dce_gauche)

# 最適化構造を出力
print(dce_gauche.save_string_xyz())

# 最適化構造をXYZファイルとして保存
dce_gauche.save_xyz_file('dce_gauche.xyz', True)