# psi4のimport
import psi4

# 計算ログファイルの設定
psi4.set_output_file('h2o_psiapi.log')

# 水分子の構造定義
h2o = psi4.geometry('''
0 1
O       -0.1176269719      0.7387773605      0.0000000000
H        0.8523730281      0.7387773605      0.0000000000
H       -0.4409567836      1.4770262439     -0.5397651517
''')

# エネルギー計算の実行
psi4.energy('hf/sto-3g', molecule=h2o)