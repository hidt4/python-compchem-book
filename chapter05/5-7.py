"""アルゴン二量体のエネルギー最小値を求める
"""

import psi4


def get_minima(theory: str,
               geom: str,
               distance_min: float,
               distance_max: float) -> float:
    """
    指定レベルで探索範囲内でエネルギー最小値を与える原子間距離を求める。

    Args:
        theory: 計算レベル
        geom: 計算対象分子の構造
        distance_min: 探索する最短原子間距離
        distance_max: 探索する最長原子間距離

    Returns:
        float: エネルギー最小値を与える原子間距離
    """
    from scipy import optimize

    level = f'{theory}/cc-pvtz'


    def object_function(length: float) -> float:
        """
        最小化する目的関数。指定原子間距離でエネルギー計算を行う。

        引数：
            length: 原子間距離
        戻り値：
            float: エネルギー値
        """
        mol = psi4.geometry(geom.format(length))
        energy = psi4.energy(level, molecule=mol)

        return energy


    distance = optimize.fminbound(object_function,
                                  x1=distance_min,
                                  x2=distance_max)

    return distance


Ar2_geom = '''
0 1
Ar  0   0   0
Ar  0   0   {}
'''

methods = ['mp2', 'b3lyp', 'b3lyp-d3bj', 'wb97x-d']

for method in methods:
    dist = get_minima(method, Ar2_geom, 1.5, 5.5)
    print(f'{method.upper()}/cc-pVTZ: {dist: .2f} angstrom')