"""按元素守恒配平化学方程式。

把每种物质在各元素上的原子数排成矩阵（反应物取正、产物取负），
用 SymPy 求零空间。零空间维数恰为 1 时比例唯一，取向量中的
有理数通分为严格正整数并约去公因子。
"""

from __future__ import annotations

from functools import reduce
from math import gcd, lcm

import sympy as sp


class BalanceError(Exception):
    """配平失败，code 为机器可读的错误类别。"""

    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code


def balance_equation(
    reactants: list[dict[str, int]], products: list[dict[str, int]]
) -> tuple[list[int], list[int]]:
    """返回 (反应物系数, 产物系数)，均为严格正整数且整体互质。"""
    species = list(reactants) + list(products)
    n_reactants = len(reactants)

    elements: list[str] = []
    seen: set[str] = set()
    for counts in species:
        for symbol in counts:
            if symbol not in seen:
                seen.add(symbol)
                elements.append(symbol)

    matrix = sp.zeros(len(elements), len(species))
    for col, counts in enumerate(species):
        sign = 1 if col < n_reactants else -1
        for row, symbol in enumerate(elements):
            amount = counts.get(symbol)
            if amount:
                matrix[row, col] = sign * amount

    nullspace = matrix.nullspace()
    if not nullspace:
        raise BalanceError(
            "no_solution", "no non-trivial balance exists for these species"
        )
    if len(nullspace) > 1:
        raise BalanceError(
            "non_unique_solution",
            "the coefficient ratio is not unique for these species",
        )

    vector = list(nullspace[0])
    if any(v == 0 for v in vector):
        raise BalanceError(
            "no_positive_solution",
            "no strictly positive coefficient assignment balances the equation",
        )
    if all(v < 0 for v in vector):
        vector = [-v for v in vector]
    if not all(v > 0 for v in vector):
        raise BalanceError(
            "no_positive_solution",
            "species cannot all be balanced on their given sides",
        )

    common_denominator = reduce(lcm, (int(v.q) for v in vector), 1)
    coefficients = [int(v * common_denominator) for v in vector]
    common_factor = reduce(gcd, coefficients)
    coefficients = [c // common_factor for c in coefficients]

    return coefficients[:n_reactants], coefficients[n_reactants:]
