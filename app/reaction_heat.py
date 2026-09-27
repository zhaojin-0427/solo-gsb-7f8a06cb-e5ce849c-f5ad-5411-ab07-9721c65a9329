"""由已知反应热推导目标反应热（Hess 定律）。

每个反应表示为物种净系数向量（消耗为负、生成为正）。目标反应若能写成
已知反应的线性组合，则其反应热等于各已知反应热按同一权重的组合
（Hess 定律）。全部运算使用精确有理数（:class:`fractions.Fraction` 与
SymPy Rational），不存在浮点误差。

判定顺序：

1. 若已知反应之间存在「净物质变化为零但热量非零」的线性组合，则数据
   自相矛盾（违反 Hess 定律），任何目标热量都不唯一，报
   ``inconsistent_known_reactions`` 并给出一组见证权重与其热量；
2. 否则若目标反应不在已知反应的张成空间内，报 ``target_not_derivable``；
3. 否则返回一组权重（允许负数与零；不唯一时任取一组）与目标热量。
"""

from __future__ import annotations

import re
from decimal import Decimal
from fractions import Fraction
from typing import NamedTuple

import sympy as sp

# 可带正负号、不含指数、整数段无前导零的有限十进制字符串
_HEAT_PATTERN = re.compile(r"^[+-]?(?:0|[1-9]\d*)(?:\.\d+)?$")


class HeatError(ValueError):
    """heat_kj 不是合法的有限十进制字符串。"""


class ReactionHeatError(Exception):
    """反应热推导失败，code 为机器可读的错误类别。"""

    def __init__(self, code: str, message: str, extra: dict | None = None):
        super().__init__(message)
        self.code = code
        self.extra = extra or {}


class ReactionHeatResult(NamedTuple):
    weights: list[Fraction]
    heat_kj: Fraction


def parse_heat(text: str) -> Fraction:
    """把 kJ 热量字符串解析为精确分数，要求有限纯十进制（可带正负号）。"""
    if not isinstance(text, str) or not _HEAT_PATTERN.fullmatch(text):
        raise HeatError(
            f"heat {text!r} is not a finite plain decimal string "
            "(optional sign, digits with an optional decimal point, no exponent)"
        )
    return Fraction(Decimal(text))


def format_rational(value: Fraction) -> str:
    """把分数格式化为分母为正的最简 p/q 字符串（整数也写作 p/1）。"""
    return f"{value.numerator}/{value.denominator}"


def solve_reaction_heat(
    known_nets: list[dict[str, int]],
    heats_kj: list[Fraction],
    target_net: dict[str, int],
) -> ReactionHeatResult:
    """求一组权重 w 使 Σ w_i·known_i = target，并给出 Σ w_i·heat_i。

    known_nets 与 target_net 为 {物种化学式原文: 净系数}，heats_kj 与
    known_nets 等长、顺序一致。物种以化学式原文为键。
    """
    species: list[str] = []
    seen: set[str] = set()
    for net in list(known_nets) + [target_net]:
        for name in net:
            if name not in seen:
                seen.add(name)
                species.append(name)

    n_known = len(known_nets)
    matrix = sp.Matrix(
        len(species),
        n_known,
        lambda row, col: known_nets[col].get(species[row], 0),
    )

    # 1) 一致性：零空间（净物质变化为零的组合）中的热量必须全部为零
    for basis in matrix.nullspace():
        cycle_weights = [_to_fraction(component) for component in basis]
        cycle_heat = sum(
            weight * heat for weight, heat in zip(cycle_weights, heats_kj)
        )
        if cycle_heat != 0:
            raise ReactionHeatError(
                "inconsistent_known_reactions",
                "known reactions admit a zero-net combination with non-zero "
                "heat, so the given reaction heats are inconsistent",
                {
                    "weights": [format_rational(w) for w in cycle_weights],
                    "heat_kj": format_rational(cycle_heat),
                },
            )

    # 2) 可推导性：精确求解 K·w = target
    target_vector = sp.Matrix([target_net.get(name, 0) for name in species])
    try:
        solution, params = matrix.gauss_jordan_solve(target_vector)
    except ValueError:
        raise ReactionHeatError(
            "target_not_derivable",
            "target reaction is not a linear combination of the known reactions",
        )
    # 解不唯一时把自由变量取 0，得到一组特解（权重不唯一时任一均可）
    particular = solution.xreplace({param: 0 for param in params})
    weights = [_to_fraction(component) for component in particular]

    heat = sum(
        weight * known_heat for weight, known_heat in zip(weights, heats_kj)
    )
    return ReactionHeatResult(weights=weights, heat_kj=heat)


def _to_fraction(value: sp.Expr) -> Fraction:
    """把 SymPy 有理数转成 fractions.Fraction（自动约分、分母为正）。"""
    rational = sp.Rational(value)
    return Fraction(int(rational.p), int(rational.q))
