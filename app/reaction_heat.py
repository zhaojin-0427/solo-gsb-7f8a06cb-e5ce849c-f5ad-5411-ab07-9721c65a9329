"""Hess 定律反应热推导。

已知反应以 net 映射给出（化学式原文 → 净系数：负整数 = 消耗，正整数 = 生成），
连同反应热（kJ，有限十进制字符串）。把每个已知反应看作物种空间中的向量，
用 :class:`fractions.Fraction` 精确有理数做高斯消元：

* 若已知反应的某个线性组合净物质变化为零而热量非零，则数据违反 Hess 定律
  （同一「零反应」对应两个热量），返回该组合的权重与热量；
* 否则若目标反应不在已知反应的张成空间内，则不可推导；
* 可推导时给出一组权重（允许负数与零，不唯一时任取其一）与目标反应热。
"""

from __future__ import annotations

import re
from fractions import Fraction

from .chemistry import FormulaError, parse_formula

MAX_COEFFICIENT = 99  # net 系数绝对值上限
MAX_SPECIES = 12      # 全请求不同物种总数上限

# 可带正负号、不含指数的有限十进制；整数段不允许前导零（与 parse_mass 一致）
_HEAT_PATTERN = re.compile(r"^[+-]?(?:0|[1-9]\d*)(?:\.\d+)?$")


class HeatError(ValueError):
    """heat_kj 不是可带符号的有限十进制字符串。"""


class ReactionHeatError(Exception):
    """反应热语义错误，code 为机器可读类别，extra 携带定位/组合信息。"""

    def __init__(self, code: str, message: str, extra: dict | None = None):
        super().__init__(message)
        self.code = code
        self.extra = extra or {}


def parse_heat(text: str) -> Fraction:
    """把 heat_kj 字符串解析为精确分数；拒绝指数、非有限值与前导零。"""
    if not isinstance(text, str) or not _HEAT_PATTERN.fullmatch(text):
        raise HeatError(
            f"heat {text!r} is not a finite decimal string "
            "(optional sign, digits with an optional decimal point, no exponent)"
        )
    return Fraction(text)


def format_rational(value: Fraction) -> str:
    """格式化为分母为正的最简 p/q 字符串（Fraction 本身即保证）。"""
    return f"{value.numerator}/{value.denominator}"


def validate_net(net: dict[str, int]) -> None:
    """校验单个 net 映射：系数范围、化学式可解析、正负两侧非空、元素守恒。

    错误抛 :class:`ReactionHeatError`，字段路径由调用方补充。
    """
    for species, coefficient in net.items():
        if coefficient == 0 or abs(coefficient) > MAX_COEFFICIENT:
            raise ReactionHeatError(
                "invalid_coefficient",
                f"coefficient for species {species!r} must be a non-zero "
                f"integer with absolute value 1-{MAX_COEFFICIENT}",
                {"species": species, "value": coefficient},
            )

    parsed: dict[str, dict[str, int]] = {}
    for species in net:
        try:
            parsed[species] = parse_formula(species)
        except FormulaError as exc:
            raise ReactionHeatError(
                "invalid_formula",
                f"species {species!r} is not a valid formula: {exc}",
                {"species": species},
            ) from exc

    has_consumed = any(coefficient < 0 for coefficient in net.values())
    has_produced = any(coefficient > 0 for coefficient in net.values())
    if not has_consumed or not has_produced:
        missing = "consumed" if not has_consumed else "produced"
        raise ReactionHeatError(
            "empty_side",
            "net must contain at least one consumed (negative coefficient) "
            f"and one produced (positive coefficient) species; "
            f"the {missing} side is empty",
            {"missing_side": missing},
        )

    imbalance: dict[str, int] = {}
    for species, coefficient in net.items():
        for element, count in parsed[species].items():
            imbalance[element] = imbalance.get(element, 0) + coefficient * count
    nonzero = {element: n for element, n in imbalance.items() if n != 0}
    if nonzero:
        raise ReactionHeatError(
            "unbalanced_reaction",
            "net reaction does not conserve elements (net atom change: "
            + ", ".join(f"{element} {n:+d}" for element, n in nonzero.items())
            + ")",
            {"imbalance": nonzero},
        )


def solve_reaction_heat(
    known_nets: list[dict[str, int]],
    heats_kj: list[Fraction],
    target_net: dict[str, int],
) -> tuple[list[Fraction], Fraction]:
    """求 target = Σ weights[i]·known[i] 的一组权重与目标反应热（精确分数）。

    已知数据自相矛盾时抛 ``inconsistent_known``（附组合权重与热量），
    目标不可由已知反应组合时抛 ``not_derivable``。
    """
    # 物种顺序：按 known 顺序首次出现，最后是 target 独有的物种
    species: list[str] = []
    seen: set[str] = set()
    for net in [*known_nets, target_net]:
        for name in net:
            if name not in seen:
                seen.add(name)
                species.append(name)

    n_known = len(known_nets)
    # 每行一个物种、每列一个已知反应的系数矩阵
    rows = [
        [Fraction(net.get(name, 0)) for net in known_nets] for name in species
    ]
    target_column = [Fraction(target_net.get(name, 0)) for name in species]

    # 一致性：系数矩阵零空间中的组合净物质变化为零，其热量必须也为零
    for combination in _nullspace(rows, n_known):
        heat = sum(
            (w * h for w, h in zip(combination, heats_kj)), Fraction(0)
        )
        if heat != 0:
            raise ReactionHeatError(
                "inconsistent_known",
                "a linear combination of the known reactions has zero net "
                "species change but non-zero heat; the known data is "
                "self-contradictory",
                {
                    "weights": [format_rational(w) for w in combination],
                    "heat_kj": format_rational(heat),
                },
            )

    weights = _solve(rows, n_known, target_column)
    if weights is None:
        raise ReactionHeatError(
            "not_derivable",
            "target reaction cannot be expressed as a linear combination "
            "of the known reactions",
        )

    target_heat = sum(
        (w * h for w, h in zip(weights, heats_kj)), Fraction(0)
    )
    return weights, target_heat


def _rref(
    rows: list[list[Fraction]], n_cols: int
) -> tuple[list[list[Fraction]], list[int]]:
    """精确行最简形，返回 (化简后的行, 主元列下标)。"""
    matrix = [list(row) for row in rows]
    n_rows = len(matrix)
    pivots: list[int] = []
    r = 0
    for c in range(n_cols):
        if r >= n_rows:
            break
        pivot = next(
            (i for i in range(r, n_rows) if matrix[i][c] != 0), None
        )
        if pivot is None:
            continue
        matrix[r], matrix[pivot] = matrix[pivot], matrix[r]
        scale = matrix[r][c]
        matrix[r] = [x / scale for x in matrix[r]]
        for i in range(n_rows):
            if i != r and matrix[i][c] != 0:
                factor = matrix[i][c]
                matrix[i] = [a - factor * b for a, b in zip(matrix[i], matrix[r])]
        pivots.append(c)
        r += 1
    return matrix, pivots


def _nullspace(rows: list[list[Fraction]], n_cols: int) -> list[list[Fraction]]:
    """矩阵（rows，n_cols 列）零空间的一组基，自由列依次取 1。"""
    reduced, pivots = _rref(rows, n_cols)
    free_columns = [c for c in range(n_cols) if c not in set(pivots)]
    basis: list[list[Fraction]] = []
    for free in free_columns:
        vector = [Fraction(0)] * n_cols
        vector[free] = Fraction(1)
        for row_index, pivot_col in enumerate(pivots):
            vector[pivot_col] = -reduced[row_index][free]
        basis.append(vector)
    return basis


def _solve(
    rows: list[list[Fraction]], n_cols: int, rhs: list[Fraction]
) -> list[Fraction] | None:
    """解 A·x = rhs（精确），无解返回 None；有自由变量时取零得一组特解。"""
    augmented = [list(row) + [b] for row, b in zip(rows, rhs)]
    reduced, pivots = _rref(augmented, n_cols + 1)
    if pivots and pivots[-1] == n_cols:
        # 增广列成为主元列：存在 [0 ... 0 | 非零] 行，方程组不相容
        return None
    particular = [Fraction(0)] * n_cols
    for row_index, pivot_col in enumerate(pivots):
        particular[pivot_col] = reduced[row_index][n_cols]
    return particular
