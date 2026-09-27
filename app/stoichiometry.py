"""化学计量计算：摩尔质量、限量试剂、目标产物质量与反应物剩余质量。

原子量与反应物质量均为有限小数 :class:`decimal.Decimal`，摩尔质量用
高精度 Decimal 上下文计算；由于这些量都是有限小数，将其换算成
:class:`fractions.Fraction` 后，物质的量之比（n/系数）与最终质量完全
精确，不会因浮点或有限位除法把真正的并列限量试剂判出高低。仅在最终
输出时把分数转回 Decimal，按 ROUND_HALF_UP 保留六位小数。
"""

from __future__ import annotations

import re
from decimal import Decimal, ROUND_HALF_UP, localcontext
from fractions import Fraction
from typing import NamedTuple

# 本服务支持的元素及其相对原子质量
ATOMIC_MASSES: dict[str, Decimal] = {
    "H": Decimal("1.008"),
    "C": Decimal("12.011"),
    "N": Decimal("14.007"),
    "O": Decimal("15.999"),
    "S": Decimal("32.06"),
    "Cl": Decimal("35.45"),
    "Na": Decimal("22.99"),
    "K": Decimal("39.098"),
    "Ca": Decimal("40.078"),
    "P": Decimal("30.974"),
}
SUPPORTED_ELEMENTS: list[str] = sorted(ATOMIC_MASSES)

MIN_MASS_G = Decimal("0.001")
MAX_MASS_G = Decimal("1000000000")  # 1e9

# 仅接受不带指数、不带符号、整数段无前导零的十进制字符串
_DECIMAL_PATTERN = re.compile(r"^(?:0|[1-9]\d*)(?:\.\d+)?$")
_SIX_PLACES = Decimal("0.000001")
_NINE_PLACES = Decimal("0.000000001")
_CALC_PRECISION = 100


class UnsupportedElementError(ValueError):
    """化学式含有本服务未提供原子量的元素。"""

    def __init__(self, element: str):
        super().__init__(f"element {element!r} has no atomic weight configured")
        self.element = element


class MassError(ValueError):
    """质量字符串不是取值范围内的合法十进制数。"""


class StoichiometryError(Exception):
    """化学计量语义错误（如并列限量试剂），code 为机器可读类别。"""

    def __init__(self, code: str, message: str, extra: dict | None = None):
        super().__init__(message)
        self.code = code
        self.extra = extra or {}


def parse_mass(text: str) -> Decimal:
    """把克质量字符串解析为 Decimal，要求 0.001～1e9 的纯十进制数。"""
    if not isinstance(text, str) or not _DECIMAL_PATTERN.fullmatch(text):
        raise MassError(
            f"mass {text!r} is not a plain decimal string "
            "(digits with an optional decimal point, no exponent/sign)"
        )
    value = Decimal(text)
    if value < MIN_MASS_G or value > MAX_MASS_G:
        raise MassError(
            f"mass {text!r} is out of range [{MIN_MASS_G}, {MAX_MASS_G}] grams"
        )
    return value


def molar_mass(counts: dict[str, int]) -> Decimal:
    """按 {元素: 原子个数} 用 Decimal 求摩尔质量（g/mol）。"""
    if not counts:
        raise ValueError("a species must contain at least one element")
    with localcontext() as ctx:
        ctx.prec = _CALC_PRECISION
        total = Decimal(0)
        for element, n in counts.items():
            atomic = ATOMIC_MASSES.get(element)
            if atomic is None:
                # 正常情况下 API 层已先行校验，这里做防御
                raise UnsupportedElementError(element)
            total += atomic * n
        return +total  # 显式按当前上下文取整规范


class StoichiometryResult(NamedTuple):
    limiting_reactant_index: int
    target_product_mass_g: Fraction
    remaining_reactant_masses_g: list[Fraction]


def solve_stoichiometry(
    coefficients: list[int],
    masses_g: list[Decimal],
    reactant_molar_masses: list[Decimal],
    product_molar_masses: list[Decimal],
    target_product: int,
) -> StoichiometryResult:
    """按 n/系数比确定唯一限量试剂并计算各质量（分数精确值）。

    coefficients 为已校验的配平结果，顺序为「反应物系数 + 产物系数」。
    """
    n_reactants = len(masses_g)

    # 物质的量 n = m / M；M 与 m 均为有限小数，分数形式即精确的 Decimal 商
    mole_amounts = [
        Fraction(masses_g[i]) / Fraction(reactant_molar_masses[i])
        for i in range(n_reactants)
    ]
    # 单位反应进度份额 = n_i / 系数_i，最小者为限量试剂
    extent_ratios = [
        mole_amounts[i] / coefficients[i] for i in range(n_reactants)
    ]

    smallest = min(extent_ratios)
    tied = [i for i, ratio in enumerate(extent_ratios) if ratio == smallest]
    if len(tied) > 1:
        raise StoichiometryError(
            "non_unique_limiting_reagent",
            "two or more reactants run out simultaneously; "
            "the limiting reagent is not unique",
            {
                "field": "reactant_masses_g",
                "tied_indices": tied,
                "extent_ratios": [_format_fraction(r, _NINE_PLACES) for r in extent_ratios],
            },
        )

    limiting = tied[0]
    extent = smallest  # 反应进度（按一单位配平反应计）

    target_mass = (
        extent
        * coefficients[n_reactants + target_product]
        * Fraction(product_molar_masses[target_product])
    )

    remaining: list[Fraction] = []
    for i in range(n_reactants):
        if i == limiting:
            remaining.append(Fraction(0))
            continue
        reacted_mass = (
            extent * coefficients[i] * Fraction(reactant_molar_masses[i])
        )
        remaining.append(Fraction(masses_g[i]) - reacted_mass)

    return StoichiometryResult(
        limiting_reactant_index=limiting,
        target_product_mass_g=target_mass,
        remaining_reactant_masses_g=remaining,
    )


def format_mass(value: Fraction) -> str:
    """把精确质量分数格式化为固定六位小数字符串（ROUND_HALF_UP）。"""
    return _format_fraction(value, _SIX_PLACES)


def _format_fraction(value: Fraction, quantum: Decimal) -> str:
    with localcontext() as ctx:
        ctx.prec = _CALC_PRECISION
        decimal_value = Decimal(value.numerator) / Decimal(value.denominator)
        rounded = decimal_value.quantize(quantum, rounding=ROUND_HALF_UP)
    if rounded == 0:
        # 避免出现 "-0.000000"
        rounded = Decimal(0).quantize(quantum)
    return str(rounded)
