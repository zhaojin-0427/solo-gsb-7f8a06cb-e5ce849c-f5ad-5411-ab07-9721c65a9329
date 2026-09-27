"""化学方程式配平与化学计量 API。"""

from __future__ import annotations

from fractions import Fraction

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field, StrictInt

from .balancer import BalanceError, balance_equation
from .chemistry import FormulaError, parse_formula
from .reaction_heat import (
    HeatError,
    ReactionHeatError,
    format_rational,
    parse_heat,
    solve_reaction_heat,
)
from .stoichiometry import (
    SUPPORTED_ELEMENTS,
    MassError,
    StoichiometryError,
    UnsupportedElementError,
    format_mass,
    molar_mass,
    parse_mass,
    solve_stoichiometry,
)

app = FastAPI(
    title="Chemical Equation Balancer",
    version="1.2.0",
    description="按元素守恒配平化学方程式，支持限量试剂与质量关系的化学计量计算，"
    "并可由已知反应热按 Hess 定律推导目标反应热。",
)


class BalanceRequest(BaseModel):
    reactants: list[str] = Field(min_length=1, max_length=4)
    products: list[str] = Field(min_length=1, max_length=4)


class BalanceResponse(BaseModel):
    reactants: list[int]
    products: list[int]
    equation: str


class StoichiometryRequest(BaseModel):
    reactants: list[str] = Field(min_length=1, max_length=4)
    products: list[str] = Field(min_length=1, max_length=4)
    coefficients: list[StrictInt] = Field(
        min_length=2,
        max_length=8,
        description="严格正整数系数，顺序为「反应物系数 + 产物系数」，"
        "须与 /balance 返回的唯一最简结果一致。",
    )
    reactant_masses_g: list[str] = Field(
        min_length=1,
        max_length=4,
        description="与 reactants 等长的克质量数组，取值 0.001～1e9 的十进制字符串。",
    )
    target_product: StrictInt = Field(
        ge=0,
        le=3,
        description="目标产物在 products 中的下标（从 0 开始）。",
    )


class StoichiometryResponse(BaseModel):
    limiting_reactant_index: int
    target_product_mass_g: str = Field(description="目标产物质量（克，六位小数）")
    remaining_reactant_masses_g: list[str] = Field(
        description="各反应物剩余质量（克，六位小数），顺序同 reactants"
    )


# /reaction-heat 的约束：已知反应 1～8 个、全部反应涉及的不同物质总计最多 12 种、
# net 系数为绝对值 1～99 的非零整数
MAX_KNOWN_REACTIONS = 8
MAX_NET_SPECIES = 12
MAX_TOTAL_SPECIES = 12
MAX_NET_COEFFICIENT = 99


class KnownReactionHeat(BaseModel):
    net: dict[str, StrictInt] = Field(
        min_length=1,
        max_length=MAX_NET_SPECIES,
        description="物种净系数映射：化学式原文 -> 系数"
        "（负整数 = 消耗，正整数 = 生成，绝对值 1～99）",
    )
    heat_kj: str = Field(
        description="反应热（kJ），可带正负号、不含指数的有限十进制字符串"
    )


class TargetReaction(BaseModel):
    net: dict[str, StrictInt] = Field(
        min_length=1,
        max_length=MAX_NET_SPECIES,
        description="目标反应的物种净系数映射，约定同 known[].net",
    )


class ReactionHeatRequest(BaseModel):
    known: list[KnownReactionHeat] = Field(
        min_length=1, max_length=MAX_KNOWN_REACTIONS
    )
    target: TargetReaction


class ReactionHeatResponse(BaseModel):
    weights: list[str] = Field(
        description="各已知反应的权重（按 known 顺序），分母为正的最简 p/q 字符串"
    )
    heat_kj: str = Field(description="目标反应热（kJ），分母为正的最简 p/q 字符串")


@app.get("/")
def root() -> dict[str, str]:
    return {
        "service": "chemical-equation-balancer",
        "usage": (
            "POST /balance to balance an equation; "
            "POST /stoichiometry for limiting reagent and masses; "
            "POST /reaction-heat to derive a reaction heat from known reactions"
        ),
        "docs": "/docs",
    }


@app.post("/balance", response_model=BalanceResponse)
def balance(request: BalanceRequest) -> BalanceResponse:
    parsed_reactants, parsed_products = _parse_species(
        request.reactants, request.products
    )

    try:
        reactant_coeffs, product_coeffs = balance_equation(
            parsed_reactants, parsed_products
        )
    except BalanceError as exc:
        raise HTTPException(
            status_code=422,
            detail={"error": exc.code, "message": str(exc)},
        ) from exc

    return BalanceResponse(
        reactants=reactant_coeffs,
        products=product_coeffs,
        equation=_format_equation(
            request.reactants, reactant_coeffs, request.products, product_coeffs
        ),
    )


@app.post("/stoichiometry", response_model=StoichiometryResponse)
def stoichiometry(request: StoichiometryRequest) -> StoichiometryResponse:
    parsed_reactants, parsed_products = _parse_species(
        request.reactants,
        request.products,
        require_supported_elements=True,
    )

    n_reactants = len(request.reactants)
    n_products = len(request.products)

    # 系数与反应物质量的数组长度定位错误（400）
    if len(request.coefficients) != n_reactants + n_products:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "invalid_array_length",
                "field": "coefficients",
                "expected_length": n_reactants + n_products,
                "actual_length": len(request.coefficients),
                "message": (
                    f"coefficients must have {n_reactants + n_products} entries "
                    f"({n_reactants} reactants + {n_products} products)"
                ),
            },
        )
    if len(request.reactant_masses_g) != n_reactants:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "invalid_array_length",
                "field": "reactant_masses_g",
                "expected_length": n_reactants,
                "actual_length": len(request.reactant_masses_g),
                "message": (
                    "reactant_masses_g must have the same length as reactants"
                ),
            },
        )

    # target_product 上界随 products 实际长度收紧（400，带字段定位）
    if not 0 <= request.target_product < n_products:
        raise HTTPException(
            status_code=400,
            detail={
                "error": "index_out_of_range",
                "field": "target_product",
                "index": request.target_product,
                "valid_range": [0, n_products - 1],
                "message": (
                    f"target_product index {request.target_product} is out of range "
                    f"for {n_products} products"
                ),
            },
        )

    # 系数必须全部为严格正整数（类型由 schema 保证，此处校验取值，400 定位下标）
    for flat_index, coefficient in enumerate(request.coefficients):
        if coefficient <= 0:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": "invalid_coefficient",
                    "field": "coefficients",
                    "index": flat_index,
                    "value": coefficient,
                    "message": f"coefficient at index {flat_index} must be >= 1",
                },
            )

    # 解析并校验各反应物质量（400，定位到具体下标）
    masses = []
    for index, mass_text in enumerate(request.reactant_masses_g):
        try:
            masses.append(parse_mass(mass_text))
        except MassError as exc:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": "invalid_mass",
                    "field": "reactant_masses_g",
                    "index": index,
                    "value": mass_text,
                    "message": str(exc),
                },
            ) from exc

    # 配平：/stoichiometry 同样要求存在唯一正整数最简解（422）
    try:
        balanced_reactants, balanced_products = balance_equation(
            parsed_reactants, parsed_products
        )
    except BalanceError as exc:
        raise HTTPException(
            status_code=422,
            detail={"error": exc.code, "message": str(exc)},
        ) from exc
    expected = list(balanced_reactants) + list(balanced_products)

    # 提交的系数必须就是 /balance 的唯一结果（422，逐项定位差异）
    if request.coefficients != expected:
        mismatches = [
            {
                "index": i,
                "side": "reactants" if i < n_reactants else "products",
                "side_index": i if i < n_reactants else i - n_reactants,
                "given": given,
                "expected": wanted,
            }
            for i, (given, wanted) in enumerate(
                zip(request.coefficients, expected)
            )
            if given != wanted
        ]
        raise HTTPException(
            status_code=422,
            detail={
                "error": "coefficient_mismatch",
                "field": "coefficients",
                "expected": expected,
                "given": request.coefficients,
                "mismatches": mismatches,
                "message": (
                    "coefficients must equal the unique /balance result "
                    "(strictly positive, coprime integers in input order)"
                ),
            },
        )

    reactant_molar_masses = [molar_mass(c) for c in parsed_reactants]
    product_molar_masses = [molar_mass(c) for c in parsed_products]

    try:
        result = solve_stoichiometry(
            coefficients=request.coefficients,
            masses_g=masses,
            reactant_molar_masses=reactant_molar_masses,
            product_molar_masses=product_molar_masses,
            target_product=request.target_product,
        )
    except UnsupportedElementError as exc:  # 防御性：正常在 _parse_species 已拦截
        raise HTTPException(
            status_code=400,
            detail={
                "error": "unsupported_element",
                "element": exc.element,
                "supported_elements": SUPPORTED_ELEMENTS,
                "message": str(exc),
            },
        ) from exc
    except StoichiometryError as exc:
        detail = {"error": exc.code, "message": str(exc)}
        detail.update(exc.extra)
        raise HTTPException(status_code=422, detail=detail) from exc

    return StoichiometryResponse(
        limiting_reactant_index=result.limiting_reactant_index,
        target_product_mass_g=format_mass(result.target_product_mass_g),
        remaining_reactant_masses_g=[
            format_mass(m) for m in result.remaining_reactant_masses_g
        ],
    )


@app.post("/reaction-heat", response_model=ReactionHeatResponse)
def reaction_heat(request: ReactionHeatRequest) -> ReactionHeatResponse:
    # 逐项解析 known：net 系数取值与化学式语法、heat_kj 格式（422，带字段路径）
    parsed_known: list[dict[str, dict[str, int]]] = []
    heats: list[Fraction] = []
    for index, item in enumerate(request.known):
        parsed_known.append(_parse_net(item.net, f"known[{index}].net"))
        try:
            heats.append(parse_heat(item.heat_kj))
        except HeatError as exc:
            raise HTTPException(
                status_code=422,
                detail={
                    "error": "invalid_heat",
                    "field": f"known[{index}].heat_kj",
                    "value": item.heat_kj,
                    "message": str(exc),
                },
            ) from exc
    parsed_target = _parse_net(request.target.net, "target.net")

    # 全部反应涉及的不同物质（以化学式原文计）总计不得超过 12 种
    all_species: set[str] = set()
    for item in request.known:
        all_species.update(item.net)
    all_species.update(request.target.net)
    if len(all_species) > MAX_TOTAL_SPECIES:
        raise HTTPException(
            status_code=422,
            detail={
                "error": "too_many_species",
                "species_count": len(all_species),
                "max_species": MAX_TOTAL_SPECIES,
                "message": (
                    f"reactions involve {len(all_species)} distinct species, "
                    f"exceeding the limit of {MAX_TOTAL_SPECIES}"
                ),
            },
        )

    # 逐项校验正负两侧非空与元素守恒（422，带字段路径）
    nets = [item.net for item in request.known] + [request.target.net]
    parsed_nets = parsed_known + [parsed_target]
    fields = [f"known[{i}].net" for i in range(len(request.known))]
    fields.append("target.net")
    for net, parsed, field in zip(nets, parsed_nets, fields):
        _check_net_reaction(net, parsed, field)

    try:
        result = solve_reaction_heat(nets[:-1], heats, request.target.net)
    except ReactionHeatError as exc:
        detail = {"error": exc.code, "message": str(exc)}
        detail.update(exc.extra)
        raise HTTPException(status_code=422, detail=detail) from exc

    return ReactionHeatResponse(
        weights=[format_rational(w) for w in result.weights],
        heat_kj=format_rational(result.heat_kj),
    )


def _parse_net(net: dict[str, int], field: str) -> dict[str, dict[str, int]]:
    """校验 net 映射的系数取值与化学式语法，返回 {物种原文: 元素计数}。"""
    parsed: dict[str, dict[str, int]] = {}
    for formula, coefficient in net.items():
        if coefficient == 0 or abs(coefficient) > MAX_NET_COEFFICIENT:
            raise HTTPException(
                status_code=422,
                detail={
                    "error": "invalid_coefficient",
                    "field": field,
                    "formula": formula,
                    "value": coefficient,
                    "message": (
                        f"coefficient for {formula!r} must be a non-zero "
                        "integer with absolute value in 1.."
                        f"{MAX_NET_COEFFICIENT} "
                        "(negative = consumed, positive = produced)"
                    ),
                },
            )
        try:
            parsed[formula] = parse_formula(formula)
        except FormulaError as exc:
            raise HTTPException(
                status_code=422,
                detail={
                    "error": "invalid_formula",
                    "field": field,
                    "formula": formula,
                    "message": str(exc),
                },
            ) from exc
    return parsed


def _check_net_reaction(
    net: dict[str, int], parsed: dict[str, dict[str, int]], field: str
) -> None:
    """校验单个反应的正负两侧非空与元素守恒（422，带字段路径）。

    先查两侧非空：全正或全负的 net 必然不守恒，若先查守恒则该错误不可达。
    """
    if all(coefficient > 0 for coefficient in net.values()):
        raise HTTPException(
            status_code=422,
            detail={
                "error": "empty_side",
                "field": field,
                "side": "reactants",
                "message": (
                    "reaction has no consumed species "
                    "(negative coefficients)"
                ),
            },
        )
    if all(coefficient < 0 for coefficient in net.values()):
        raise HTTPException(
            status_code=422,
            detail={
                "error": "empty_side",
                "field": field,
                "side": "products",
                "message": (
                    "reaction has no produced species "
                    "(positive coefficients)"
                ),
            },
        )

    element_changes: dict[str, int] = {}
    for formula, coefficient in net.items():
        for element, count in parsed[formula].items():
            element_changes[element] = (
                element_changes.get(element, 0) + coefficient * count
            )
    imbalanced = {e: n for e, n in element_changes.items() if n != 0}
    if imbalanced:
        raise HTTPException(
            status_code=422,
            detail={
                "error": "unbalanced_reaction",
                "field": field,
                "element_changes": imbalanced,
                "message": (
                    "net reaction does not conserve elements: "
                    + ", ".join(
                        f"{element} {change:+d}"
                        for element, change in imbalanced.items()
                    )
                ),
            },
        )


def _parse_species(
    reactant_formulas: list[str],
    product_formulas: list[str],
    require_supported_elements: bool = False,
) -> tuple[list[dict[str, int]], list[dict[str, int]]]:
    """解析两侧化学式，返回元素计数；错误一律带数组名与下标（400）。"""
    parsed_sides: list[list[dict[str, int]]] = []
    for array_name, formulas in (
        ("reactants", reactant_formulas),
        ("products", product_formulas),
    ):
        parsed_side: list[dict[str, int]] = []
        for index, formula in enumerate(formulas):
            try:
                counts = parse_formula(formula)
            except FormulaError as exc:
                raise HTTPException(
                    status_code=400,
                    detail={
                        "error": "invalid_formula",
                        "array": array_name,
                        "index": index,
                        "formula": formula,
                        "message": str(exc),
                    },
                ) from exc

            if require_supported_elements:
                for element in counts:
                    if element not in SUPPORTED_ELEMENTS:
                        raise HTTPException(
                            status_code=400,
                            detail={
                                "error": "unsupported_element",
                                "array": array_name,
                                "index": index,
                                "formula": formula,
                                "element": element,
                                "supported_elements": SUPPORTED_ELEMENTS,
                                "message": (
                                    f"formula {formula!r} contains element "
                                    f"{element!r}, which has no configured "
                                    "atomic weight"
                                ),
                            },
                        )
            parsed_side.append(counts)
        parsed_sides.append(parsed_side)
    return parsed_sides[0], parsed_sides[1]


def _format_equation(
    reactants: list[str],
    reactant_coeffs: list[int],
    products: list[str],
    product_coeffs: list[int],
) -> str:
    def side(formulas: list[str], coeffs: list[int]) -> str:
        terms = [
            formula if coeff == 1 else f"{coeff} {formula}"
            for formula, coeff in zip(formulas, coeffs)
        ]
        return " + ".join(terms)

    return f"{side(reactants, reactant_coeffs)} -> {side(products, product_coeffs)}"
