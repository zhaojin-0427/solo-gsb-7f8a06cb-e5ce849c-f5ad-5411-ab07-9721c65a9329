"""化学方程式配平与化学计量 API。"""

from __future__ import annotations

from fractions import Fraction

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field, StrictInt

from .balancer import BalanceError, balance_equation
from .chemistry import FormulaError, parse_formula
from .reaction_heat import (
    MAX_SPECIES,
    HeatError,
    ReactionHeatError,
    format_rational,
    parse_heat,
    solve_reaction_heat,
    validate_net,
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
    description="按元素守恒配平化学方程式，并支持限量试剂与质量关系的化学计量计算，"
    "以及由已知反应热线性推导目标反应热（Hess 定律）。",
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


class KnownReaction(BaseModel):
    net: dict[str, StrictInt] = Field(
        min_length=1,
        description="物种化学式原文 → 净系数（负=消耗，正=生成，绝对值 1～99）",
    )
    heat_kj: str = Field(
        description="反应热（kJ），可带正负号、不含指数的有限十进制字符串"
    )


class TargetReaction(BaseModel):
    net: dict[str, StrictInt] = Field(
        min_length=1,
        description="目标反应的 net 映射，格式同 known 项",
    )


class ReactionHeatRequest(BaseModel):
    known: list[KnownReaction] = Field(min_length=1, max_length=8)
    target: TargetReaction


class ReactionHeatResponse(BaseModel):
    weights: list[str] = Field(
        description="组合权重（与 known 同序），分母为正的最简 p/q 字符串"
    )
    target_heat_kj: str = Field(
        description="目标反应热（kJ），分母为正的最简 p/q 字符串"
    )


@app.get("/")
def root() -> dict[str, str]:
    return {
        "service": "chemical-equation-balancer",
        "usage": (
            "POST /balance to balance an equation; "
            "POST /stoichiometry for limiting reagent and masses; "
            "POST /reaction-heat to derive a target reaction heat "
            "from known reactions"
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
    known_nets: list[dict[str, int]] = []
    heats_kj: list[Fraction] = []
    species_seen: set[str] = set()

    def check_species_budget(net: dict[str, int], field: str) -> None:
        species_seen.update(net)
        if len(species_seen) > MAX_SPECIES:
            raise HTTPException(
                status_code=422,
                detail={
                    "error": "too_many_species",
                    "field": field,
                    "species_total": len(species_seen),
                    "max_species": MAX_SPECIES,
                    "message": (
                        f"the request involves {len(species_seen)} distinct "
                        f"species, exceeding the limit of {MAX_SPECIES}"
                    ),
                },
            )

    # 逐项校验已知反应：net 语义、heat_kj 格式、物种总数预算（422 带字段路径）
    for index, item in enumerate(request.known):
        net_field = f"known[{index}].net"
        try:
            validate_net(item.net)
        except ReactionHeatError as exc:
            raise _reaction_heat_http_error(net_field, exc) from exc
        try:
            heats_kj.append(parse_heat(item.heat_kj))
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
        known_nets.append(item.net)
        check_species_budget(item.net, net_field)

    # 目标反应做同样的 net 校验与物种预算检查
    try:
        validate_net(request.target.net)
    except ReactionHeatError as exc:
        raise _reaction_heat_http_error("target.net", exc) from exc
    check_species_budget(request.target.net, "target.net")

    # 精确有理数推导：先查已知数据自相矛盾，再查目标可推导性（均为 422）
    try:
        weights, target_heat = solve_reaction_heat(
            known_nets, heats_kj, request.target.net
        )
    except ReactionHeatError as exc:
        detail = {"error": exc.code, "message": str(exc)}
        detail.update(exc.extra)
        raise HTTPException(status_code=422, detail=detail) from exc

    return ReactionHeatResponse(
        weights=[format_rational(w) for w in weights],
        target_heat_kj=format_rational(target_heat),
    )


def _reaction_heat_http_error(
    field: str, exc: ReactionHeatError
) -> HTTPException:
    """把 net 校验错误包装成带字段路径的 422。"""
    detail = {"error": exc.code, "field": field, "message": str(exc)}
    detail.update(exc.extra)
    return HTTPException(status_code=422, detail=detail)


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
