"""POST /stoichiometry 的单元与 HTTP 测试。"""

import json
import warnings
from decimal import Decimal
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.chemistry import parse_formula
from app.main import app
from app.stoichiometry import (
    ATOMIC_MASSES,
    MassError,
    StoichiometryError,
    UnsupportedElementError,
    format_mass,
    molar_mass,
    parse_mass,
    solve_stoichiometry,
)

warnings.filterwarnings("ignore")
client = TestClient(app)


# ---------- 摩尔质量（Decimal） ----------

@pytest.mark.parametrize(
    "formula, expected",
    [
        ("H2O", "18.015"),
        ("O2", "31.998"),
        ("H2", "2.016"),
        ("NaCl", "58.44"),
        ("Cl2", "70.90"),
        ("Ca(OH)2", "74.092"),
        ("H3PO4", "97.994"),
        ("Ca3(PO4)2", "310.174"),
        ("(NH4)2SO4", "132.134"),
    ],
)
def test_molar_mass(formula, expected):
    assert molar_mass(parse_formula(formula)) == Decimal(expected)


def test_molar_mass_unsupported_element():
    with pytest.raises(UnsupportedElementError) as exc:
        molar_mass(parse_formula("Fe2(SO4)3"))
    assert exc.value.element == "Fe"


def test_atomic_weight_table():
    assert ATOMIC_MASSES == {
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


def test_nested_parentheses_molar_mass():
    # K4(ON(SO3)2)2 的括号嵌套：4K + 2N + 14O + 4S
    # 4*39.098 + 2*14.007 + 14*15.999 + 4*32.06 = 536.632
    assert molar_mass(parse_formula("K4(ON(SO3)2)2")) == Decimal("536.632")


# ---------- 质量字符串解析 ----------

@pytest.mark.parametrize(
    "text",
    ["0.001", "1", "1.0", "123.456789", "1000000000", "1000000000.000", "0.0010"],
)
def test_parse_mass_valid(text):
    value = parse_mass(text)
    assert value == Decimal(text)
    assert Decimal("0.001") <= value <= Decimal("1000000000")


@pytest.mark.parametrize(
    "text",
    [
        "1e-3", "1E3", "NaN", "Infinity", "-Infinity",
        "-1", "+2.0", "01", ".5", "1.", " 2", "2 ", "",
        "0", "0.0009", "1000000001", "1_000",
    ],
)
def test_parse_mass_invalid(text):
    with pytest.raises(MassError):
        parse_mass(text)


# ---------- 精确求解 ----------

def test_solve_simple_limiting():
    # 2H2 + O2 -> 2H2O；2g H2 + 20g O2，H2 限量
    result = solve_stoichiometry(
        coefficients=[2, 1, 2],
        masses_g=[Decimal("2"), Decimal("20")],
        reactant_molar_masses=[Decimal("2.016"), Decimal("31.998")],
        product_molar_masses=[Decimal("18.015")],
        target_product=0,
    )
    assert result.limiting_reactant_index == 0
    assert format_mass(result.target_product_mass_g) == "17.872024"
    assert [format_mass(m) for m in result.remaining_reactant_masses_g] == [
        "0.000000",
        "4.127976",
    ]


def test_solve_other_reagent_limiting():
    # 1.2g H2 + 20g O2，仍是 H2 限量；换成 0.1g H2 + 20g O2
    result = solve_stoichiometry(
        coefficients=[2, 1, 2],
        masses_g=[Decimal("0.1"), Decimal("20")],
        reactant_molar_masses=[Decimal("2.016"), Decimal("31.998")],
        product_molar_masses=[Decimal("18.015")],
        target_product=0,
    )
    assert result.limiting_reactant_index == 0
    # 大过量 O2：反应 0.1g H2 只需约 0.793601g O2
    assert format_mass(result.remaining_reactant_masses_g[1]) == "19.206399"
    assert format_mass(result.target_product_mass_g) == "0.893601"


def test_solve_single_reactant():
    # 2KClO3 -> 2KCl + 3O2，122.545g KClO3（恰为 1mol，半单位反应）
    result = solve_stoichiometry(
        coefficients=[2, 2, 3],
        masses_g=[Decimal("122.545")],
        reactant_molar_masses=[Decimal("122.545")],
        product_molar_masses=[Decimal("74.548"), Decimal("31.998")],
        target_product=1,
    )
    assert result.limiting_reactant_index == 0
    assert format_mass(result.target_product_mass_g) == "47.997000"
    assert format_mass(result.remaining_reactant_masses_g[0]) == "0.000000"


def test_solve_tie_is_exact():
    # CH4 + 2O2 -> ...；1mol CH4(16.043g) 与 2mol O2(63.996g) 同时耗尽
    with pytest.raises(StoichiometryError) as exc:
        solve_stoichiometry(
            coefficients=[1, 2, 1, 2],
            masses_g=[Decimal("16.043"), Decimal("63.996")],
            reactant_molar_masses=[Decimal("16.043"), Decimal("31.998")],
            product_molar_masses=[Decimal("44.009"), Decimal("18.015")],
            target_product=0,
        )
    assert exc.value.code == "non_unique_limiting_reagent"
    assert exc.value.extra["tied_indices"] == [0, 1]


def test_format_mass_six_places_and_half_up():
    assert format_mass(Fraction(0)) == "0.000000"
    assert format_mass(Fraction(1, 400000)) == "0.000003"  # 0.0000025 四舍五入
    assert format_mass(Fraction(17872024, 1000000)) == "17.872024"


# ---------- HTTP API：成功 ----------

BASE = {"reactants": ["H2", "O2"], "products": ["H2O"]}


def post(payload, raw=False):
    if raw:
        return client.post(
            "/stoichiometry",
            content=json.dumps(payload),
            headers={"content-type": "application/json"},
        )
    return client.post("/stoichiometry", json=payload)


def test_api_ok_h2_limiting():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body == {
        "limiting_reactant_index": 0,
        "target_product_mass_g": "17.872024",
        "remaining_reactant_masses_g": ["0.000000", "4.127976"],
    }


def test_api_complex_equation_and_target_switch():
    payload = {
        "reactants": ["Ca(OH)2", "H3PO4"],
        "products": ["Ca3(PO4)2", "H2O"],
        "coefficients": [3, 2, 1, 6],
        "reactant_masses_g": ["74.094", "100"],
    }
    resp = post({**payload, "target_product": 0})
    assert resp.status_code == 200
    body = resp.json()
    assert body["limiting_reactant_index"] == 0
    assert body["target_product_mass_g"] == "103.394124"
    assert body["remaining_reactant_masses_g"] == ["0.000000", "34.668903"]

    resp = post({**payload, "target_product": 1})
    assert resp.status_code == 200
    assert resp.json()["target_product_mass_g"] == "36.030973"


def test_api_boundary_masses():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["0.001", "1000000000"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["limiting_reactant_index"] == 0
    assert body["target_product_mass_g"] == "0.008936"
    assert body["remaining_reactant_masses_g"] == [
        "0.000000",
        "999999999.992064",
    ]


def test_api_second_reactant_limiting():
    # Na2CO3 + 2HCl -> 2NaCl + H2O + CO2，HCl 限量
    resp = post(
        {
            "reactants": ["Na2CO3", "HCl"],
            "products": ["NaCl", "H2O", "CO2"],
            "coefficients": [1, 2, 2, 1, 1],
            "reactant_masses_g": ["105.988", "70"],
            "target_product": 2,
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["limiting_reactant_index"] == 1
    assert body["target_product_mass_g"] == "42.249026"
    assert body["remaining_reactant_masses_g"] == ["4.238590", "0.000000"]


# ---------- HTTP API：并列限量 ----------

def test_api_tie_returns_422_with_indices():
    resp = post(
        {
            "reactants": ["CH4", "O2"],
            "products": ["CO2", "H2O"],
            "coefficients": [1, 2, 1, 2],
            "reactant_masses_g": ["16.043", "63.996"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "non_unique_limiting_reagent"
    assert detail["field"] == "reactant_masses_g"
    assert detail["tied_indices"] == [0, 1]
    assert detail["extent_ratios"] == ["1.000000000", "1.000000000"]


# ---------- HTTP API：系数问题 ----------

def test_api_coefficient_mismatch_locates_entries():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 1],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "coefficient_mismatch"
    assert detail["field"] == "coefficients"
    assert detail["expected"] == [2, 1, 2]
    assert detail["given"] == [2, 1, 1]
    assert detail["mismatches"] == [
        {
            "index": 2,
            "side": "products",
            "side_index": 0,
            "given": 1,
            "expected": 2,
        }
    ]


def test_api_non_coprime_coefficients_rejected():
    # [4,2,4] 与最简结果 [2,1,2] 比例相同但不是 /balance 的返回值
    resp = post(
        {
            **BASE,
            "coefficients": [4, 2, 4],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "coefficient_mismatch"
    assert {m["index"] for m in detail["mismatches"]} == {0, 1, 2}


def test_api_non_positive_coefficient_located():
    for bad in (0, -2):
        resp = post(
            {
                **BASE,
                "coefficients": [bad, 1, 2],
                "reactant_masses_g": ["2", "20"],
                "target_product": 0,
            }
        )
        assert resp.status_code == 400
        detail = resp.json()["detail"]
        assert detail["error"] == "invalid_coefficient"
        assert detail["field"] == "coefficients"
        assert detail["index"] == 0
        assert detail["value"] == bad


def test_api_coefficient_length_mismatch():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_array_length"
    assert detail["field"] == "coefficients"
    assert detail["expected_length"] == 3
    assert detail["actual_length"] == 2


def test_api_boolean_and_float_coefficients_rejected_by_schema():
    for bad in (True, 2.0):
        resp = post(
            {
                **BASE,
                "coefficients": [bad, 1, 2],
                "reactant_masses_g": ["2", "20"],
                "target_product": 0,
            },
            raw=True,
        )
        assert 400 <= resp.status_code < 500
        loc = resp.json()["detail"][0]["loc"]
        assert loc[:2] == ["body", "coefficients"]


# ---------- HTTP API：质量问题 ----------

@pytest.mark.parametrize("bad", ["1e-3", "NaN", "Infinity", "0.0009", "1000000001",
                                 "-1", "01", ".5", "1.", "", "+2"])
def test_api_invalid_mass_located(bad):
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": [bad, "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_mass"
    assert detail["field"] == "reactant_masses_g"
    assert detail["index"] == 0
    assert detail["value"] == bad


def test_api_invalid_mass_second_index():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "0"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["field"] == "reactant_masses_g"
    assert detail["index"] == 1


def test_api_mass_length_must_equal_reactants():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_array_length"
    assert detail["field"] == "reactant_masses_g"
    assert detail["expected_length"] == 2
    assert detail["actual_length"] == 1


def test_api_mass_wrong_type_has_body_loc():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": [2, "20"],
            "target_product": 0,
        },
        raw=True,
    )
    assert resp.status_code == 422
    loc = resp.json()["detail"][0]["loc"]
    assert loc == ["body", "reactant_masses_g", 0]


# ---------- HTTP API：target_product ----------

def test_api_target_product_out_of_range_located():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": 1,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "index_out_of_range"
    assert detail["field"] == "target_product"
    assert detail["index"] == 1
    assert detail["valid_range"] == [0, 0]


def test_api_target_product_negative_is_4xx():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": -1,
        },
        raw=True,
    )
    assert 400 <= resp.status_code < 500
    assert resp.json()["detail"][0]["loc"] == ["body", "target_product"]


# ---------- HTTP API：化学式 ----------

def test_api_invalid_formula_located():
    resp = post(
        {
            "reactants": ["H2", "O2"],
            "products": ["H100"],
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["array"] == "products"
    assert detail["index"] == 0


def test_api_unsupported_element_located():
    resp = post(
        {
            "reactants": ["Fe", "O2"],
            "products": ["Fe2O3"],
            "coefficients": [4, 3, 2],
            "reactant_masses_g": ["10", "10"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "unsupported_element"
    assert detail["array"] == "reactants"
    assert detail["index"] == 0
    assert detail["element"] == "Fe"
    assert set(detail["supported_elements"]) == set(ATOMIC_MASSES)


def test_api_unsupported_element_in_second_reactant():
    resp = post(
        {
            "reactants": ["H2", "Mg"],
            "products": ["H2O"],
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["array"] == "reactants"
    assert detail["index"] == 1
    assert detail["element"] == "Mg"


# ---------- HTTP API：配平失败 ----------

def test_api_unbalanceable_equation():
    resp = post(
        {
            "reactants": ["H2", "O2"],
            "products": ["O2"],
            "coefficients": [2, 1, 2],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "no_positive_solution"


def test_api_non_unique_equation():
    resp = post(
        {
            "reactants": ["C", "O2"],
            "products": ["CO", "CO2"],
            "coefficients": [1, 1, 1, 1],
            "reactant_masses_g": ["2", "20"],
            "target_product": 0,
        }
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "non_unique_solution"


# ---------- HTTP API：整体请求结构 ----------

def test_api_missing_field_is_4xx():
    resp = post(
        {
            **BASE,
            "coefficients": [2, 1, 2],
            "target_product": 0,
            # 缺少 reactant_masses_g
        }
    )
    assert 400 <= resp.status_code < 500
    assert resp.json()["detail"][0]["loc"][:2] == ["body", "reactant_masses_g"]


@pytest.mark.parametrize(
    "payload",
    [
        {"reactants": [], "products": ["H2O"], "coefficients": [1, 1],
         "reactant_masses_g": [], "target_product": 0},
        {"reactants": ["H2"], "products": ["H2O"], "coefficients": [],
         "reactant_masses_g": ["1"], "target_product": 0},
    ],
)
def test_api_schema_limits(payload):
    resp = post(payload)
    assert 400 <= resp.status_code < 500
