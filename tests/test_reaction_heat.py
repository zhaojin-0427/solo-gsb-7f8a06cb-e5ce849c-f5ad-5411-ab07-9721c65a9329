"""POST /reaction-heat 的单元与 HTTP 测试。"""

import json
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.reaction_heat import (
    HeatError,
    ReactionHeatError,
    format_rational,
    parse_heat,
    solve_reaction_heat,
)

client = TestClient(app)


# ---------- 热量字符串解析 ----------

@pytest.mark.parametrize(
    "text, expected",
    [
        ("0", Fraction(0)),
        ("-0", Fraction(0)),
        ("1", Fraction(1)),
        ("+1", Fraction(1)),
        ("-285.83", Fraction(-28583, 100)),
        ("+0.5", Fraction(1, 2)),
        ("0.001", Fraction(1, 1000)),
        ("123456789.987654321", Fraction(123456789987654321, 1000000000)),
    ],
)
def test_parse_heat_valid(text, expected):
    assert parse_heat(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "1e3", "1E-2", "NaN", "Infinity", "-Infinity", "abc", "",
        ".5", "1.", "01", "-01", "+", "-", "1.2.3", " 1", "1 ", "1_0",
    ],
)
def test_parse_heat_invalid(text):
    with pytest.raises(HeatError):
        parse_heat(text)


def test_parse_heat_rejects_non_string():
    for bad in (1.5, 2, None, True):
        with pytest.raises(HeatError):
            parse_heat(bad)


# ---------- p/q 格式化 ----------

@pytest.mark.parametrize(
    "value, expected",
    [
        (Fraction(0), "0/1"),
        (Fraction(2), "2/1"),
        (Fraction(-2), "-2/1"),
        (Fraction(3, 2), "3/2"),
        (Fraction(-4, 2), "-2/1"),   # 自动约分
        (Fraction(1, -3), "-1/3"),   # 分母恒为正
    ],
)
def test_format_rational(value, expected):
    assert format_rational(value) == expected


# ---------- 精确求解 ----------

def test_solve_identity():
    result = solve_reaction_heat(
        known_nets=[{"C": -1, "O2": -1, "CO2": 1}],
        heats_kj=[Fraction(-787, 2)],
        target_net={"C": -1, "O2": -1, "CO2": 1},
    )
    assert result.weights == [Fraction(1)]
    assert result.heat_kj == Fraction(-787, 2)


def test_solve_hess_combination():
    # 2*(C+O2->CO2) - (2CO+O2->2CO2) = 2C+O2->2CO
    result = solve_reaction_heat(
        known_nets=[
            {"C": -1, "O2": -1, "CO2": 1},
            {"CO": -2, "O2": -1, "CO2": 2},
        ],
        heats_kj=[Fraction(-787, 2), Fraction(-566)],
        target_net={"C": -2, "O2": -1, "CO": 2},
    )
    assert result.weights == [Fraction(2), Fraction(-1)]
    assert result.heat_kj == Fraction(-221)


def test_solve_non_unique_weights_still_succeeds():
    # 两个已知反应成比例，权重不唯一，但热量一致，必须成功
    result = solve_reaction_heat(
        known_nets=[
            {"H2": -2, "O2": -1, "H2O": 2},
            {"H2": -4, "O2": -2, "H2O": 4},
        ],
        heats_kj=[Fraction(-28583, 50), Fraction(-28583, 25)],
        target_net={"H2": -2, "O2": -1, "H2O": 2},
    )
    # 任一组合都必须满足 w0 + 2*w1 = 1 且热量唯一
    assert result.weights[0] + 2 * result.weights[1] == 1
    assert result.heat_kj == Fraction(-28583, 50)


def test_solve_inconsistent_cycle():
    # known[1] 是 known[0] 的逆反应，但热量不互为相反数
    with pytest.raises(ReactionHeatError) as exc:
        solve_reaction_heat(
            known_nets=[
                {"H2": -2, "O2": -1, "H2O": 2},
                {"H2": 2, "O2": 1, "H2O": -2},
            ],
            heats_kj=[Fraction(-28583, 50), Fraction(-500)],
            target_net={"H2": -2, "O2": -1, "H2O": 2},
        )
    assert exc.value.code == "inconsistent_known_reactions"
    weights = exc.value.extra["weights"]
    heat = exc.value.extra["heat_kj"]
    # 见证组合：净物质变化为零但热量非零
    assert weights in (["1/1", "1/1"], ["-1/1", "-1/1"])
    assert heat in ("-53583/50", "53583/50")


def test_solve_consistent_cycle_is_allowed():
    # 逆反应热量互为相反数：零净变化组合的热量也为零，数据一致
    result = solve_reaction_heat(
        known_nets=[
            {"H2": -2, "O2": -1, "H2O": 2},
            {"H2": 2, "O2": 1, "H2O": -2},
        ],
        heats_kj=[Fraction(-28583, 50), Fraction(28583, 50)],
        target_net={"H2": -2, "O2": -1, "H2O": 2},
    )
    assert result.heat_kj == Fraction(-28583, 50)


def test_solve_target_not_derivable():
    with pytest.raises(ReactionHeatError) as exc:
        solve_reaction_heat(
            known_nets=[{"H2": -2, "O2": -1, "H2O": 2}],
            heats_kj=[Fraction(-28583, 50)],
            target_net={"C": -1, "O2": -1, "CO2": 1},
        )
    assert exc.value.code == "target_not_derivable"


# ---------- HTTP API：成功 ----------

def post(payload, raw=False):
    if raw:
        return client.post(
            "/reaction-heat",
            content=json.dumps(payload),
            headers={"content-type": "application/json"},
        )
    return client.post("/reaction-heat", json=payload)


def hess_payload():
    return {
        "known": [
            {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "-393.5"},
            {"net": {"CO": -2, "O2": -1, "CO2": 2}, "heat_kj": "-566"},
        ],
        "target": {"net": {"C": -2, "O2": -1, "CO": 2}},
    }


def test_api_hess_combination():
    resp = post(hess_payload())
    assert resp.status_code == 200
    assert resp.json() == {"weights": ["2/1", "-1/1"], "heat_kj": "-221/1"}


def test_api_weights_follow_known_order_and_allow_zero():
    # 第三个已知反应与目标无关，权重为 0/1；权重顺序与 known 一致
    payload = hess_payload()
    payload["known"].append(
        {"net": {"Na": -1, "Cl": -1, "NaCl": 1}, "heat_kj": "30"}
    )
    resp = post(payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "weights": ["2/1", "-1/1", "0/1"],
        "heat_kj": "-221/1",
    }


def test_api_non_unique_weights_succeeds():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.66"},
                {"net": {"H2": -4, "O2": -2, "H2O": 4}, "heat_kj": "-1143.32"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["heat_kj"] == "-28583/50"
    # 返回的权重必须是合法组合：w0 + 2*w1 = 1
    w0, w1 = (Fraction(w) for w in body["weights"])
    assert w0 + 2 * w1 == 1


def test_api_fractional_weights():
    # known: 2A+2B->2AB（即 2 倍），target: A+B->AB，权重 1/2
    resp = post(
        {
            "known": [
                {"net": {"A": -2, "B": -2, "AB": 2}, "heat_kj": "-80"},
            ],
            "target": {"net": {"A": -1, "B": -1, "AB": 1}},
        }
    )
    assert resp.status_code == 200
    assert resp.json() == {"weights": ["1/2"], "heat_kj": "-40/1"}


def test_api_twelve_species_boundary_ok():
    known = [
        {"net": {"A": -1, "B": -1, "AB": 1}, "heat_kj": "1"},
        {"net": {"C": -1, "D": -1, "CD": 1}, "heat_kj": "2"},
        {"net": {"E": -1, "F": -1, "EF": 1}, "heat_kj": "3"},
        {"net": {"G": -1, "H": -1, "GH": 1}, "heat_kj": "4"},
    ]
    target = {
        "net": {
            "A": -1, "B": -1, "AB": 1,
            "C": -1, "D": -1, "CD": 1,
            "E": -1, "F": -1, "EF": 1,
            "G": -1, "H": -1, "GH": 1,
        }
    }
    resp = post({"known": known, "target": target})
    assert resp.status_code == 200
    assert resp.json() == {
        "weights": ["1/1", "1/1", "1/1", "1/1"],
        "heat_kj": "10/1",
    }


def test_api_coefficient_boundaries_99_ok():
    resp = post(
        {
            "known": [
                {"net": {"A": -99, "B": -99, "AB": 99}, "heat_kj": "+99"},
            ],
            "target": {"net": {"A": -99, "B": -99, "AB": 99}},
        }
    )
    assert resp.status_code == 200
    assert resp.json() == {"weights": ["1/1"], "heat_kj": "99/1"}


# ---------- HTTP API：数据矛盾 / 不可推导 ----------

def test_api_inconsistent_known_reactions():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.66"},
                {"net": {"H2": 2, "O2": 1, "H2O": -2}, "heat_kj": "-500"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "inconsistent_known_reactions"
    assert detail["weights"] == ["1/1", "1/1"]
    assert detail["heat_kj"] == "-53583/50"


def test_api_target_not_derivable():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.66"},
            ],
            "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
        }
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "target_not_derivable"


# ---------- HTTP API：net 校验（422，带字段路径） ----------

BASE_KNOWN = [{"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.66"}]
BASE_TARGET = {"net": {"H2": -2, "O2": -1, "H2O": 2}}


def test_api_invalid_formula_located_in_known():
    resp = post(
        {
            "known": [
                {"net": {"H2": -1, "O2": -1, "H2O": 1}, "heat_kj": "1"},
                {"net": {"H2": -1, "H0": 1}, "heat_kj": "1"},
            ],
            "target": BASE_TARGET,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["field"] == "known[1].net"
    assert detail["formula"] == "H0"


def test_api_invalid_formula_located_in_target():
    resp = post({"known": BASE_KNOWN, "target": {"net": {"H2": -1, "(OH": 1}}})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["field"] == "target.net"
    assert detail["formula"] == "(OH"


@pytest.mark.parametrize("bad", [0, 100, -100])
def test_api_coefficient_out_of_range_located(bad):
    resp = post(
        {
            "known": [{"net": {"H2": bad, "H2O": 1}, "heat_kj": "1"}],
            "target": BASE_TARGET,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_coefficient"
    assert detail["field"] == "known[0].net"
    assert detail["formula"] == "H2"
    assert detail["value"] == bad


@pytest.mark.parametrize(
    "bad", ["1e3", "abc", "", ".5", "1.", "01", "NaN", "Infinity", "+-1"]
)
def test_api_invalid_heat_located(bad):
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": bad},
            ],
            "target": BASE_TARGET,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_heat"
    assert detail["field"] == "known[0].heat_kj"
    assert detail["value"] == bad


def test_api_unbalanced_known_located():
    resp = post(
        {
            "known": [{"net": {"H2": -1, "O2": -1, "H2O": 1}, "heat_kj": "1"}],
            "target": BASE_TARGET,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "unbalanced_reaction"
    assert detail["field"] == "known[0].net"
    assert detail["element_changes"] == {"O": -1}


def test_api_unbalanced_target_located():
    resp = post(
        {"known": BASE_KNOWN, "target": {"net": {"H2": -1, "O2": -1, "H2O": 1}}}
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "unbalanced_reaction"
    assert detail["field"] == "target.net"


def test_api_empty_reactant_side():
    resp = post(
        {
            "known": [{"net": {"H2": 2, "H2O": 1}, "heat_kj": "1"}],
            "target": BASE_TARGET,
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "empty_side"
    assert detail["field"] == "known[0].net"
    assert detail["side"] == "reactants"


def test_api_empty_product_side():
    resp = post({"known": BASE_KNOWN, "target": {"net": {"H2": -2, "H2O": -1}}})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "empty_side"
    assert detail["field"] == "target.net"
    assert detail["side"] == "products"


def test_api_too_many_species():
    # 每个 net 都不超过 12 项，但不同物质总计 3 + 12 = 15 种
    known = [{"net": {"A": -1, "B": -1, "AB": 1}, "heat_kj": "10"}]
    target = {
        "net": {
            "C": -1, "D": -1, "CD": 1,
            "E": -1, "F": -1, "EF": 1,
            "G": -1, "H": -1, "GH": 1,
            "I": -1, "J": -1, "IJ": 1,
        }
    }
    resp = post({"known": known, "target": target})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "too_many_species"
    assert detail["species_count"] == 15
    assert detail["max_species"] == 12


# ---------- HTTP API：整体请求结构（Pydantic 422） ----------

@pytest.mark.parametrize(
    "payload, loc_prefix",
    [
        ({"known": [], "target": BASE_TARGET}, ["body", "known"]),
        ({"known": BASE_KNOWN * 9, "target": BASE_TARGET}, ["body", "known"]),
        ({"known": [{"net": {}, "heat_kj": "1"}], "target": BASE_TARGET},
         ["body", "known", 0, "net"]),
        ({"known": [{"net": {"H2": -1, "O2": -1}}], "target": BASE_TARGET},
         ["body", "known", 0, "heat_kj"]),
        ({"known": BASE_KNOWN}, ["body", "target"]),
    ],
)
def test_api_schema_limits(payload, loc_prefix):
    resp = post(payload)
    assert resp.status_code == 422
    assert resp.json()["detail"][0]["loc"] == loc_prefix


@pytest.mark.parametrize("bad", [True, 1.5, "2"])
def test_api_non_integer_coefficient_rejected_by_schema(bad):
    resp = post(
        {
            "known": [{"net": {"H2": bad, "O2": -1}, "heat_kj": "1"}],
            "target": BASE_TARGET,
        },
        raw=True,
    )
    assert resp.status_code == 422
    loc = resp.json()["detail"][0]["loc"]
    assert loc == ["body", "known", 0, "net", "H2"]


def test_api_heat_must_be_string():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": 1.5},
            ],
            "target": BASE_TARGET,
        },
        raw=True,
    )
    assert resp.status_code == 422
    loc = resp.json()["detail"][0]["loc"]
    assert loc == ["body", "known", 0, "heat_kj"]


def test_api_existing_endpoints_untouched():
    resp = client.post(
        "/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    assert resp.status_code == 200
    assert resp.json()["reactants"] == [2, 1]
    resp = client.get("/")
    assert resp.status_code == 200
