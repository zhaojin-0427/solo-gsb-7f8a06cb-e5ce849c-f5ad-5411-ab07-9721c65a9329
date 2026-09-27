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
    validate_net,
)

client = TestClient(app)


# ---------- heat_kj 字符串解析 ----------

@pytest.mark.parametrize(
    "text, expected",
    [
        ("0", Fraction(0)),
        ("-0", Fraction(0)),
        ("285.75", Fraction(28575, 100)),
        ("-285.75", Fraction(-28575, 100)),
        ("+571.6", Fraction(5716, 10)),
        ("-1143.2", Fraction(-11432, 10)),
        ("1000000000000", Fraction(10**12)),
        ("0.000000001", Fraction(1, 10**9)),
    ],
)
def test_parse_heat_valid(text, expected):
    assert parse_heat(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "", "abc", "NaN", "Infinity", "-Infinity",
        "1e5", "1E-3", "2.5e2",   # 不允许指数
        ".5", "1.",               # 小数点两侧须完整
        "01", "-01",              # 前导零
        "1_000", "1,5", " 2", "2 ", "--1", "++1", "+-1",
    ],
)
def test_parse_heat_invalid(text):
    with pytest.raises(HeatError):
        parse_heat(text)


def test_format_rational_reduced_positive_denominator():
    assert format_rational(Fraction(6, 4)) == "3/2"
    assert format_rational(Fraction(-6, 4)) == "-3/2"
    assert format_rational(Fraction(6, -4)) == "-3/2"
    assert format_rational(Fraction(0)) == "0/1"
    assert format_rational(Fraction(-221)) == "-221/1"


# ---------- net 校验 ----------

def test_validate_net_ok():
    validate_net({"C": -1, "O2": -1, "CO2": 1})
    validate_net({"Ca(OH)2": -3, "H3PO4": -2, "Ca3(PO4)2": 1, "H2O": 6})


@pytest.mark.parametrize("coefficient", [0, 100, -100, 250])
def test_validate_net_coefficient_range(coefficient):
    with pytest.raises(ReactionHeatError) as exc:
        validate_net({"H2": coefficient, "O2": -1, "H2O": 1})
    assert exc.value.code == "invalid_coefficient"
    assert exc.value.extra["species"] == "H2"
    assert exc.value.extra["value"] == coefficient


def test_validate_net_invalid_formula():
    with pytest.raises(ReactionHeatError) as exc:
        validate_net({"H2": -1, "H0": 1})
    assert exc.value.code == "invalid_formula"
    assert exc.value.extra["species"] == "H0"


def test_validate_net_empty_side():
    with pytest.raises(ReactionHeatError) as exc:
        validate_net({"H2": -1, "O2": -1})
    assert exc.value.code == "empty_side"
    assert exc.value.extra["missing_side"] == "produced"

    with pytest.raises(ReactionHeatError) as exc:
        validate_net({"H2O": 1})
    assert exc.value.code == "empty_side"
    assert exc.value.extra["missing_side"] == "consumed"


def test_validate_net_unbalanced_locates_elements():
    with pytest.raises(ReactionHeatError) as exc:
        validate_net({"H2": -1, "H2O": 1})
    assert exc.value.code == "unbalanced_reaction"
    assert exc.value.extra["imbalance"] == {"O": 1}


# ---------- 精确求解核心 ----------

def test_solve_basic_hess_cycle():
    # 2C + O2 -> 2CO 由 2×(C + O2 -> CO2) 减去 (2CO + O2 -> 2CO2) 得到
    weights, heat = solve_reaction_heat(
        known_nets=[
            {"C": -1, "O2": -1, "CO2": 1},
            {"CO": -2, "O2": -1, "CO2": 2},
        ],
        heats_kj=[Fraction("-393.5"), Fraction("-566")],
        target_net={"C": -2, "O2": -1, "CO": 2},
    )
    assert weights == [Fraction(2), Fraction(-1)]
    assert heat == Fraction(-221)


def test_solve_fractional_weight():
    # 目标恰为已知反应的一半
    weights, heat = solve_reaction_heat(
        known_nets=[{"H2": -4, "O2": -2, "H2O": 4}],
        heats_kj=[Fraction("-1143.2")],
        target_net={"H2": -2, "O2": -1, "H2O": 2},
    )
    assert weights == [Fraction(1, 2)]
    assert heat == Fraction("-571.6")


def test_solve_reverse_reaction_negative_weight():
    weights, heat = solve_reaction_heat(
        known_nets=[{"H2": -2, "O2": -1, "H2O": 2}],
        heats_kj=[Fraction("-571.6")],
        target_net={"H2O": -2, "H2": 2, "O2": 1},
    )
    assert weights == [Fraction(-1)]
    assert heat == Fraction("571.6")


def test_solve_non_unique_weights_succeeds():
    # known[1] 是 known[0] 的二倍且热量一致：权重不唯一但必须成功
    weights, heat = solve_reaction_heat(
        known_nets=[
            {"H2": -2, "O2": -1, "H2O": 2},
            {"H2": -4, "O2": -2, "H2O": 4},
        ],
        heats_kj=[Fraction("-571.6"), Fraction("-1143.2")],
        target_net={"H2": -2, "O2": -1, "H2O": 2},
    )
    assert len(weights) == 2
    # 任一组合法解：净变化与热量都必须吻合
    assert weights[0] + 2 * weights[1] == 1
    assert weights[0] * Fraction("-571.6") + weights[1] * Fraction("-1143.2") == heat
    assert heat == Fraction("-571.6")


def test_solve_inconsistent_known_reports_combination():
    # known[1] 的净变化是 known[0] 的二倍，热量却不是：数据自相矛盾
    with pytest.raises(ReactionHeatError) as exc:
        solve_reaction_heat(
            known_nets=[
                {"H2": -2, "O2": -1, "H2O": 2},
                {"H2": -4, "O2": -2, "H2O": 4},
            ],
            heats_kj=[Fraction("-571.6"), Fraction("-1000")],
            target_net={"H2": -2, "O2": -1, "H2O": 2},
        )
    assert exc.value.code == "inconsistent_known"
    weights = exc.value.extra["weights"]
    assert len(weights) == 2
    # 组合的净物质变化为零、热量非零
    assert exc.value.extra["heat_kj"] != "0/1"


def test_solve_not_derivable():
    with pytest.raises(ReactionHeatError) as exc:
        solve_reaction_heat(
            known_nets=[{"H2": -2, "O2": -1, "H2O": 2}],
            heats_kj=[Fraction("-571.6")],
            target_net={"C": -1, "O2": -1, "CO2": 1},
        )
    assert exc.value.code == "not_derivable"


# ---------- HTTP API ----------

HESS_PAYLOAD = {
    "known": [
        {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "-393.5"},
        {"net": {"CO": -2, "O2": -1, "CO2": 2}, "heat_kj": "-566"},
    ],
    "target": {"net": {"C": -2, "O2": -1, "CO": 2}},
}


def post(payload, raw=False):
    if raw:
        return client.post(
            "/reaction-heat",
            content=json.dumps(payload),
            headers={"content-type": "application/json"},
        )
    return client.post("/reaction-heat", json=payload)


def test_api_ok_basic():
    resp = post(HESS_PAYLOAD)
    assert resp.status_code == 200
    assert resp.json() == {
        "weights": ["2/1", "-1/1"],
        "target_heat_kj": "-221/1",
    }


def test_api_ok_fractional_result():
    resp = post(
        {
            "known": [
                {"net": {"H2": -4, "O2": -2, "H2O": 4}, "heat_kj": "-1143.2"}
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "weights": ["1/2"],
        "target_heat_kj": "-2858/5",
    }


def test_api_ok_zero_weight_and_non_unique():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.6"},
                {"net": {"H2": -4, "O2": -2, "H2O": 4}, "heat_kj": "-1143.2"},
                {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "-393.5"},
            ],
            "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["weights"] == ["0/1", "0/1", "1/1"]
    assert body["target_heat_kj"] == "-787/2"


def test_api_ok_signed_heat_strings():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "+571.6"},
            ],
            "target": {"net": {"H2O": -2, "H2": 2, "O2": 1}},
        }
    )
    assert resp.status_code == 200
    assert resp.json() == {"weights": ["-1/1"], "target_heat_kj": "-2858/5"}


def test_api_inconsistent_known_422_with_weights_and_heat():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.6"},
                {"net": {"H2": -4, "O2": -2, "H2O": 4}, "heat_kj": "-1000"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "inconsistent_known"
    assert detail["weights"] == ["-2/1", "1/1"]
    assert detail["heat_kj"] == "716/5"


def test_api_not_derivable_422():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "-571.6"}
            ],
            "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
        }
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "not_derivable"


# ---------- HTTP API：net 逐项校验（422 带字段路径） ----------

def test_api_invalid_formula_located_in_known():
    payload = {
        "known": [
            {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "-393.5"},
            {"net": {"H2": -1, "H0": 1}, "heat_kj": "0"},
        ],
        "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
    }
    resp = post(payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["field"] == "known[1].net"
    assert detail["species"] == "H0"


def test_api_invalid_formula_located_in_target():
    resp = post({**HESS_PAYLOAD, "target": {"net": {"C": -1, "CO2)": 1}}})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["field"] == "target.net"


@pytest.mark.parametrize("coefficient", [0, 100, -100])
def test_api_coefficient_out_of_range_located(coefficient):
    payload = {
        "known": [
            {
                "net": {"H2": coefficient, "O2": -1, "H2O": 2},
                "heat_kj": "-571.6",
            }
        ],
        "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
    }
    resp = post(payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_coefficient"
    assert detail["field"] == "known[0].net"
    assert detail["species"] == "H2"
    assert detail["value"] == coefficient


def test_api_coefficient_boundary_99_ok():
    resp = post(
        {
            "known": [
                {"net": {"H2": -98, "O2": -49, "H2O": 98}, "heat_kj": "1"},
                {"net": {"H2O": -98, "H2": 98, "O2": 49}, "heat_kj": "-1"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "weights": ["1/49", "0/1"],
        "target_heat_kj": "1/49",
    }


def test_api_empty_side_located():
    resp = post(
        {
            "known": [
                {"net": {"H2": -2, "O2": -1}, "heat_kj": "-571.6"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "empty_side"
    assert detail["field"] == "known[0].net"
    assert detail["missing_side"] == "produced"


def test_api_unbalanced_reaction_located_with_imbalance():
    resp = post(
        {
            "known": [
                {"net": {"H2": -1, "H2O": 1}, "heat_kj": "0"},
            ],
            "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
        }
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "unbalanced_reaction"
    assert detail["field"] == "known[0].net"
    assert detail["imbalance"] == {"O": 1}


def test_api_unbalanced_target_located():
    resp = post({**HESS_PAYLOAD, "target": {"net": {"C": -1, "CO": 1}}})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "unbalanced_reaction"
    assert detail["field"] == "target.net"


# ---------- HTTP API：heat_kj ----------

@pytest.mark.parametrize(
    "bad", ["", "abc", "1e5", "NaN", "Infinity", ".5", "1.", "01", "1_000"]
)
def test_api_invalid_heat_located(bad):
    payload = {
        "known": [
            {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": bad},
        ],
        "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
    }
    resp = post(payload)
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_heat"
    assert detail["field"] == "known[0].heat_kj"
    assert detail["value"] == bad


def test_api_heat_wrong_type_is_4xx():
    payload = {
        "known": [{"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": -393.5}],
        "target": {"net": {"C": -1, "O2": -1, "CO2": 1}},
    }
    resp = post(payload, raw=True)
    assert 400 <= resp.status_code < 500
    loc = resp.json()["detail"][0]["loc"]
    assert loc[:3] == ["body", "known", 0]


# ---------- HTTP API：物种总数 ----------

def test_api_too_many_species_located():
    # 12 种以内合法；target 引入第 13 种物种时定位到 target.net
    known = [
        {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "1"},
        {"net": {"N2": -1, "H2": -3, "NH3": 2}, "heat_kj": "2"},
        {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "3"},
        {"net": {"S": -1, "O2": -1, "SO2": 1}, "heat_kj": "4"},
        {"net": {"Na": -2, "Cl2": -1, "NaCl": 2}, "heat_kj": "5"},
    ]  # 已有 H2 O2 H2O N2 NH3 C CO2 S SO2 Na Cl2 NaCl 共 12 种
    resp = post(
        {**dict(known=known), "target": {"net": {"K": -1, "NaCl": -1, "Na": 1, "KCl": 1}}}
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "too_many_species"
    assert detail["field"] == "target.net"
    assert detail["max_species"] == 12
    assert detail["species_total"] > 12


def test_api_exactly_12_species_ok():
    known = [
        {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "1"},
        {"net": {"N2": -1, "H2": -3, "NH3": 2}, "heat_kj": "2"},
        {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "3"},
        {"net": {"S": -1, "O2": -1, "SO2": 1}, "heat_kj": "4"},
        {"net": {"Na": -2, "Cl2": -1, "NaCl": 2}, "heat_kj": "5"},
    ]  # 恰好 12 种
    resp = post(
        {**dict(known=known), "target": {"net": {"H2O": -2, "H2": 2, "O2": 1}}}
    )
    assert resp.status_code == 200


# ---------- HTTP API：请求结构（schema 422） ----------

@pytest.mark.parametrize(
    "payload",
    [
        {},                                            # 缺 known 与 target
        {"known": [], "target": {"net": {"H2": -1, "H2O": 1}}},   # known 为空
        {"target": {"net": {"H2": -1, "H2O": 1}}},     # 缺 known
        {"known": [{"net": {"H2": -1, "H2O": 1}, "heat_kj": "1"}]},  # 缺 target
        {
            "known": [{"net": {}, "heat_kj": "1"}],
            "target": {"net": {"H2": -1, "H2O": 1}},
        },                                             # net 为空映射
        {
            "known": [{"net": {"H2": -1, "H2O": 1}}],  # 缺 heat_kj
            "target": {"net": {"H2": -1, "H2O": 1}},
        },
        {
            "known": [{"net": {"H2": -1.5, "H2O": 1}, "heat_kj": "1"}],
            "target": {"net": {"H2": -1, "H2O": 1}},
        },                                             # 非整数系数
        {
            "known": [{"net": {"H2": True, "H2O": 1}, "heat_kj": "1"}],
            "target": {"net": {"H2": -1, "H2O": 1}},
        },                                             # 布尔系数
    ],
)
def test_api_schema_validation(payload):
    resp = post(payload, raw=True)
    assert resp.status_code == 422


def test_api_known_max_8_items():
    item = {"net": {"H2": -2, "O2": -1, "H2O": 2}, "heat_kj": "1"}
    payload = {
        "known": [dict(item) for _ in range(9)],
        "target": {"net": {"H2": -2, "O2": -1, "H2O": 2}},
    }
    resp = post(payload)
    assert resp.status_code == 422

    payload["known"] = [dict(item) for _ in range(8)]
    resp = post(payload)
    assert resp.status_code == 200


def test_api_existing_endpoints_untouched():
    resp = client.post(
        "/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    assert resp.status_code == 200
    assert resp.json()["reactants"] == [2, 1]

    resp = client.get("/")
    assert resp.status_code == 200
    assert "reaction-heat" in resp.json()["usage"]
