import pytest
from fastapi.testclient import TestClient

from app.balancer import BalanceError, balance_equation
from app.chemistry import FormulaError, parse_formula
from app.main import app

client = TestClient(app)


# ---------- 解析器 ----------

@pytest.mark.parametrize(
    "formula, expected",
    [
        ("H2O", {"H": 2, "O": 1}),
        ("O2", {"O": 2}),
        ("C", {"C": 1}),
        ("Fe", {"Fe": 1}),
        ("H99", {"H": 99}),
        ("C60", {"C": 60}),
        ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
        ("Fe2(SO4)3", {"Fe": 2, "S": 3, "O": 12}),
        ("Mg(OH)2", {"Mg": 1, "O": 2, "H": 2}),
        ("(NH4)2SO4", {"N": 2, "H": 8, "S": 1, "O": 4}),
        ("K4(ON(SO3)2)2", {"K": 4, "O": 14, "N": 2, "S": 4}),
        ("((H)2)3", {"H": 6}),
        ("CuSO4(H2O)5", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
    ],
)
def test_parse_valid(formula, expected):
    assert parse_formula(formula) == expected


@pytest.mark.parametrize(
    "formula",
    [
        "",
        "   ",
        "2H",        # 数字开头
        "h2O",       # 小写开头
        "H0",        # 下标 0
        "H01",       # 前导零
        "H100",      # 下标超过 99
        "H2O)",      # 多余的右括号
        "(H2O",      # 左括号未闭合
        "H()",       # 空括号
        "()",        # 空括号
        "H 2 O",     # 含空白
        "H2O$",      # 非法字符
        "H[2]",      # 只支持圆括号
        "H2O1x",     # 小写字母不能单独出现
    ],
)
def test_parse_invalid(formula):
    with pytest.raises(FormulaError):
        parse_formula(formula)


# ---------- 配平核心 ----------

def balance_strings(reactants, products):
    parsed = [parse_formula(f) for f in reactants + products]
    return balance_equation(parsed[: len(reactants)], parsed[len(reactants):])


def test_balance_simple():
    assert balance_strings(["H2", "O2"], ["H2O"]) == ([2, 1], [2])


def test_balance_combustion():
    assert balance_strings(["CH4", "O2"], ["CO2", "H2O"]) == ([1, 2], [1, 2])


def test_balance_nested_parentheses():
    assert balance_strings(
        ["Ca(OH)2", "H3PO4"], ["Ca3(PO4)2", "H2O"]
    ) == ([3, 2], [1, 6])


def test_balance_fraction_reduced():
    # 零空间向量带分数，通分后要约去公因子
    assert balance_strings(["C2H6", "O2"], ["CO2", "H2O"]) == ([2, 7], [4, 6])


def test_balance_single_species_each_side():
    assert balance_strings(["H2"], ["H2"]) == ([1], [1])


def test_balance_four_species_per_side():
    # 构造的 4+4 方程，零空间维数为 1，比例唯一
    reactants = ["A5B6", "C3D7", "E7F8", "G2"]
    products = ["A", "BC", "D2E3", "F3G"]
    assert balance_strings(reactants, products) == ([1, 2, 3, 4], [5, 6, 7, 8])


def test_balance_no_solution():
    with pytest.raises(BalanceError) as exc:
        balance_strings(["H2"], ["O2"])
    assert exc.value.code == "no_solution"


def test_balance_non_unique():
    with pytest.raises(BalanceError) as exc:
        balance_strings(["C", "O2"], ["CO", "CO2"])
    assert exc.value.code == "non_unique_solution"


def test_balance_no_positive_solution():
    with pytest.raises(BalanceError) as exc:
        balance_strings(["H2", "H2O"], ["O2"])
    assert exc.value.code == "no_positive_solution"


# ---------- HTTP API ----------

def test_api_balance_ok():
    resp = client.post(
        "/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["reactants"] == [2, 1]
    assert body["products"] == [2]
    assert body["equation"] == "2 H2 + O2 -> 2 H2O"


def test_api_balance_nested():
    resp = client.post(
        "/balance",
        json={
            "reactants": ["Ca(OH)2", "H3PO4"],
            "products": ["Ca3(PO4)2", "H2O"],
        },
    )
    assert resp.status_code == 200
    assert resp.json()["reactants"] == [3, 2]
    assert resp.json()["products"] == [1, 6]


def test_api_invalid_formula_reports_array_and_index():
    resp = client.post(
        "/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O", "H0"]},
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["array"] == "products"
    assert detail["index"] == 1
    assert detail["formula"] == "H0"


def test_api_invalid_formula_in_reactants():
    resp = client.post(
        "/balance",
        json={"reactants": ["(H2"], "products": ["H2O"]},
    )
    assert resp.status_code == 400
    detail = resp.json()["detail"]
    assert detail["array"] == "reactants"
    assert detail["index"] == 0


def test_api_no_solution():
    resp = client.post(
        "/balance", json={"reactants": ["H2"], "products": ["O2"]}
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "no_solution"


def test_api_non_unique():
    resp = client.post(
        "/balance",
        json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "non_unique_solution"


def test_api_no_positive_solution():
    resp = client.post(
        "/balance",
        json={"reactants": ["H2", "H2O"], "products": ["O2"]},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "no_positive_solution"


@pytest.mark.parametrize(
    "payload",
    [
        {"reactants": [], "products": ["H2O"]},                       # 少于 1 个
        {"reactants": ["H2"], "products": []},
        {"reactants": ["H2"] * 5, "products": ["H2O"]},               # 多于 4 个
        {"reactants": ["H2"], "products": ["H2O"] * 5},
        {"reactants": ["H2"]},                                        # 缺字段
        {"reactants": [2], "products": ["H2O"]},                      # 非字符串
    ],
)
def test_api_schema_validation(payload):
    resp = client.post("/balance", json=payload)
    assert 400 <= resp.status_code < 500


def test_api_root():
    resp = client.get("/")
    assert resp.status_code == 200
