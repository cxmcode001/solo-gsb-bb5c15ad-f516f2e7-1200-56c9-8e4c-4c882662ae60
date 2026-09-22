"""/api/v1/balance/variants 枚举接口测试。"""

from math import gcd

import pytest
from fastapi.testclient import TestClient

from app.balancer import enumerate_variants, parse_side
from app.main import app

client = TestClient(app)


def payload(reactants, products, *, max_coefficient=10, limit=100):
    return {
        "reactants": reactants,
        "products": products,
        "max_coefficient": max_coefficient,
        "limit": limit,
    }


# ------------------------------------------------------------ 核心枚举
def test_enumerate_c_o2_vectors():
    reactants = parse_side(["C", "O2"], "反应物")
    products = parse_side(["CO", "CO2"], "生成物")
    solutions, elements, nullity = enumerate_variants(
        reactants, products, max_coefficient=10
    )
    assert nullity == 2
    assert elements == ("C", "O")
    vectors = [s.coefficients for s in solutions]
    # 最小系数总和方案：3C + 2O2 -> 2CO + CO2
    assert vectors[0] == (3, 2, 2, 1)
    assert solutions[0].total == 8
    # 方案按 (总和, 字典序) 排序
    keys = [(s.total, s.coefficients) for s in solutions]
    assert keys == sorted(keys)


def test_enumerate_all_primitives_in_bounds_and_balanced():
    reactants = parse_side(["H2", "O2", "N2"], "反应物")
    products = parse_side(["H2O", "NH3", "H2O2"], "生成物")
    solutions, _, nullity = enumerate_variants(
        reactants, products, max_coefficient=30
    )
    assert nullity == 3
    assert solutions  # 存在多个方案
    for sol in solutions:
        x = sol.coefficients
        assert all(1 <= c <= 30 for c in x)
        assert gcd(*x) == 1
        assert sol.total == sum(x)
        # 守恒复核（物种顺序：H2,N2,O2 | H2O,H2O2,NH3）
        a, b, c, d, e, f = x  # noqa: E741
        assert 2 * a == 2 * d + 2 * e + 3 * f  # H
        assert 2 * b == f                       # N
        assert 2 * c == d + 2 * e               # O


def test_unique_reaction_returns_single_variant():
    reactants = parse_side(["H2", "O2"], "反应物")
    products = parse_side(["H2O"], "生成物")
    solutions, _, nullity = enumerate_variants(
        reactants, products, max_coefficient=10
    )
    assert nullity == 1
    assert [s.coefficients for s in solutions] == [(2, 1, 2)]


def test_unique_variant_above_bound_is_empty():
    reactants = parse_side(["KMnO4", "HCl"], "反应物")
    products = parse_side(["KCl", "MnCl2", "Cl2", "H2O"], "生成物")
    solutions, _, _ = enumerate_variants(
        reactants, products, max_coefficient=2
    )
    assert solutions == []


def test_nullity_zero_impossible():
    from app.balancer import BalanceError

    reactants = parse_side(["H2O"], "反应物")
    products = parse_side(["Au"], "生成物")
    with pytest.raises(BalanceError) as exc:
        enumerate_variants(reactants, products, max_coefficient=10)
    assert exc.value.code == "IMPOSSIBLE_BALANCE"


def test_too_many_free_variables():
    from app.balancer import BalanceError

    reactants = parse_side(["H2", "O2", "N2", "C"], "反应物")
    products = parse_side(["H2O", "CO", "NO", "CH4"], "生成物")
    with pytest.raises(BalanceError) as exc:
        enumerate_variants(reactants, products, max_coefficient=10)
    assert exc.value.code == "VARIANT_SPACE_TOO_LARGE"
    assert exc.value.details["nullity"] == 4
    assert exc.value.details["max_free_variables"] == 3


def test_nullity_one_zero_coefficient_error():
    from app.balancer import BalanceError

    reactants = parse_side(["H2", "O2", "N2"], "反应物")
    products = parse_side(["H2O"], "生成物")
    with pytest.raises(BalanceError) as exc:
        enumerate_variants(reactants, products, max_coefficient=10)
    assert exc.value.code == "ZERO_COEFFICIENT"
    assert "N2" in exc.value.details["zero_coefficient_species"]


# ------------------------------------------------------------ HTTP 接口
def test_variants_success_structure_and_order():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(["C", "O2"], ["CO", "CO2"], max_coefficient=5),
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["max_coefficient"] == 5
    assert data["nullity"] == 2
    assert data["degrees_of_freedom"] == 1
    assert data["elements"] == ["C", "O"]
    assert data["normalized"]["reactants"] == ["C", "O2"]
    assert data["normalized"]["products"] == ["CO", "CO2"]
    totals = [v["total"] for v in data["variants"]]
    assert totals == sorted(totals)
    first = data["variants"][0]
    assert first["coefficients"]["vector"] == [3, 2, 2, 1]
    assert first["equation"] == "3C + 2O2 -> 2CO + CO2"
    for variant in data["variants"]:
        x = variant["coefficients"]["vector"]
        assert all(1 <= c <= 5 for c in x)
        assert gcd(*x) == 1


def test_variants_limit_truncation():
    full = client.post(
        "/api/v1/balance/variants",
        json=payload(["C", "O2"], ["CO", "CO2"], max_coefficient=10, limit=100),
    ).json()
    total = full["count"]
    assert total == 15
    assert full["truncated"] is False

    capped = client.post(
        "/api/v1/balance/variants",
        json=payload(["C", "O2"], ["CO", "CO2"], max_coefficient=10, limit=3),
    ).json()
    assert capped["count"] == 3
    assert capped["truncated"] is True
    assert capped["variants"] == full["variants"][:3]

    exact = client.post(
        "/api/v1/balance/variants",
        json=payload(["C", "O2"], ["CO", "CO2"], max_coefficient=10, limit=15),
    ).json()
    assert exact["count"] == 15
    assert exact["truncated"] is False


def test_variants_unique_balancing_single_item():
    data = client.post(
        "/api/v1/balance/variants",
        json=payload(["H2", "O2"], ["H2O"]),
    ).json()
    assert data["count"] == 1
    assert data["truncated"] is False
    assert data["nullity"] == 1
    assert data["variants"][0]["equation"] == "2H2 + O2 -> 2H2O"
    assert data["variants"][0]["coefficients"]["vector"] == [2, 1, 2]


def test_variants_empty_when_bound_excludes_unique_solution():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(
            ["KMnO4", "HCl"], ["KCl", "MnCl2", "Cl2", "H2O"],
            max_coefficient=2,
        ),
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["count"] == 0
    assert data["variants"] == []
    assert data["truncated"] is False


def test_variants_space_too_large_422():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(
            ["H2", "O2", "N2", "C"], ["H2O", "CO", "NO", "CH4"],
        ),
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "VARIANT_SPACE_TOO_LARGE"
    assert body["error"]["details"]["nullity"] == 4


def test_variants_impossible_422():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(["H2O"], ["Au"]),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


@pytest.mark.parametrize(
    "max_coefficient,limit",
    [(1, 10), (31, 10), (2, 0), (2, 101)],
)
def test_variants_parameter_ranges(max_coefficient, limit):
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(["H2"], ["H2O"], max_coefficient=max_coefficient, limit=limit),
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"


def test_variants_parse_error_400():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(["Xx"], ["H2O"]),
    )
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "FORMULA_PARSE_ERROR"


def test_variants_duplicate_species_409():
    resp = client.post(
        "/api/v1/balance/variants",
        json=payload(["H2", "H2"], ["H2O"]),
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_variants_charge_conserved():
    data = client.post(
        "/api/v1/balance/variants",
        json=payload(
            ["MnO4-", "H+", "e-"], ["Mn2+", "H2O"],
            max_coefficient=10,
        ),
    ).json()
    assert data["count"] == 1
    # 物种按既有规范键排序（反应物 e- < H^+ < MnO4^-，生成物 H2O < Mn^2+）
    assert data["variants"][0]["equation"] == (
        "5e^- + 8H^+ + MnO4^- -> 4H2O + Mn^2+"
    )
    vector = data["variants"][0]["coefficients"]["vector"]
    assert vector == [5, 8, 1, 4, 1]
