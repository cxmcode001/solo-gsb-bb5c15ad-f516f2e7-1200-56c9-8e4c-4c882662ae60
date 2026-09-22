"""HTTP API 测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_balance_success():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2", "O2"], "products": ["H2O"]}
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["balanced"] is True
    vector = data["coefficients"]["vector"]
    assert vector == [2, 1, 2]
    assert data["normalized"]["equation"] == "2H2 + O2 -> 2H2O"
    assert data["charge_balance"]["balanced"] is True
    assert data["cached"] is False


def test_idempotent_same_normalized_reaction():
    payload_a = {"reactants": ["Ca(OH)2", "HCl"], "products": ["CaCl2", "H2O"]}
    # 等价写法：方括号替代圆括号、内部空白不同
    payload_b = {
        "reactants": ["Ca [OH]2", " HCl"],
        "products": ["CaCl2", "H2 O"],
    }

    r1 = client.post("/api/v1/balance", json=payload_a)
    r2 = client.post("/api/v1/balance", json=payload_b)
    assert r1.status_code == 200 and r2.status_code == 200
    d1, d2 = r1.json(), r2.json()
    assert d1["reaction_id"] == d2["reaction_id"]
    assert d2["cached"] is True
    assert d1["coefficients"]["vector"] == d2["coefficients"]["vector"]


def test_idempotent_order_independent():
    r1 = client.post(
        "/api/v1/balance",
        json={"reactants": ["O2", "H2"], "products": ["H2O"]},
    )
    r2 = client.post(
        "/api/v1/balance",
        json={"reactants": ["H2", "O2"], "products": ["H2O"]},
    )
    assert r1.json()["reaction_id"] == r2.json()["reaction_id"]


def test_get_reaction_by_id():
    created = client.post(
        "/api/v1/balance", json={"reactants": ["Fe", "O2"], "products": ["Fe2O3"]}
    ).json()
    rid = created["reaction_id"]
    resp = client.get(f"/api/v1/reactions/{rid}")
    assert resp.status_code == 200
    assert resp.json()["reaction_id"] == rid
    assert resp.json()["cached"] is True


def test_get_reaction_not_found():
    resp = client.get("/api/v1/reactions/rxn_doesnotexist")
    assert resp.status_code == 404
    assert resp.json()["error"]["code"] == "NOT_FOUND"


def test_balance_formula_parse_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["Xx"], "products": ["H2O"]}
    )
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "FORMULA_PARSE_ERROR"
    assert "reaction_id" in body


def test_balance_non_unique_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]}
    )
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "NON_UNIQUE_BALANCE"
    assert body["error"]["details"]["nullity"] == 2


def test_balance_impossible_error():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2O"], "products": ["Au"]}
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "IMPOSSIBLE_BALANCE"


def test_duplicate_species_conflict():
    resp = client.post(
        "/api/v1/balance", json={"reactants": ["H2", "H2"], "products": ["H2O"]}
    )
    assert resp.status_code == 409
    assert resp.json()["error"]["code"] == "DUPLICATE_SPECIES"


def test_invalid_request_body():
    resp = client.post(
        "/api/v1/balance", json={"reactants": [], "products": ["H2O"]}
    )
    assert resp.status_code == 422
    assert resp.json()["error"]["code"] == "INVALID_REQUEST"

    resp = client.post("/api/v1/balance", json={"reactants": ["H2"]})
    assert resp.status_code == 422


def test_analyze_shows_raw_and_balanced():
    resp = client.post(
        "/api/v1/analyze", json={"reactants": ["H2", "O2"], "products": ["H2O"]}
    )
    assert resp.status_code == 200
    data = resp.json()
    raw_rows = {row["element"]: row for row in data["raw"]["elements"]}
    assert raw_rows["H"]["reactant_total"] == 2
    assert raw_rows["H"]["product_total"] == 2
    assert raw_rows["H"]["difference"] == 0
    assert raw_rows["O"]["reactant_total"] == 2
    assert raw_rows["O"]["product_total"] == 1
    assert raw_rows["O"]["difference"] == 1
    assert data["raw"]["coefficients_assumed"] == 1
    assert data["balanced_ok"] is True
    balanced_rows = {row["element"]: row for row in data["balanced"]["elements"]}
    assert balanced_rows["O"]["reactant_total"] == 2
    assert balanced_rows["O"]["product_total"] == 2
    assert all(row["balanced"] for row in balanced_rows.values())


def test_analyze_non_unique_returns_error_in_body():
    resp = client.post(
        "/api/v1/analyze",
        json={"reactants": ["C", "O2"], "products": ["CO", "CO2"]},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["balanced_ok"] is False
    assert data["balanced"] is None
    assert data["error"]["code"] == "NON_UNIQUE_BALANCE"
    # 原始视图仍然给出差异
    assert {row["element"] for row in data["raw"]["elements"]} == {"C", "O"}


def test_analyze_charge_rows():
    resp = client.post(
        "/api/v1/analyze",
        json={"reactants": ["MnO4-", "H+", "e-"], "products": ["Mn2+", "H2O"]},
    )
    data = resp.json()
    assert data["balanced"]["charge"] == {
        "reactant_total": 2,
        "product_total": 2,
        "difference": 0,
        "balanced": True,
    }


def test_redox_full_reaction():
    resp = client.post(
        "/api/v1/balance",
        json={
            "reactants": ["Cu", "NO3-", "H+"],
            "products": ["Cu2+", "NO", "H2O"],
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    coef = {item["normalized"]: item["coefficient"]
            for item in data["coefficients"]["reactants"]
            + data["coefficients"]["products"]}
    assert coef == {"Cu": 3, "NO3^-": 2, "H^+": 8, "Cu^2+": 3, "NO": 2, "H2O": 4}
