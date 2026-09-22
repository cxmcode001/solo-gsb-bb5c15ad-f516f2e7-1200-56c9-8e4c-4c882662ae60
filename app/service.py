"""业务编排层：解析、配平、缓存与响应组装。"""

from __future__ import annotations

from typing import Any

from .balancer import (
    BalanceError,
    Species,
    build_balanced_view,
    build_raw_view,
    enumerate_variants,
    parse_side,
    solve_reaction,
)
from .cache import ReactionRecord, build_signature, cache, species_key
from .parser import FormulaError


# 错误码 -> HTTP 状态码
STATUS_CODES = {
    "FORMULA_PARSE_ERROR": 400,
    "DUPLICATE_SPECIES": 409,
    "INVALID_REQUEST": 422,
    "IMPOSSIBLE_BALANCE": 422,
    "NON_UNIQUE_BALANCE": 422,
    "ZERO_COEFFICIENT": 422,
    "INVALID_STOICHIOMETRY": 422,
    "CHARGE_NOT_CONSERVED": 422,
    "VARIANT_SPACE_TOO_LARGE": 422,
}


def error_payload(code: str, message: str, details: dict[str, Any]) -> dict[str, Any]:
    return {"code": code, "message": message, "details": details}


def _parse_all(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[list[Species], list[Species]]:
    try:
        return list(parse_side(reactant_texts, "反应物")), list(
            parse_side(product_texts, "生成物")
        )
    except FormulaError as exc:
        raise BalanceError("FORMULA_PARSE_ERROR", str(exc)) from exc


def _equation(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> str:
    def term(sp: Species) -> str:
        coef = coefficient_map.get(species_key(sp), 1) if coefficient_map else 1
        return ("" if coef == 1 else str(coef)) + sp.normalized

    return " + ".join(map(term, reactants)) + " -> " + " + ".join(map(term, products))


def _normalized_block(
    reactants: list[Species],
    products: list[Species],
    coefficient_map: dict[str, int] | None,
) -> dict[str, Any]:
    return {
        "reactants": [sp.normalized for sp in reactants],
        "products": [sp.normalized for sp in products],
        "equation": _equation(reactants, products, coefficient_map),
    }


def _species_coefficients(
    side: list[Species], coefficient_map: dict[str, int]
) -> list[dict[str, Any]]:
    return [
        {
            "input": sp.raw,
            "normalized": sp.normalized,
            "coefficient": coefficient_map[species_key(sp)],
        }
        for sp in side
    ]


def _coef_map(solved) -> dict[str, int]:
    return {
        species_key(sp): coef
        for sp, coef in zip(solved.all_species, solved.coefficients)
    }


# ------------------------------------------------------------ 解析失败处理
def _text_signature(reactant_texts: list[str], product_texts: list[str]) -> str:
    """解析失败时的退化幂等键：压缩空白后排序。"""
    def norm(texts: list[str]) -> list[str]:
        return sorted("".join(t.split()) for t in texts)

    return f"raw:r={norm(reactant_texts)!r};p={norm(product_texts)!r}"


def _record_parse_error(
    reactant_texts: list[str], product_texts: list[str], exc: BalanceError
) -> tuple[ReactionRecord, bool]:
    signature = _text_signature(reactant_texts, product_texts)
    existing = cache.get_by_signature(signature)
    if existing is not None:
        return existing, True
    record = cache.put(
        signature, [], [], None, [],
        error_payload(exc.code, exc.message, dict(exc.details)),
    )
    return record, False


# ------------------------------------------------------------ 配平主流程
def process_balance(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """配平入口，返回 (http_status, payload)。成功与错误结果都按规范签名缓存。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    if record.error is not None:
        return STATUS_CODES[record.error["code"]], {
            "error": record.error,
            "reaction_id": record.reaction_id,
            "cached": cached,
        }
    return 200, build_balance_payload(record, reactants, products, cached)


def build_balance_payload(
    record: ReactionRecord,
    reactants: list[Species],
    products: list[Species],
    cached: bool,
) -> dict[str, Any]:
    coefficient_map = record.coefficients
    charge_left = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in reactants
    )
    charge_right = sum(
        coefficient_map[species_key(sp)] * sp.formula.charge for sp in products
    )
    return {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(reactants, products, coefficient_map),
        "coefficients": {
            "reactants": _species_coefficients(reactants, coefficient_map),
            "products": _species_coefficients(products, coefficient_map),
            "vector": [
                coefficient_map[species_key(sp)] for sp in (*reactants, *products)
            ],
        },
        "balanced": True,
        "elements": list(record.elements),
        "charge_balance": {
            "reactant_total": charge_left,
            "product_total": charge_right,
            "difference": charge_left - charge_right,
            "balanced": charge_left == charge_right,
        },
    }


# ------------------------------------------------------------ 变体枚举
def process_variants(
    reactant_texts: list[str],
    product_texts: list[str],
    max_coefficient: int,
    limit: int,
) -> tuple[int, dict[str, Any]]:
    """枚举多种配平方案，返回 (http_status, payload)。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        return STATUS_CODES[exc.code], {"error": error_payload(
            exc.code, exc.message, dict(exc.details)
        )}

    if len(reactants) + len(products) > 100:
        exc = BalanceError("INVALID_REQUEST", "物种总数不能超过 100")
        return STATUS_CODES[exc.code], {"error": error_payload(
            exc.code, exc.message, dict(exc.details)
        )}

    try:
        solutions, elements, nullity = enumerate_variants(
            reactants, products, max_coefficient=max_coefficient
        )
    except BalanceError as exc:
        return STATUS_CODES[exc.code], {"error": error_payload(
            exc.code, exc.message, dict(exc.details)
        )}

    total_matched = len(solutions)
    truncated = total_matched > limit
    selected = solutions[:limit]

    variants = [
        _variant_block(reactants, products, solution.coefficients)
        for solution in selected
    ]
    return 200, {
        "normalized": _normalized_block(reactants, products, None),
        "max_coefficient": max_coefficient,
        "count": len(variants),
        "truncated": truncated,
        "variants": variants,
        "elements": list(elements),
        "nullity": nullity,
        "degrees_of_freedom": nullity - 1,
    }


def _variant_block(
    reactants: list[Species],
    products: list[Species],
    coefficients: tuple[int, ...],
) -> dict[str, Any]:
    coefficient_map = {
        species_key(sp): coef
        for sp, coef in zip((*reactants, *products), coefficients)
    }
    return {
        "coefficients": {
            "reactants": _species_coefficients(reactants, coefficient_map),
            "products": _species_coefficients(products, coefficient_map),
            "vector": list(coefficients),
        },
        "total": sum(coefficients),
        "equation": _equation(reactants, products, coefficient_map),
    }


# ------------------------------------------------------------ 逐元素分析
def process_analyze(
    reactant_texts: list[str], product_texts: list[str]
) -> tuple[int, dict[str, Any]]:
    """逐元素差异分析。可解析即返回 200（含原始/平衡两个视图）；无法解析返回 400。"""
    try:
        reactants, products = _parse_all(reactant_texts, product_texts)
    except BalanceError as exc:
        record, cached = _record_parse_error(reactant_texts, product_texts, exc)
        return STATUS_CODES[record.error["code"]], {
            "reaction_id": record.reaction_id,
            "cached": cached,
            "error": record.error,
        }

    signature = build_signature(reactants, products)
    record = cache.get_by_signature(signature)
    cached = record is not None
    elements = sorted(
        {el for sp in (*reactants, *products) for el in sp.formula.composition}
    )

    if record is None:
        try:
            solved = solve_reaction(reactant_texts, product_texts)
            record = cache.put(
                signature,
                list(solved.reactants),
                list(solved.products),
                _coef_map(solved),
                list(solved.elements),
                None,
            )
        except BalanceError as exc:
            record = cache.put(
                signature,
                reactants,
                products,
                None,
                list(exc.details.get("elements", [])),
                error_payload(exc.code, exc.message, dict(exc.details)),
            )

    raw = build_raw_view(reactants, products, elements)
    payload: dict[str, Any] = {
        "reaction_id": record.reaction_id,
        "cached": cached,
        "normalized": _normalized_block(
            reactants, products, None if record.error else record.coefficients
        ),
        "raw": _view(reactants, products, raw, coefficient_map=None, assumed=1),
        "balanced_ok": record.error is None,
        "error": record.error,
    }
    if record.error is None:
        solved = solve_reaction(reactant_texts, product_texts)
        balanced_view = build_balanced_view(solved)
        payload["balanced"] = _view(
            list(solved.reactants),
            list(solved.products),
            balanced_view,
            coefficient_map=record.coefficients,
            assumed=None,
        )
    else:
        payload["balanced"] = None
    return 200, payload


def _view(
    reactants: list[Species],
    products: list[Species],
    view: dict[str, Any],
    *,
    coefficient_map: dict[str, int] | None,
    assumed: int | None,
) -> dict[str, Any]:
    def side_block(side: list[Species]) -> list[dict[str, Any]]:
        block = []
        for sp in side:
            coef = (
                coefficient_map[species_key(sp)] if coefficient_map else (assumed or 1)
            )
            block.append(
                {"input": sp.raw, "normalized": sp.normalized, "coefficient": coef}
            )
        return block

    result: dict[str, Any] = {
        "species": {
            "reactants": side_block(reactants),
            "products": side_block(products),
        },
        "elements": view["elements"],
        "charge": view["charge"],
    }
    if assumed is not None:
        result["coefficients_assumed"] = assumed
    return result


def get_reaction(reaction_id: str) -> dict[str, Any] | None:
    """按 reaction_id 取回缓存结果（幂等查询）。"""
    record = cache.get(reaction_id)
    if record is None:
        return None
    reactants = _record_species(record, "canonical_reactants")
    products = _record_species(record, "canonical_products")
    if record.error is not None:
        return {
            "reaction_id": record.reaction_id,
            "cached": True,
            "error": record.error,
        }
    return build_balance_payload(record, reactants, products, True)


def _record_species(record: ReactionRecord, attr: str) -> list[Species]:
    """缓存里只存了 normalized 文本；查询回显时构造轻量 Species。"""
    from .parser import parse_formula

    return [
        Species(
            raw=item["normalized"],
            normalized=item["normalized"],
            formula=parse_formula(item["normalized"]),
        )
        for item in getattr(record, attr)
    ]
