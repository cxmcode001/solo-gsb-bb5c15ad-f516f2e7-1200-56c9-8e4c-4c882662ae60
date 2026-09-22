"""化学反应配平核心。

方法：把各元素原子守恒与总电荷守恒写成齐次线性方程组 ``A x = 0``
（反应物列取正、生成物列取负），用 :class:`fractions.Fraction` 做精确
高斯-若尔当消元（RREF）求得零空间：

* 零空间维数为 1 且解向量所有分量同号非零时，把解化为最小互质正整数，
  得到唯一配平；
* 维数为 0 → 不可能守恒；维数 > 1 → 配平不唯一（自由度 ``d-1``）；
* 维数为 1 但存在 0 分量或符号不一致 → 该反应无法成立。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm
from typing import Iterable

from .parser import ParsedFormula, parse_formula


class BalanceError(Exception):
    """无法配平时抛出，携带机器可读的错误码。"""

    def __init__(self, code: str, message: str, **details: object) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


@dataclass(frozen=True)
class Species:
    """一侧中的一个化学物种。"""

    raw: str              # 用户输入（去首尾空白）
    normalized: str       # 规范化化学式
    formula: ParsedFormula


@dataclass(frozen=True)
class SolvedReaction:
    """配平结果。系数按 reactants + products 的顺序给出。"""

    reactants: tuple[Species, ...]
    products: tuple[Species, ...]
    coefficients: tuple[int, ...]
    elements: tuple[str, ...]

    @property
    def all_species(self) -> tuple[Species, ...]:
        return self.reactants + self.products

    @property
    def reactant_coefficients(self) -> tuple[int, ...]:
        return self.coefficients[: len(self.reactants)]

    @property
    def product_coefficients(self) -> tuple[int, ...]:
        return self.coefficients[len(self.reactants) :]


# ---------------------------------------------------------------- 解析一侧
def parse_side(texts: Iterable[str], side_name: str) -> tuple[Species, ...]:
    """解析一侧的化学式，去重并按规范键排序。"""
    formulas: list[Species] = []
    seen_keys: set[tuple] = set()
    for raw in texts:
        formula = parse_formula(raw)
        if formula.key in seen_keys:
            raise BalanceError(
                "DUPLICATE_SPECIES",
                f"{side_name}中存在重复物种：{formula.normalized}",
                species=formula.normalized,
                side="reactants" if side_name == "反应物" else "products",
            )
        seen_keys.add(formula.key)
        formulas.append(Species(raw=raw.strip(), normalized=formula.normalized, formula=formula))
    formulas.sort(key=lambda sp: sp.formula.key)
    return tuple(formulas)


# ---------------------------------------------------------------- 配平入口
def solve_reaction(
    reactant_texts: Iterable[str],
    product_texts: Iterable[str],
) -> SolvedReaction:
    reactants = parse_side(reactant_texts, "反应物")
    products = parse_side(product_texts, "生成物")

    if not reactants:
        raise BalanceError("INVALID_REQUEST", "反应物列表不能为空")
    if not products:
        raise BalanceError("INVALID_REQUEST", "生成物列表不能为空")
    if len(reactants) + len(products) > 100:
        raise BalanceError("INVALID_REQUEST", "物种总数不能超过 100")

    species = reactants + products

    element_set = sorted(
        {el for sp in species for el in sp.formula.composition}
    )
    # 电子没有任何元素，需要确保矩阵至少包含电荷行
    n_cols = len(species)
    matrix: list[list[Fraction]] = []
    for element in element_set:
        row = [Fraction(sp.formula.composition.get(element, 0)) for sp in species]
        # 守恒式（左侧原子数 = 右侧原子数）写成 A x = 0：反应物列取正、
        # 生成物列取负。物理上成立的正系数解在该约定下整体同号。
        for j in range(len(reactants), n_cols):
            row[j] = -row[j]
        matrix.append(row)
    # 电荷守恒行
    charge_row = [Fraction(sp.formula.charge) for sp in species]
    for j in range(len(reactants), n_cols):
        charge_row[j] = -charge_row[j]
    matrix.append(charge_row)

    basis = _null_space(matrix)

    if len(basis) == 0:
        raise BalanceError(
            "IMPOSSIBLE_BALANCE",
            "元素或电荷守恒方程无解：该反应在给定反应物与生成物下无法配平",
            elements=element_set,
        )
    if len(basis) > 1:
        degrees = len(basis) - 1
        raise BalanceError(
            "NON_UNIQUE_BALANCE",
            f"配平不唯一：零空间维数为 {len(basis)}，存在 {degrees} 个自由度。"
            "通常意味着缺少物种（如电子 e-、H2O、H+ / OH-）或两侧物种不独立。",
            nullity=len(basis),
            degrees_of_freedom=degrees,
            elements=element_set,
        )

    vec = basis[0]
    # 合法解要求所有物种非零且整体同号（矩阵列已带方向，见上）
    if any(x == 0 for x in vec):
        zero_species = [
            species[i].normalized for i, x in enumerate(vec) if x == 0
        ]
        raise BalanceError(
            "ZERO_COEFFICIENT",
            "守恒解要求部分物种的系数为 0，该反应无法作为一个整体配平；"
            "请把无关物种从反应式中去掉。",
            zero_coefficient_species=zero_species,
            elements=element_set,
        )
    signs = {1 if x > 0 else -1 for x in vec}
    if len(signs) != 1:
        raise BalanceError(
            "INVALID_STOICHIOMETRY",
            "守恒解中各物种的系数符号不一致：部分物种必须移到另一侧，"
            "或该反应在化学上不成立。",
            elements=element_set,
        )
    if -1 in signs:
        vec = [-x for x in vec]

    ints = _to_positive_integers(vec)

    # 最后做一次整数复核，防止任何意外
    for element in element_set:
        left = sum(
            ints[i] * species[i].formula.composition.get(element, 0)
            for i in range(len(reactants))
        )
        right = sum(
            ints[i] * species[i].formula.composition.get(element, 0)
            for i in range(len(reactants), n_cols)
        )
        if left != right:
            raise BalanceError(  # 理论上不可达
                "IMPOSSIBLE_BALANCE",
                f"内部校验失败：元素 {element} 不守恒",
                elements=element_set,
            )
    charge_left = sum(
        ints[i] * species[i].formula.charge for i in range(len(reactants))
    )
    charge_right = sum(
        ints[i] * species[i].formula.charge
        for i in range(len(reactants), n_cols)
    )
    if charge_left != charge_right:
        raise BalanceError(
            "CHARGE_NOT_CONSERVED",
            f"总电荷不守恒：反应物侧 {charge_left:+d}，生成物侧 {charge_right:+d}。"
            "氧化还原半反应请记得加入电子 e-。",
            reactant_charge=charge_left,
            product_charge=charge_right,
        )

    return SolvedReaction(
        reactants=reactants,
        products=products,
        coefficients=tuple(ints),
        elements=tuple(element_set),
    )


# ---------------------------------------------------------------- 线性代数
def _rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    """对有理数矩阵做高斯-若尔当消元，返回 (RREF, 主元列)。"""
    m = [row[:] for row in matrix]
    if not m:
        return m, []
    rows, cols = len(m), len(m[0])
    pivot_cols: list[int] = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        pivot = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [x - factor * y for x, y in zip(m[i], m[r])]
        pivot_cols.append(c)
        r += 1
    return m, pivot_cols


def _null_space(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    """返回矩阵零空间的一组基（每列对应一个基向量）。"""
    if not matrix:
        n = 0
    else:
        n = len(matrix[0])
    rref, pivots = _rref(matrix)
    free_cols = [c for c in range(n) if c not in pivots]
    basis: list[list[Fraction]] = []
    for free in free_cols:
        vec = [Fraction(0) for _ in range(n)]
        vec[free] = Fraction(1)
        for row_idx, col_idx in enumerate(pivots):
            vec[col_idx] = -rref[row_idx][free]
        basis.append(vec)
    return basis


def _to_positive_integers(vec: list[Fraction]) -> list[int]:
    """把有理向量按最小公倍数放大，再除以最大公约数，得到互质正整数。"""
    denominators = [x.denominator for x in vec]
    scale = lcm(*denominators) if denominators else 1
    scaled = [int(x * scale) for x in vec]
    common = 0
    for value in scaled:
        common = gcd(common, abs(value))
    if common == 0:
        raise BalanceError("IMPOSSIBLE_BALANCE", "得到全零系数，无法配平")
    return [v // common for v in scaled]


# ---------------------------------------------------------------- 守恒视图
def side_totals(
    species: Iterable[Species],
    coefficients: Iterable[int],
    elements: Iterable[str],
) -> dict[str, int]:
    """计算一侧每种元素的原子总数与总电荷（键 ``__charge__``）。"""
    species = list(species)
    coefficients = list(coefficients)
    totals = {el: 0 for el in elements}
    charge = 0
    for sp, coef in zip(species, coefficients):
        for el in elements:
            totals[el] += coef * sp.formula.composition.get(el, 0)
        charge += coef * sp.formula.charge
    totals["__charge__"] = charge
    return totals


def build_balanced_view(solved: SolvedReaction) -> dict:
    """配平后的逐元素守恒表。"""
    r_tot = side_totals(
        solved.reactants, solved.reactant_coefficients, solved.elements
    )
    p_tot = side_totals(
        solved.products, solved.product_coefficients, solved.elements
    )
    rows = [
        {
            "element": el,
            "reactant_total": r_tot[el],
            "product_total": p_tot[el],
            "difference": r_tot[el] - p_tot[el],
            "balanced": r_tot[el] == p_tot[el],
        }
        for el in solved.elements
    ]
    return {
        "elements": rows,
        "charge": {
            "reactant_total": r_tot["__charge__"],
            "product_total": p_tot["__charge__"],
            "difference": r_tot["__charge__"] - p_tot["__charge__"],
            "balanced": r_tot["__charge__"] == p_tot["__charge__"],
        },
        "all_balanced": all(row["balanced"] for row in rows)
        and r_tot["__charge__"] == p_tot["__charge__"],
    }


def build_raw_view(
    reactants: Iterable[Species],
    products: Iterable[Species],
    elements: Iterable[str] | None = None,
) -> dict:
    """按单位系数（用户原始写法，不做配平）计算逐元素差异。"""
    reactants = tuple(reactants)
    products = tuple(products)
    if elements is None:
        elements = sorted(
            {el for sp in (*reactants, *products) for el in sp.formula.composition}
        )
    elements = tuple(elements)
    ones_r = [1] * len(reactants)
    ones_p = [1] * len(products)
    r_tot = side_totals(reactants, ones_r, elements)
    p_tot = side_totals(products, ones_p, elements)
    return {
        "coefficients_assumed": 1,
        "elements": [
            {
                "element": el,
                "reactant_total": r_tot[el],
                "product_total": p_tot[el],
                "difference": r_tot[el] - p_tot[el],
                "balanced": r_tot[el] == p_tot[el],
            }
            for el in elements
        ],
        "charge": {
            "reactant_total": r_tot["__charge__"],
            "product_total": p_tot["__charge__"],
            "difference": r_tot["__charge__"] - p_tot["__charge__"],
            "balanced": r_tot["__charge__"] == p_tot["__charge__"],
        },
    }
