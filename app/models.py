"""Pydantic 请求/响应模型（同时驱动 OpenAPI 文档）。"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field, field_validator


class ReactionRequest(BaseModel):
    reactants: list[str] = Field(
        ...,
        min_length=1,
        description="反应物化学式，如 [\"KMnO4\", \"HCl\"]；电荷写作 SO4^2-、Fe3+、e-",
        examples=[["H2", "O2"]],
    )
    products: list[str] = Field(
        ...,
        min_length=1,
        description="生成物化学式，如 [\"KCl\", \"MnCl2\", \"Cl2\", \"H2O\"]",
        examples=[["H2O"]],
    )

    @field_validator("reactants", "products")
    @classmethod
    def _strip(cls, values: list[str]) -> list[str]:
        cleaned: list[str] = []
        for value in values:
            if not isinstance(value, str):
                raise ValueError("化学式必须是字符串")
            stripped = value.strip()
            if not stripped:
                raise ValueError("化学式不能为空字符串")
            cleaned.append(stripped)
        return cleaned


class SpeciesCoefficient(BaseModel):
    input: str = Field(description="客户端提交的原始化学式")
    normalized: str = Field(description="规范化后的化学式")
    coefficient: int = Field(description="最小互质正整数系数")


class CoefficientBlock(BaseModel):
    reactants: list[SpeciesCoefficient]
    products: list[SpeciesCoefficient]
    vector: list[int] = Field(description="按 reactants + products 顺序的系数向量")


class NormalizedReaction(BaseModel):
    reactants: list[str]
    products: list[str]
    equation: str = Field(description="规范化的化学方程式文本")


class BalanceResponse(BaseModel):
    reaction_id: str
    cached: bool = Field(description="本次结果是否命中幂等缓存")
    normalized: NormalizedReaction
    coefficients: CoefficientBlock
    balanced: bool
    elements: list[str]
    charge_balance: dict[str, Any]


class ElementRow(BaseModel):
    element: str
    reactant_total: int
    product_total: int
    difference: int = Field(description="反应物合计 - 生成物合计；0 表示守恒")
    balanced: bool


class ChargeRow(BaseModel):
    reactant_total: int
    product_total: int
    difference: int
    balanced: bool


class SideSpecies(BaseModel):
    input: str
    normalized: str
    coefficient: int


class StoichView(BaseModel):
    coefficients_assumed: int | None = Field(
        description="该视图使用的统一假设系数；平衡视图省略此字段"
    )
    species: dict[str, list[SideSpecies]]
    elements: list[ElementRow]
    charge: ChargeRow


class AnalyzeResponse(BaseModel):
    reaction_id: str
    cached: bool
    normalized: NormalizedReaction
    raw: StoichView = Field(description="按单位系数统计的原始差异（配平前）")
    balanced: StoichView | None = Field(
        default=None, description="配平后的守恒视图；无法配平时为 null"
    )
    balanced_ok: bool
    error: dict[str, Any] | None = None


class ErrorResponse(BaseModel):
    error: dict[str, Any]
