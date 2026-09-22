# 化学反应计量平衡服务

基于 **Python + FastAPI** 的化学反应方程式配平微服务。客户端提交包含
**元素、括号与电荷**的反应物 / 生成物分子式，服务端：

1. 用递归下降解析器解析化学式（支持三种括号、嵌套、下标乘数、离子电荷、电子 `e-`）；
2. 建立 **元素守恒 + 总电荷守恒** 矩阵（反应物列取正、生成物列取负）；
3. 用 `fractions.Fraction` 做精确 **高斯-若尔当消元（RREF）** 求零空间，
   把唯一的有理数解放大、约分得到 **最小互质正整数系数**；
4. 对同一规范化反应做幂等缓存，重复请求返回相同的 `reaction_id` 与结果。

---

## 1. 一条命令启动

需要已安装 Docker 与 Docker Compose v2：

```bash
docker compose up --build -d
```

服务监听 **宿主机端口 `8000`**（映射容器 8000）。

| 入口 | 地址 |
| --- | --- |
| 健康检查 | <http://localhost:8000/health> |
| Swagger 交互文档 | <http://localhost:8000/docs> |
| ReDoc 文档 | <http://localhost:8000/redoc> |

> 想换宿主机端口：编辑 `docker-compose.yml` 的 `ports`，例如
> `"9000:8000"` 后访问 9000 端口。容器内固定监听 8000。

停止服务：

```bash
docker compose down
```

### 不用 Docker 的本地运行方式

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

运行测试：

```bash
pip install -r requirements-dev.txt
pytest          # 101 个测试
```

---

## 2. API 一览（前缀 `/api/v1`）

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/v1/balance` | 配平反应，返回最小正整数系数与电荷守恒结果 |
| POST | `/api/v1/balance/variants` | 枚举多种互质正整数配平方案（自由变量 ≤ 3） |
| POST | `/api/v1/analyze` | 逐元素差异说明：配平前（单位系数）与配平后两个视图 |
| GET | `/api/v1/reactions/{reaction_id}` | 按幂等 ID 取回之前的结果 |
| GET | `/health` | 健康检查 |

### 2.1 配平 `POST /api/v1/balance`

请求体：

```json
{
  "reactants": ["KMnO4", "HCl"],
  "products": ["KCl", "MnCl2", "Cl2", "H2O"]
}
```

响应 `200`：

```json
{
  "reaction_id": "rxn_7ed134d0c6ec5699",
  "cached": false,
  "normalized": {
    "reactants": ["HCl", "KMnO4"],
    "products": ["KCl", "Cl2", "MnCl2", "H2O"],
    "equation": "16HCl + 2KMnO4 -> 2KCl + 5Cl2 + 2MnCl2 + 8H2O"
  },
  "coefficients": {
    "reactants": [
      {"input": "HCl", "normalized": "HCl", "coefficient": 16},
      {"input": "KMnO4", "normalized": "KMnO4", "coefficient": 2}
    ],
    "products": [
      {"input": "KCl", "normalized": "KCl", "coefficient": 2},
      {"input": "Cl2", "normalized": "Cl2", "coefficient": 5},
      {"input": "MnCl2", "normalized": "MnCl2", "coefficient": 2},
      {"input": "H2O", "normalized": "H2O", "coefficient": 8}
    ],
    "vector": [16, 2, 2, 5, 2, 8]
  },
  "balanced": true,
  "elements": ["Cl", "H", "K", "Mn", "O"],
  "charge_balance": {"reactant_total": 0, "product_total": 0, "difference": 0, "balanced": true}
}
```

`curl`：

```bash
curl -X POST http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["KMnO4","HCl"],"products":["KCl","MnCl2","Cl2","H2O"]}'
```

### 2.2 枚举多种配平方案 `POST /api/v1/balance/variants`

`/balance` 只处理**唯一配平**的反应；当守恒方程存在多个自由度（如
`C + O2 -> CO + CO2`）时，可用本接口枚举给定系数上界内的全部互质正整数解。

请求体在 `reactants` / `products` 之外增加：

| 字段 | 范围 | 说明 |
| --- | --- | --- |
| `max_coefficient` | 2～30 | 所有系数都必须落在 `1..max_coefficient` |
| `limit` | 1～100 | 最多返回的方案数；超出时 `truncated=true` |

算法沿用同一个有理数 RREF：当**自由变量（零空间维数）不超过 3 个**时，
把自由分量作为枚举参数精确计算主元分量，只保留元素与总电荷守恒、系数全正、
最大公约数为 1 的解；按 **系数总和升序、再按规范化物种顺序的向量字典序** 排序。

```bash
curl -X POST http://localhost:8000/api/v1/balance/variants \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["C","O2"],"products":["CO","CO2"],"max_coefficient":5,"limit":100}'
```

响应 `200`（节选）：

```json
{
  "normalized": {
    "reactants": ["C", "O2"],
    "products": ["CO", "CO2"],
    "equation": "C + O2 -> CO + CO2"
  },
  "max_coefficient": 5,
  "count": 2,
  "truncated": false,
  "nullity": 2,
  "degrees_of_freedom": 1,
  "elements": ["C", "O"],
  "variants": [
    {
      "coefficients": {
        "reactants": [
          {"input": "C", "normalized": "C", "coefficient": 3},
          {"input": "O2", "normalized": "O2", "coefficient": 2}
        ],
        "products": [
          {"input": "CO", "normalized": "CO", "coefficient": 2},
          {"input": "CO2", "normalized": "CO2", "coefficient": 1}
        ],
        "vector": [3, 2, 2, 1]
      },
      "total": 8,
      "equation": "3C + 2O2 -> 2CO + CO2"
    },
    {
      "coefficients": {
        "reactants": [
          {"input": "C", "normalized": "C", "coefficient": 4},
          {"input": "O2", "normalized": "O2", "coefficient": 3}
        ],
        "products": [
          {"input": "CO", "normalized": "CO", "coefficient": 2},
          {"input": "CO2", "normalized": "CO2", "coefficient": 2}
        ],
        "vector": [4, 3, 2, 2]
      },
      "total": 11,
      "equation": "4C + 3O2 -> 2CO + 2CO2"
    }
  ]
}
```

说明：

* **唯一配平的反应**（零空间维数 1）也返回单个方案；若其最小系数超过
  `max_coefficient`，则返回空列表、`count=0`、`truncated=false`；
* **自由变量超过 3 个**（零空间维数 ≥ 4）时返回 `422`，错误码
  `VARIANT_SPACE_TOO_LARGE`，`details` 中给出 `nullity` 与 `degrees_of_freedom`；
* `truncated` 仅在符合条件的方案数严格大于 `limit` 时为 `true`；
* `vector` 的分量顺序与 `normalized.reactants + normalized.products` 一致。

### 2.3 逐元素差异 `POST /api/v1/analyze`

返回两个视图：

* `raw`：所有物种按系数 1 统计（即“配平前”的差异，便于解释哪里不平衡）；
* `balanced`：用求得的系数统计，此时每个元素与电荷的 `difference` 均为 0。

以高锰酸根酸性还原半反应为例：

```bash
curl -X POST http://localhost:8000/api/v1/analyze \
  -H 'Content-Type: application/json' \
  -d '{"reactants":["MnO4-","H+","e-"],"products":["Mn2+","H2O"]}'
```

`raw.charge` 为 `{"reactant_total": -1, "product_total": 2, "difference": -3}`，
`balanced` 视图给出 `MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O`，电荷两侧均为 +2。

即使反应**无法唯一配平**，`/analyze` 也返回 `200`，同时带 `raw` 视图与
`error` 字段（`balanced_ok=false`、`balanced=null`），方便客户端展示诊断信息。

### 2.4 幂等与查询

* 幂等键与书写顺序、括号风格、空白以及电荷的等价写法无关，只取决于每一侧
  物种的 **（元素组成多重集, 电荷）**。
* 因此下面两个请求得到 **同一个 `reaction_id`**，第二次 `cached=true`：

```bash
curl -X POST .../api/v1/balance -d '{"reactants":["Ca(OH)2","HCl"],"products":["CaCl2","H2O"]}'
curl -X POST .../api/v1/balance -d '{"reactants":["Ca [OH]2"," HCl"],"products":["CaCl2","H2 O"]}'
```

* 用 `GET /api/v1/reactions/{reaction_id}` 可随时取回缓存结果（LRU，默认 1 万条）。

---

## 3. 化学式书写约定

| 特性 | 写法示例 |
| --- | --- |
| 元素（118 种） | `H`、`Fe`、`O`、`Uuo`（必须大写开头） |
| 下标 | `H2O`、`Al2(SO4)3` |
| 三种括号（等价，可嵌套） | `Ca(OH)2`、`K4[Fe(CN)6]`、`Mg3{Fe(CN)6}2` |
| 显式电荷 | `SO4^2-`、`Fe^3+`、`NH4^+`、`O^{2-}`、`[Cu(NH3)4]^2+` |
| 纯文本电荷 | `Fe3+`、`O2-`、`NH4+`、`MnO4-`、`OH-` |
| 电荷前空格消歧 | `SO4 2-` 等价于 `SO4^2-` |
| 电子 / 正电子 | `e-`、`e`（默认 −1）、`e+` |
| 空白 | 化学式内部空白一律忽略 |

**末尾 `数字±`（无 `^`）消歧规则**（沿用常见化学软件约定）：

* 物种只含 **一个原子/基团** 时，数字是电荷：`Fe3+` → Fe³⁺、`O2-` → O²⁻；
* 含 **多个原子** 时，数字是下标、符号单独表示 ±1：`NH4+` → NH₄⁺、`MnO4-` → MnO₄⁻；
* 数字前有空格时一定是电荷：`SO4 2-` → SO₄²⁻；
* 任何不确定的场景都建议写显式 `^`，无歧义。

---

## 4. 配平原理

设所有物种系数向量为 `x`（未知）。对每种元素写一行守恒方程，再加一行
总电荷守恒；反应物一侧的列取正，生成物一侧的列取负，得到齐次方程组
`A x = 0`。

对整数矩阵的每个元素套上 `Fraction` 做 RREF：

* **零空间维数 = 1**：基向量所有分量非零且同号 → 整体取正，按分母最小公倍数
  放大、再除以最大公约数，得到最小互质正整数系数；
* **维数 = 0**：`IMPOSSIBLE_BALANCE`（元素/电荷守恒本身无解）；
* **维数 > 1**：`NON_UNIQUE_BALANCE`（自由度 = 维数 − 1，常见于缺少电子、H₂O、
  H⁺/OH⁻ 等物种）；
* 维数为 1 但有零分量 → `ZERO_COEFFICIENT`；符号不一致 → `INVALID_STOICHIOMETRY`。

`/balance/variants` 对 **维数 2～4（自由变量 ≤ 3）** 的情形做枚举：RREF 给出
参数化 ``x_主元 = Σ (-rref 元素)·x_自由``，把各自由分量在 1～`max_coefficient`
内取值（最多 30³ = 27000 个组合），用整数除法筛出主元分量也为正整数且在上界内
的向量，再去掉最大公约数 > 1 的倍数解，按 (系数总和, 向量字典序) 排序后截取
`limit`；维数 ≥ 4 时枚举空间过大，直接返回 `VARIANT_SPACE_TOO_LARGE`。

全程使用有理数，不会出现浮点误差，大系数反应（如
`10 K4Fe(CN)6 + 122 KMnO4 + 299 H2SO4 -> …`）也能精确求解。

---

## 5. 错误码

错误响应统一形如 `{"error": {"code": ..., "message": ..., "details": {...}}, "reaction_id": ..., "cached": ...}`。

| code | HTTP | 含义 |
| --- | --- | --- |
| `FORMULA_PARSE_ERROR` | 400 | 化学式无法解析（未知元素、括号不匹配、电荷非法等） |
| `DUPLICATE_SPECIES` | 409 | 同一侧出现重复物种 |
| `INVALID_REQUEST` | 422 | 请求体缺失/为空/物种过多（>100） |
| `IMPOSSIBLE_BALANCE` | 422 | 守恒方程组无解（含电荷无法守恒） |
| `NON_UNIQUE_BALANCE` | 422 | 零空间维数 > 1，配平不唯一（返回 `nullity` 与自由度） |
| `ZERO_COEFFICIENT` | 422 | 唯一解要求某些物种系数为 0 |
| `INVALID_STOICHIOMETRY` | 422 | 唯一解符号不一致，物种可能需要换边 |
| `VARIANT_SPACE_TOO_LARGE` | 422 | `/balance/variants` 零空间自由变量 > 3（返回 `nullity` 与自由度） |
| `NOT_FOUND` | 404 | `reaction_id` 不在缓存中 |

不唯一示例：`C + O2 -> CO + CO2` 缺少约束，返回

```json
{
  "error": {
    "code": "NON_UNIQUE_BALANCE",
    "message": "配平不唯一：零空间维数为 2，存在 1 个自由度……",
    "details": {"nullity": 2, "degrees_of_freedom": 1, "elements": ["C", "O"]}
  }
}
```

---

## 6. 项目结构

```
app/
  parser.py    递归下降化学式解析（元素/括号/电荷/电子）
  balancer.py  守恒矩阵、有理数 RREF、零空间、最小正整数、变体枚举、守恒视图
  cache.py     规范化签名（SHA-256）与 LRU 幂等缓存
  service.py   业务编排与响应组装
  models.py    Pydantic 请求/响应模型
  main.py      FastAPI 路由与异常处理
tests/         解析器、配平核心（含变体枚举）、HTTP API 共 101 个测试
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```

## 7. 已验证的反应示例

```
2H2 + O2 -> 2H2O
4Fe + 3O2 -> 2Fe2O3
C3H8 + 5O2 -> 3CO2 + 4H2O
2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 5Cl2 + 8H2O
3Ca(OH)2 + 2H3PO4 -> Ca3(PO4)2 + 6H2O
(NH4)2Cr2O7 -> Cr2O3 + N2 + 4H2O
Fe^2+ -> Fe^3+ + e^-
MnO4^- + 8H^+ + 5e^- -> Mn^2+ + 4H2O
3Cu + 2NO3^- + 8H^+ -> 3Cu^2+ + 2NO + 4H2O
10K4Fe(CN)6 + 122KMnO4 + 299H2SO4
  -> 162KHSO4 + 5Fe2(SO4)3 + 122MnSO4 + 60HNO3 + 60CO2 + 188H2O
```
