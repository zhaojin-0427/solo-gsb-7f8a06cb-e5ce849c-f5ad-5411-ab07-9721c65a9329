# 化学方程式配平 API

基于 FastAPI + SymPy 的化学方程式配平与化学计量服务。按元素守恒建立齐次线性方程组，
用有理数精确求解零空间，得到唯一比例的严格正整数系数（约去公因子后返回）；
`/stoichiometry` 在此基础上用 `Decimal`/精确分数完成限量试剂与质量关系计算；
`/reaction-heat` 按 Hess 定律用精确有理数线性组合由已知反应热推导目标反应热。

## 运行

```bash
# 构建镜像并启动（监听宿主机 8000 端口）
docker compose up --build

# 后台运行
docker compose up --build -d

# 不依赖 compose 直接运行容器
docker build -t chemistry-api .
docker run --rm -p 8000:8000 chemistry-api
```

服务默认监听 **8000 端口**（容器内由 `uvicorn --port 8000` 指定）。更换访问端口：

- compose：修改 `docker-compose.yml` 中端口映射左侧，如 `"9000:8000"`；
- docker run：`docker run --rm -p 9000:8000 chemistry-api`（宿主机 9000 → 容器 8000）；
- 本地开发：`uvicorn app.main:app --port 9000`。

启动后访问 <http://localhost:8000>：

- 交互式文档（Swagger UI）：<http://localhost:8000/docs>
- 配平接口：`POST http://localhost:8000/balance`
- 化学计量接口：`POST http://localhost:8000/stoichiometry`
- 反应热推导接口：`POST http://localhost:8000/reaction-heat`

服务完全自包含，不依赖任何外部服务。

## 接口说明

### POST /balance

请求体：`reactants`、`products` 各为 1～4 个化学式组成的数组。

```bash
curl -X POST http://localhost:8000/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["Ca(OH)2", "H3PO4"], "products": ["Ca3(PO4)2", "H2O"]}'
```

成功响应（200），系数按输入顺序返回两侧：

```json
{
  "reactants": [3, 2],
  "products": [1, 6],
  "equation": "3 Ca(OH)2 + 2 H3PO4 -> Ca3(PO4)2 + 6 H2O"
}
```

### POST /stoichiometry

在配平的基础上计算限量试剂、目标产物质量与反应物剩余质量。

请求字段：

| 字段                 | 类型           | 说明                                                                 |
| -------------------- | -------------- | -------------------------------------------------------------------- |
| `reactants`          | string[1..4]   | 反应物化学式                                                         |
| `products`           | string[1..4]   | 生成物化学式                                                         |
| `coefficients`       | int[2..8]      | **严格正整数**系数，顺序为「反应物系数 + 产物系数」，须与 `/balance` 的唯一最简结果逐位一致 |
| `reactant_masses_g`  | string[]       | 与 `reactants` **等长**的克质量数组，取值 **0.001～1e9** 的纯十进制字符串（不接受指数、符号、前导零） |
| `target_product`     | int            | 目标产物在 `products` 中的下标（从 0 开始）                          |

```bash
curl -X POST http://localhost:8000/stoichiometry \
  -H 'Content-Type: application/json' \
  -d '{
    "reactants": ["H2", "O2"],
    "products": ["H2O"],
    "coefficients": [2, 1, 2],
    "reactant_masses_g": ["2", "20"],
    "target_product": 0
  }'
```

成功响应（200），所有质量为固定六位小数字符串（`ROUND_HALF_UP`）：

```json
{
  "limiting_reactant_index": 0,
  "target_product_mass_g": "17.872024",
  "remaining_reactant_masses_g": ["0.000000", "4.127976"]
}
```

判定规则：计算各反应物的「物质的量 / 系数」比值，比值最小且唯一者为限量试剂；
并列最小（两种及以上反应物同时耗尽）返回 422。质量计算全程使用
`Decimal`（摩尔质量、物质的量）与精确分数（比值与并列判定），
仅在输出时保留六位小数。

支持的元素与原子量：

| 元素 | H | C | N | O | S | Cl | Na | K | Ca | P |
| ---- | ------ | ------ | ------ | ------ | ----- | ----- | ----- | ------ | ------ | ------ |
| 原子量 | 1.008 | 12.011 | 14.007 | 15.999 | 32.06 | 35.45 | 22.99 | 39.098 | 40.078 | 30.974 |

### POST /reaction-heat

由已知反应热推导目标反应热（Hess 定律）。每个反应用一个 `net` 映射表示：
**化学式原文 → 净系数**（负整数 = 消耗，正整数 = 生成，系数绝对值 1～99）。
若目标反应是已知反应的线性组合，则目标反应热等于各已知反应热按同一权重的
组合。全程使用精确有理数运算，所有数值结果统一为**分母为正的最简 `p/q`
字符串**（整数也写作 `p/1`，如 `2/1`、`-221/1`）。

请求字段：

| 字段                | 类型           | 说明                                                                 |
| ------------------- | -------------- | -------------------------------------------------------------------- |
| `known`             | object[1..8]   | 已知反应数组，每项含 `net` 与 `heat_kj`                               |
| `known[].net`       | object         | 物种净系数映射（1～12 项），键为化学式原文，值为绝对值 1～99 的非零整数 |
| `known[].heat_kj`   | string         | 反应热（kJ），可带正负号、不含指数的有限十进制字符串                  |
| `target.net`        | object         | 目标反应的净系数映射，约定同 `known[].net`                            |

约束：`known` 与 `target` 涉及的不同物质（按化学式原文计）**总计最多 12 种**；
每个反应（含 target）必须**正负两侧均非空**且**满足元素守恒**。

```bash
# 已知：C + O2 -> CO2（-393.5 kJ）、2CO + O2 -> 2CO2（-566 kJ）
# 求：  2C + O2 -> 2CO 的反应热
curl -X POST http://localhost:8000/reaction-heat \
  -H 'Content-Type: application/json' \
  -d '{
    "known": [
      {"net": {"C": -1, "O2": -1, "CO2": 1}, "heat_kj": "-393.5"},
      {"net": {"CO": -2, "O2": -1, "CO2": 2}, "heat_kj": "-566"}
    ],
    "target": {"net": {"C": -2, "O2": -1, "CO": 2}}
  }'
```

成功响应（200）——`weights` 按 `known` 顺序给出（允许负数和零；权重不唯一时
任取一组，同样成功），`heat_kj` 为目标反应热：

```json
{
  "weights": ["2/1", "-1/1"],
  "heat_kj": "-221/1"
}
```

即目标 = 2×known[0] − 1×known[1]，反应热 = 2×(−393.5) − (−566) = −221 kJ。

判定规则（均在精确有理数下进行）：

1. **数据矛盾（422）**：若已知反应存在「净物质变化为零但热量非零」的线性
   组合（违反 Hess 定律，任何推导结果都不唯一），返回
   `inconsistent_known_reactions`，并附该组合的 `weights` 与 `heat_kj`；
2. **不可推导（422）**：否则若目标反应不能表示为已知反应的线性组合，
   返回 `target_not_derivable`；
3. 否则返回一组权重与目标热量。

输入校验错误一律返回 **422 并带字段路径**（如 `known[0].net`、`target.net`、
`known[1].heat_kj`）：

| error                         | 定位字段                          | 含义                                         |
| ----------------------------- | --------------------------------- | -------------------------------------------- |
| `invalid_formula`             | `field` + `formula`               | net 键不是合法化学式                         |
| `invalid_coefficient`         | `field` + `formula` + `value`     | net 系数为 0 或绝对值超过 99                 |
| `invalid_heat`                | `field` + `value`                 | `heat_kj` 不是有限纯十进制字符串             |
| `empty_side`                  | `field` + `side`                  | 反应缺少消耗侧或生成侧                       |
| `unbalanced_reaction`         | `field` + `element_changes`       | 反应不满足元素守恒（附各元素净变化）         |
| `too_many_species`            | `species_count` + `max_species`   | 不同物质总计超过 12 种                       |

请求结构不合法（`known` 非 1～8 项、`net` 为空或超过 12 项、系数为布尔/浮点/
字符串、`heat_kj` 非字符串等）由 Pydantic 返回 422，`detail[].loc` 定位到具体
字段与下标。

### 化学式语法

- 元素：一个大写字母加可选小写字母（如 `H`、`O`、`Na`）
- 下标：1～99，省略时按 1（如 `H2O`、`C60`）
- 圆括号可嵌套并带下标（如 `Ca(OH)2`、`K4(ON(SO3)2)2`）

### 错误响应

所有业务错误都带**字段/数组定位**。

**无效化学式或不支持的元素（400）**——返回所在数组和序号：

```json
{
  "detail": {
    "error": "invalid_formula",
    "array": "products",
    "index": 1,
    "formula": "H100",
    "message": "subscript '100' at position 1 is out of range (1-99)"
  }
}
```

`/stoichiometry` 中出现未配置原子量的元素（如 `Fe`、`Mg`）返回
`unsupported_element`（400），并附 `array`、`index`、`element` 与
`supported_elements`。

**输入不合法（400，`/stoichiometry`）**：

| error                  | `field` / 定位                                  | 含义                             |
| ---------------------- | ----------------------------------------------- | -------------------------------- |
| `invalid_array_length` | `coefficients` / `reactant_masses_g`            | 数组长度不符（质量须与反应物等长） |
| `invalid_coefficient`  | `coefficients` + `index`                        | 系数非正或非整数                 |
| `invalid_mass`         | `reactant_masses_g` + `index`                   | 质量串格式非法或超出 0.001～1e9  |
| `index_out_of_range`   | `target_product` + `index`、`valid_range`       | 目标产物下标越界                 |

**系数与配平结果不匹配（422）**——`coefficient_mismatch` 给出期望值、
提交值以及每一处差异（含所在侧和侧内下标）：

```json
{
  "detail": {
    "error": "coefficient_mismatch",
    "field": "coefficients",
    "expected": [2, 1, 2],
    "given": [2, 1, 1],
    "mismatches": [
      {"index": 2, "side": "products", "side_index": 0, "given": 1, "expected": 2}
    ]
  }
}
```

注意：提交未约简的同比例系数（如 `[4, 2, 4]`）也会被拒绝，必须就是
`/balance` 返回的互质正整数。

**并列限量试剂（422）**：

```json
{
  "detail": {
    "error": "non_unique_limiting_reagent",
    "field": "reactant_masses_g",
    "tied_indices": [0, 1],
    "extent_ratios": ["1.000000000", "1.000000000"]
  }
}
```

**无法配平（422）**——`detail.error` 区分三种情况（两个接口一致）：

| error                  | 含义                               |
| ---------------------- | ---------------------------------- |
| `no_solution`          | 不存在非平凡配平                   |
| `non_unique_solution`  | 系数比例不唯一（零空间维数 > 1）   |
| `no_positive_solution` | 不存在严格为正的系数组合           |

请求结构不合法（数组为空、超过 4 个式子、字段类型错误、布尔值/浮点混入
整数字段等）由 Pydantic 返回 422，`detail[].loc` 定位到具体字段与下标。

## 项目结构

```
app/
  main.py            # FastAPI 应用与 /balance、/stoichiometry、/reaction-heat 接口
  chemistry.py       # 化学式解析（元素、下标、嵌套圆括号）
  balancer.py        # 基于 SymPy 零空间的配平算法
  stoichiometry.py   # 原子量、Decimal 摩尔质量、限量试剂与质量计算
  reaction_heat.py   # Hess 定律反应热推导（精确有理数线性组合）
tests/
  test_balance.py
  test_stoichiometry.py
  test_reaction_heat.py
Dockerfile           # 暴露 8000，含 HEALTHCHECK
docker-compose.yml   # 端口映射 8000:8000，含 healthcheck
requirements.txt
```

## 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
# 运行测试（需额外安装 pytest、httpx）
pytest tests/
```
