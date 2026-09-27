# 化学方程式配平 API

基于 FastAPI + SymPy 的化学方程式配平与化学计量服务。按元素守恒建立齐次线性方程组，
用有理数精确求解零空间，得到唯一比例的严格正整数系数（约去公因子后返回）；
`/stoichiometry` 在此基础上用 `Decimal`/精确分数完成限量试剂与质量关系计算。

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
  main.py            # FastAPI 应用与 /balance、/stoichiometry 接口
  chemistry.py       # 化学式解析（元素、下标、嵌套圆括号）
  balancer.py        # 基于 SymPy 零空间的配平算法
  stoichiometry.py   # 原子量、Decimal 摩尔质量、限量试剂与质量计算
tests/
  test_balance.py
  test_stoichiometry.py
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
