# 木构文物两芯联合树轮对齐服务

将同一根梁材的两支树芯（A / B）**共同**对入地区年表，避免两支分别寻找
最佳位置后给出互相矛盾的砍伐年代。

- `POST /api/dendro/align`：两芯联合对齐
- `GET /health`：健康检查

## 问题模型

单支树芯在“样芯环 × 年表项”网格上寻找路径：

| 动作 | 含义 |
| --- | --- |
| 匹配 `(i,j)→(i+1,j+1)` | 样芯第 `i` 环对入年表第 `j` 项，要求 `|sample[i]-chron[j]| ≤ tolerance[i]` |
| 伪环 `(i,j)→(i+1,j)` | 跳过一个**样芯**环（年表中无对应年） |
| 缺环 `(i,j)→(i,j+1)` | 跳过一个**年表**项（该年样芯缺轮） |

约束：

- 两支**末环都必须匹配**，且落在**同一日历年**；
- 每支至少匹配 8 环；
- 首环与末环不允许被跳过，年表窗口前后的多余项不计缺环；
- 缺环、伪环各不超过该支给定限额（0–3）；
- 所有匹配差值不得越过逐环容差。

### 年表缺测项（chronology_gaps，可选）

数字化时少量年份只有占位值、不可用于比对。请求可选用 `chronology_gaps`
给出 **1–5 个**严格递增、位于年表范围内的零基下标。缺测项：

- **不得与任一树芯环匹配**（连候选位置都不是），也不能充当共同终年；
- 位于某支首末匹配之间时可被该支**跨越**（两支共同跨过即“共同跨年”），
  跨越是**资料缺口**而非木材缺环——**不消耗**任一芯的缺环限额，也**不计入**
  第一级目标的总跳环、不影响后两级差值目标；
- 日历年份仍按原年表下标连续换算（`start_year + 下标`），缺测年份照常给出。

成功时每支映射及每份同优见证都分别给出实际跨过的缺测下标 / 年份
（`crossed_chronology_gap_indices` / `crossed_chronology_gap_years`）；
原有缺环列表只保留真实年表项。省略该字段时请求、响应与失败行为完全兼容
旧版（不出现任何 gap 字段）。

若共同终年只能落在缺测项，返回 `no_solution`（原因
`common_end_blocked_by_gaps`）；若跨越缺测后木材缺环 / 伪环仍超过真实限额，
仍按 `skip_limits_exceeded` 等返回。两类情形都在两芯证据中把**年表缺测**
（`chronology_gap_indices` / `crossable_with_eight_matches`）与**木材跳环**
（`missing` / `false` / `total_skips`）分开列示。

优化目标按字典序：

1. 两芯缺环与伪环**总数**最小；
2. 全部匹配差值的**最大绝对值**最小；
3. 全部匹配差值的**绝对值总和**最小。

对全部同优（三个目标都并列）映射判定 `unique` / `ambiguous`，并按
**映射序**（匹配对 `(样芯环, 年表项)` 序列的字典序）稳定返回前两份见证。

### 算法

- **阶段 A（求最优指标）**：帕累托动态规划。四元指标
  `(缺环, 伪环, 最大差值, 差值和)` 沿路径单调不减，严格支配者剪枝，
  指标向量相同只留一条代表（不影响最优值）。
- **阶段 B（枚举见证）**：固定最优终年与目标值后，以带可行性记忆化的
  DFS 按映射序惰性枚举，只取前两份；“唯一 / 歧义”判定精确。

算法实现在 `app/alignment.py`，并有 60 组随机算例与独立穷举参考实现
逐目标值、逐见证对照（`tests/test_alignment.py`）。

## 请求

```json
{
  "chronology": [/* 25–60 个整数生长指数 */],
  "start_year": 1990,
  "sample_a": [/* 10–24 环整数指数 */],
  "sample_b": [/* 10–24 环整数指数 */],
  "tolerance_a": [/* 与 sample_a 等长的非负逐环容差 */],
  "tolerance_b": [/* 与 sample_b 等长的非负逐环容差 */],
  "max_missing_a": 3,
  "max_false_a": 3,
  "max_missing_b": 3,
  "max_false_b": 3,
  "chronology_gaps": [/* 可选：1–5 个严格递增、位于年表范围内的零基下标 */]
}
```

## 响应

成功（HTTP 200）：

```json
{
  "status": "unique | ambiguous",
  "common_end_year": 2009,
  "mappings": {
    "a": {
      "matched": [
        {"sample_ring": 1, "chronology_index": 9, "year": 1999,
         "sample_value": 63, "chronology_value": 63, "difference": 0}
      ],
      "missing_chronology_years": [2000],
      "missing_chronology_indices": [10],
      "false_sample_rings": [],
      "_crossed_chronology_gap_indices": "仅在给出 chronology_gaps 时出现：实际跨过的缺测下标",
      "_crossed_chronology_gap_years": "对应的缺测年份（按下标连续换算）",
      "missing_count": 1, "false_count": 0, "match_count": 10,
      "end_year": 2009,
      "max_difference": 0, "difference_sum": 0
    },
    "b": { }
  },
  "objectives": {"total_skips": 1, "max_difference": 0, "difference_sum": 0},
  "witnesses": [ /* 按映射序稳定的前两份同优见证 */ ],
  "evidence": null
}
```

无解（HTTP 200，业务状态 `no_solution`）：

```json
{
  "status": "no_solution",
  "common_end_year": null,
  "mappings": null,
  "objectives": null,
  "witnesses": [],
  "evidence": {
    "reason": "last_ring_unmatchable | insufficient_matchable_rings | skip_limits_exceeded | no_common_end_year | common_end_requires_more_skips | common_end_blocked_by_gaps",
    "shared_candidate_years": [],
    "cores": {
      "a": { "feasible_within_limits": false, "...": "两支各自的终年候选、限额内外最少跳环数等" },
      "b": { }
    },
    "shared_end_details": [ ]
  }
}
```

输入越界返回 HTTP 422：

```json
{ "status": "validation_error", "message": "请求参数未通过校验", "errors": [ ] }
```

## 运行

```bash
# 宿主机端口可配置（默认 8000）
APP_PORT=9090 docker compose up --build -d
curl http://localhost:9090/health
```

## verify 一次性校验服务

名为 `verify` 的一次性服务在应用**通过健康检查后**依次执行：

1. 代码测试（pytest）；
2. 应用构建检查（字节码编译、应用导入、健康探测）；
3. 两芯联合对齐 API 冒烟（唯一解 / 歧义 / 422 校验错误 / 无解证据）。

并以自身退出码报告成败：

```bash
docker compose --profile verify up --build \
    --abort-on-container-exit --exit-code-from verify
# 成功退出码 0；任一环节失败非 0
```

也可在容器外直接运行（便于本地开发）：

```bash
pip install -r requirements.txt -r requirements-dev.txt
APP_BASE_URL=http://127.0.0.1:8000 sh scripts/verify.sh
```
