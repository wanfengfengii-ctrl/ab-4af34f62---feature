# 木构文物两芯联合树轮对齐服务

将同一根梁材的两支树芯（A / B）**共同**对入地区年表，避免两支分别寻找
最佳位置后给出互相矛盾的砍伐年代。

- `POST /api/dendro/align`：两芯联合对齐
- `GET /health`：健康检查

## 问题模型

单支树芯在“样芯环 × 年表项”网格上寻找路径：

| 动作 | 含义 |
| --- | --- |
| 匹配 `(i,j)→(i+1,j+1)` | 样芯第 `i` 环对入年表第 `j` 项，要求 `j` 不是缺测下标且 `|sample[i]-chron[j]| ≤ tolerance[i]` |
| 伪环 `(i,j)→(i+1,j)` | 跳过一个**样芯**环（年表中无对应年） |
| 缺环 `(i,j)→(i,j+1)` | 跳过一个**可用**年表项（该年样芯缺轮） |
| 跨缺测 `(i,j)→(i,j+1)` | `j` 为缺测下标：同样推进年表，但免费、不占缺环限额、不计目标 |

约束：

- 两支**末环都必须匹配**，且落在**同一日历年**；
- 每支至少匹配 8 环；
- 首环与末环不允许被跳过，年表窗口前后的多余项不计缺环；
- 缺环、伪环各不超过该支给定限额（0–3）；
- 所有匹配差值不得越过逐环容差。

年表缺测项（请求可选字段 `chronology_gaps`）：

- 一至五个**严格递增**、位于年表范围内的零基下标，只有占位值、不可比对；
- 缺测项**不得与任一样芯环匹配**（包括末环，故终年不能落在缺测项）；
- 缺测项位于某支首末匹配之间时可被该支**跨过**：两支共同对入同一窗口，
  资料缺口由两支联合跨越，推进年表但**不计缺环、不消耗缺环限额、
  不计入三级目标的总跳环**；
- 首匹配之前不允许推进年表，故首匹配之前的缺测项不会被“跨过”；
- 日历年仍按原年表下标连续换算（缺测下标也占用日历年）；
- 省略 `chronology_gaps` 时，请求、响应、裁决与失败行为完全保持兼容。

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
  "chronology_gaps": [10, 14]
}
```

`chronology_gaps` 可省略；给出时须为 1–5 个严格递增且在年表下标范围
`[0, len(chronology)-1]` 内的整数。下标重复、乱序或越界返回可定位的 422。

## 响应

成功（HTTP 200）。省略 `chronology_gaps` 时与既有响应结构完全一致；
启用缺测后，每支映射（含每一份同优见证）额外分列：

- `crossed_chronology_gap_indices` / `crossed_chronology_gap_years`：
  该芯首末匹配之间**实际跨过**的年表缺测下标与对应日历年；
- `crossed_gap_count`：跨过的缺测项数（恒不计入 `missing_count`
  与 `objectives.total_skips`）；
- 原 `missing_chronology_*` 只保留可用年表项上的木材缺环。

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
      "crossed_chronology_gap_indices": [14],
      "crossed_chronology_gap_years": [2004],
      "crossed_gap_count": 1,
      "false_sample_rings": [],
      "missing_count": 1, "false_count": 0, "match_count": 10,
      "end_year": 2009,
      "max_difference": 0, "difference_sum": 0
    },
    "b": { }
  },
  "objectives": {"total_skips": 1, "max_difference": 0, "difference_sum": 0},
  "witnesses": [ /* 按映射序稳定的前两份同优见证，各自分列跨过缺测 */ ],
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
    "reason": "last_ring_unmatchable | insufficient_matchable_rings | skip_limits_exceeded | no_common_end_year | common_end_requires_more_skips",
    "shared_candidate_years": [],
    "cores": {
      "a": { "feasible_within_limits": false, "...": "两支各自的终年候选、限额内外最少跳环数等" },
      "b": { }
    },
    "shared_end_details": [ ],
    "chronology_gaps": null
  }
}
```

启用 `chronology_gaps` 且无解时，证据会把**年表缺测**与**木材跳环**
分开：`evidence.chronology_gaps` 回显缺测下标 / 年份，
`cores.*.min_chronology_gaps_crossed_with_eight_matches` 给出达标时至少
跨过的缺测数，`shared_end_details[].a/b` 另含 `chronology_gaps_crossed`、
`crossed_gap_years`、`crossed_gap_count`，而 `missing` / `false` /
`total_skips` 只统计真实木材缺环与伪环。共同终年只能落在缺测项，或放宽
缺测后仍需超过真实跳环限额时，返回 `no_solution`。

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
