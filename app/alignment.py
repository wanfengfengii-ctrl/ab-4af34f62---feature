"""两芯联合树轮对齐核心算法。

同一根梁材的两支树芯必须落入同一地区年表，并保证两支末环对应同一日历年。

单支树芯在“样芯环 × 年表项”网格上做路径搜索：

* 匹配 (i,j) -> (i+1,j+1)：样芯第 i 环对入年表第 j 项，差值不得越过该环容差；
* 伪环 (i,j) -> (i+1,j)   ：跳过一个样芯环（年表中没有对应年）；
* 缺环 (i,j) -> (i,j+1)   ：跳过一个年表项（该年样芯缺轮）。

路径第一步与最后一步都必须是匹配，因此样芯只是放在年表窗口内，
年表前后多余项不会被误记为缺环；所有跳过都严格位于已匹配环之间。

年表缺测项（chronology_gaps，零基下标）是只有占位值、不可用于比对的年：

* 任一样芯环都不得匹配缺测下标；
* 缺测项位于该芯首末匹配之间时可被“跨过”：推进年表但不计缺环，
  不消耗缺环限额，也不计入任何优化目标的跳环计数；
* 首对匹配之前不允许推进年表，故首匹配不得落在缺测项之后；
  末环必须匹配，故终年本身也不得是缺测项。

日历年仍按原年表下标连续换算（缺测下标也占用日历年）。

一条完整映射也可等价表述为一串严格双增的匹配对
    (0, j0) < (i1, j1) < ... < (n-1, e)
相邻匹配对之间的样芯环即伪环、非缺测年表项即缺环、缺测年表项即跨过。

优化按字典序最小化：
    1. 两芯缺环与伪环总数；
    2. 全部匹配差值的最大绝对值；
    3. 全部匹配差值绝对值总和。

实现分两阶段：
    A. 帕累托动态规划求每个终年的最优指标（同指标向量只留一条代表，
       指标值与映射数量无关，剪枝不影响最优值）；
    B. 固定最优目标后，用记忆化 DFS 按映射序精确枚举前两份同优见证，
       因此“唯一 / 歧义”的判定是精确的。
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Optional


# --------------------------------------------------------------------------- #
# 阶段 A：帕累托动态规划（求指标，不求全部映射）
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class _Rec:
    miss: int
    fls: int
    mmax: int
    dsum: int
    last: str                 # 'M' / 'F' / 'X'
    j0: int = -1              # 首匹配年表下标（缺测证据用）


def _v(r: _Rec) -> tuple[int, int, int, int]:
    return (r.miss, r.fls, r.mmax, r.dsum)


def _add(bucket: list[_Rec], r: _Rec) -> None:
    """并入同一网格状态：严格支配者丢弃；指标向量完全相同只留一条代表。

    指标相同时保留首匹配下标更大（年表缺测跨过更少）的代表，供无解证据
    汇报“最少跨过缺测”的路径；指标值本身与映射数量无关，不影响裁决。
    """
    rv = _v(r)
    for k, q in enumerate(bucket):
        qv = _v(q)
        if qv == rv:
            if r.j0 > q.j0:
                bucket[k] = r
            return
        if qv[0] <= rv[0] and qv[1] <= rv[1] and qv[2] <= rv[2] and qv[3] <= rv[3]:
            return
    bucket.append(r)


def pareto_fronts(
    sample: list[int],
    chron: list[int],
    tol: list[int],
    miss_cap: int,
    false_cap: int,
    gaps: tuple[int, ...] = (),
) -> list[list[list[_Rec]]]:
    """返回 dp[i][j]：处理完前 i 个样芯环、年表推进到 j 的帕累托前沿。

    gaps 中的年表下标不可匹配；首末匹配之间跨过缺测项的推进免费
    （不计缺环、不受限额约束）。
    """
    n, m = len(sample), len(chron)
    gap_set = set(gaps)
    dp: list[list[list[_Rec]]] = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    # 首环必须匹配，不允许开头跳过：缺测下标不可匹配，首匹配之后的
    # 缺测项才允许被跨过。
    for s in range(m):
        if s in gap_set:
            continue
        d = abs(sample[0] - chron[s])
        if d <= tol[0]:
            _add(dp[1][s + 1], _Rec(0, 0, d, d, "M", s))

    for i in range(1, n):
        for j in range(m + 1):
            for r in tuple(dp[i][j]):
                # 伪环：跳过样芯第 i 环
                if r.fls < false_cap:
                    _add(dp[i + 1][j],
                         _Rec(r.miss, r.fls + 1, r.mmax, r.dsum, "F", r.j0))
                if j < m:
                    if j in gap_set:
                        # 缺测项：只能跨过，免费且不占缺环限额。
                        _add(dp[i][j + 1],
                             _Rec(r.miss, r.fls, r.mmax, r.dsum, "X", r.j0))
                    else:
                        d = abs(sample[i] - chron[j])
                        # 缺环：跳过一个可用年表项
                        if r.miss < miss_cap:
                            _add(dp[i][j + 1],
                                 _Rec(r.miss + 1, r.fls, r.mmax, r.dsum,
                                      "X", r.j0))
                        # 匹配
                        if d <= tol[i]:
                            _add(dp[i + 1][j + 1],
                                 _Rec(r.miss, r.fls, max(r.mmax, d),
                                      r.dsum + d, "M", r.j0))
    return dp


def end_records(dp: list[list[list[_Rec]]], n: int,
                min_matches: int) -> dict[int, list[_Rec]]:
    """终年表下标 -> 末步为匹配且匹配数达标的记录。"""
    out: dict[int, list[_Rec]] = {}
    for j, bucket in enumerate(dp[n]):
        rs = [r for r in bucket if r.last == "M" and n - r.fls >= min_matches]
        if rs:
            out[j - 1] = rs
    return out


def core_front(
    sample: list[int], chron: list[int], tol: list[int],
    miss_cap: int, false_cap: int, min_matches: int = 8,
    gaps: tuple[int, ...] = (),
) -> dict[int, list[_Rec]]:
    dp = pareto_fronts(sample, chron, tol, miss_cap, false_cap, gaps)
    return end_records(dp, len(sample), min_matches)


# --------------------------------------------------------------------------- #
# 阶段 B：固定目标的精确见证枚举
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Mapping:
    pairs: tuple[tuple[int, int, int], ...]   # (样芯下标, 年表下标, 差值)
    missing: tuple[int, ...]                  # 跨过的可用年表项（木材缺环）
    false: tuple[int, ...]
    crossed: tuple[int, ...] = ()             # 跨过的年表缺测下标（资料缺口）

    @property
    def miss(self) -> int:
        return len(self.missing)

    @property
    def fls(self) -> int:
        return len(self.false)

    @property
    def skips(self) -> int:
        # 缺测下标是资料缺口，不是木材缺环，不计跳环。
        return self.miss + self.fls

    @property
    def mmax(self) -> int:
        return max(d for _, _, d in self.pairs)

    @property
    def dsum(self) -> int:
        return sum(d for _, _, d in self.pairs)

    @property
    def sig(self) -> tuple[tuple[int, int], ...]:
        return tuple((i, j) for i, j, _ in self.pairs)


class _Enumerator:
    """在固定终年、跳环总数、差值上限 / 差值总和目标下枚举映射。

    映射即一串严格双增的匹配对，首对样芯下标为 0、末对样芯下标为 n-1
    且年表下标为 end。按匹配对字典序（即“映射序”）产出。

    相邻匹配对之间：样芯环为伪环；年表下标中的缺测项免费跨过
    （不计 skips），其余年表项才计缺环（计入 skips）。
    """

    def __init__(self, sample, chron, tol, end, skips, mmax_cap,
                 miss_cap, false_cap, gaps=()):
        self.s = sample
        self.c = chron
        self.t = tol
        self.n = len(sample)
        self.end = end
        self.skips = skips
        self.M = mmax_cap
        self.miss_cap = miss_cap
        self.false_cap = false_cap
        self.gaps = tuple(gaps)
        self.gaps_set = frozenset(self.gaps)

    def _gap_count(self, lo: int, hi: int) -> int:
        """区间 [lo, hi) 内的缺测下标数量。"""
        if not self.gaps or lo >= hi:
            return 0
        return bisect_left(self.gaps, hi) - bisect_left(self.gaps, lo)

    def _gap_list(self, lo: int, hi: int) -> list[int]:
        if not self.gaps or lo >= hi:
            return []
        a = bisect_left(self.gaps, lo)
        b = bisect_left(self.gaps, hi)
        return list(self.gaps[a:b])

    @lru_cache(maxsize=None)
    def _feasible(self, i: int, j: int, sk_left: int,
                  dsum_left: Optional[int]) -> bool:
        """从“下一对样芯下标 >= i、年表下标 >= j”出发能否合法收尾。"""
        if sk_left < 0 or i > self.n - 1 or j > self.end:
            return False
        if i == self.n - 1:
            # 末环必须匹配 end（end 保证不是缺测下标）。j..end-1 只能全部
            # 以缺环（可用槽，消耗 sk_left）或跨过缺测（免费）通过，
            # 末环前不允许再有伪环。
            if j > self.end:
                return False
            wood_slots = (self.end - j) - self._gap_count(j, self.end)
            if wood_slots != sk_left:
                return False
            d = abs(self.s[i] - self.c[self.end])
            return d <= self.t[i] and d <= self.M and (
                dsum_left is None or d == dsum_left)
        rings_left = self.n - i
        span = self.end - j
        # [j, end) 内可用（非缺测）年表槽数 W 与缺测槽数；缺测槽免费跨过。
        gspan = self._gap_count(j, self.end)
        wspan = span - gspan
        # 尾部几何：剩余 R 个样芯环中 f 个作伪环 -> p=R-f 对匹配；
        # 末对占 end，其余 p-1 对占可用槽，木材缺环数
        #   m = W - (p-1) = W - R + 1 + f ≥ 0，
        # 且尾部总跳环必须恰为 sk_left：f + m = sk_left。
        # 故 f 被唯一确定为 (sk_left - W + R - 1)/2，必须为区间内非负整数。
        num = sk_left - wspan + rings_left - 1
        if num < 0 or num % 2 != 0:
            return False
        f_star = num // 2
        # 伪环数有界：0 ≤ f ≤ R-1、f ≤ sk_left，
        # 且木材缺环 m = W - R + 1 + f ≥ 0。
        if not 0 <= f_star <= rings_left - 1 or f_star > sk_left:
            return False
        if wspan - rings_left + 1 + f_star < 0:
            return False
        # 枚举下一对 (ip, jp)：跳过样芯 i..ip-1（伪环）与年表 j..jp-1，
        # 其中缺测下标免费，其余为缺环。
        max_ip = min(self.n - 2, i + sk_left)
        for ip in range(i, max_ip + 1):
            f = ip - i
            rem = sk_left - f
            # jp-j 个年表槽里至多 rem 个可用缺环，其余必须是缺测；
            # 用 (j, end) 内缺测数给一个宽松上界，实际计数在循环内复核。
            jp_hi = min(self.end - 1, j + rem + self._gap_count(j, self.end))
            for jp in range(j, jp_hi + 1):
                if jp in self.gaps_set:
                    continue
                x = (jp - j) - self._gap_count(j, jp)
                used = f + x
                if used > sk_left:
                    continue
                d = abs(self.s[ip] - self.c[jp])
                if d > self.t[ip] or d > self.M:
                    continue
                ndl = None if dsum_left is None else dsum_left - d
                if ndl is not None and ndl < 0:
                    continue
                # 伪环与缺环的分顶限额只在整体计数时校验：
                # 累计伪环 / 缺环数无法仅由 sk_left 表达，保守放行，
                # 由产出端重建映射后复核（限额 ≤3，额外搜索量很小）。
                if self._feasible(ip + 1, jp + 1, sk_left - used, ndl):
                    return True
        return False

    def dfs_yield(self, i, j, sk_left, dsum_left, pairs,
                  missing, false, crossed):
        """按映射序惰性产出；dsum_left 给定时只产出差值和恰为目标值的映射。"""
        if i == self.n - 1:
            # 末环必须匹配 end；j..end-1 全部按缺环（可用槽）或跨过缺测
            # 消耗，且可用槽数必须恰为剩余 sk_left（末环前不能再跳伪环）。
            wood_slots = (self.end - j) - self._gap_count(j, self.end)
            d = abs(self.s[i] - self.c[self.end])
            if (j <= self.end and wood_slots == sk_left
                    and d <= self.t[i] and d <= self.M
                    and (dsum_left is None or d == dsum_left)):
                gap_tail = self._gap_list(j, self.end)
                miss_tail = [q for q in range(j, self.end)
                             if q not in self.gaps_set]
                pairs.append((i, self.end, d))
                mp = Mapping(
                    tuple(pairs),
                    tuple(sorted(missing + miss_tail)),
                    tuple(sorted(false)),
                    tuple(sorted(crossed + gap_tail)))
                if (mp.miss <= self.miss_cap and mp.fls <= self.false_cap
                        and len(pairs) >= 8):
                    yield mp
                pairs.pop()
            return
        max_ip = min(self.n - 2, i + sk_left)
        g_to_end = self._gap_count(j, self.end)
        for ip in range(i, max_ip + 1):
            f_skip = ip - i
            if len(false) + f_skip > self.false_cap:
                continue
            # jp-j 槽位 = 缺环（可用项）+ 缺测跨过；可用缺环至多 sk_left-f_skip。
            jp_hi = min(self.end - 1, j + (sk_left - f_skip) + g_to_end)
            for jp in range(j, jp_hi + 1):
                if jp in self.gaps_set:
                    continue
                x_skip = (jp - j) - self._gap_count(j, jp)
                used = f_skip + x_skip
                if used > sk_left:
                    continue
                if len(missing) + x_skip > self.miss_cap:
                    continue
                d = abs(self.s[ip] - self.c[jp])
                if d > self.t[ip] or d > self.M:
                    continue
                ndl = None if dsum_left is None else dsum_left - d
                if ndl is not None and ndl < 0:
                    continue
                if not self._feasible(ip + 1, jp + 1, sk_left - used, ndl):
                    continue
                gap_here = self._gap_list(j, jp)
                miss_here = [q for q in range(j, jp) if q not in self.gaps_set]
                pairs.append((ip, jp, d))
                yield from self.dfs_yield(
                    ip + 1, jp + 1, sk_left - used, ndl, pairs,
                    missing + miss_here,
                    false + list(range(i, ip)),
                    crossed + gap_here)
                pairs.pop()

    def iter_mappings(self, dsum_target, limit=None):
        """从首环对起按映射序惰性产出（至多 limit 份）。"""
        emitted = 0
        for j0 in range(0, self.end + 1):
            if j0 in self.gaps_set:
                continue
            d = abs(self.s[0] - self.c[j0])
            if d > self.t[0] or d > self.M:
                continue
            ndl = None if dsum_target is None else dsum_target - d
            if ndl is not None and ndl < 0:
                continue
            if not self._feasible(1, j0 + 1, self.skips, ndl):
                continue
            # 首对与下一对之间的缺测下标会在 dfs_yield 中按区间收集。
            for mp in self.dfs_yield(1, j0 + 1, self.skips, ndl,
                                     [(0, j0, d)], [], [], []):
                yield mp
                emitted += 1
                if limit is not None and emitted >= limit:
                    return


def _enumerate_simple(sample, chron, tol, end, skips, mmax_cap, dsum_target,
                      miss_cap, false_cap, limit=None, gaps=()):
    """枚举入口；首对之前不允许跳过，跳环全部发生在匹配对之间。

    返回按映射序排列的生成器；limit 给定时至多产出 limit 份。
    缺测下标只允许在首末匹配之间免费跨过。
    """
    en = _Enumerator(sample, chron, tol, end, skips, mmax_cap,
                     miss_cap, false_cap, gaps)
    return en.iter_mappings(dsum_target, limit)


# --------------------------------------------------------------------------- #
# 联合最优选择
# --------------------------------------------------------------------------- #

def _min_skips(front: dict[int, list[_Rec]]) -> dict[int, int]:
    return {e: min(r.miss + r.fls for r in rs) for e, rs in front.items()}


def _joint_objective(fa: dict[int, list[_Rec]],
                     fb: dict[int, list[_Rec]]) -> Optional[dict[str, Any]]:
    """跨终年求 (总跳环, 最大差值, 差值和) 最优值与可达终年集合。"""
    ka, kb = _min_skips(fa), _min_skips(fb)
    shared = set(ka) & set(kb)
    if not shared:
        return None
    s_total = min(ka[e] + kb[e] for e in shared)
    ends = [e for e in shared if ka[e] + kb[e] == s_total]

    mstar: Optional[int] = None
    sstar: Optional[int] = None
    per_end: dict[int, tuple[int, int]] = {}
    for e in ends:
        pa = [r for r in fa[e] if r.miss + r.fls == ka[e]]
        pb = [r for r in fb[e] if r.miss + r.fls == kb[e]]
        m_e = max(min(r.mmax for r in pa), min(r.mmax for r in pb))
        # 联合最大差值恰为 m_e 时，两侧差值和各自可取的最小值。
        d_a = min(r.dsum for r in pa if r.mmax <= m_e)
        d_b = min(r.dsum for r in pb if r.mmax <= m_e)
        per_end[e] = (m_e, d_a + d_b)
        if mstar is None or (m_e, d_a + d_b) < (mstar, sstar):
            mstar, sstar = m_e, d_a + d_b
    best_ends = [e for e in ends if per_end[e] == (mstar, sstar)]
    return {"skips": s_total, "mmax": mstar, "dsum": sstar,
            "ends": sorted(best_ends), "ka": ka, "kb": kb}


def _enumerate_core(sample, chron, tol, end, skips, mstar, miss_cap, false_cap,
                    dsum_target, limit, gaps=()) -> list[Mapping]:
    return _enumerate_simple(sample, chron, tol, end, skips, mstar,
                             dsum_target, miss_cap, false_cap, limit, gaps)


def find_witnesses(sa, sb, chron, ta, tb, fa, fb,
                   miss_caps, false_caps, opt, gaps=()) -> list[tuple[Mapping, Mapping]]:
    """跨所有最优终年，按映射序取前两份联合见证。

    A 侧按映射序惰性枚举；对每个 A 映射，B 侧只枚举差值和恰为
    (sstar - A 差值和) 的前两份。某终年产出两份即够用于判定歧义，
    但为给出跨终年全局映射序的前两份，每个终年各取至多两份后归并。
    """
    mstar, sstar = opt["mmax"], opt["dsum"]
    pool: list[tuple[Mapping, Mapping]] = []
    for e in opt["ends"]:
        b_cache: dict[int, list[Mapping]] = {}

        def b_for(dsum_need: int) -> list[Mapping]:
            if dsum_need not in b_cache:
                b_cache[dsum_need] = list(_enumerate_core(
                    sb, chron, tb, e, opt["kb"][e], mstar,
                    miss_caps[1], false_caps[1], dsum_need, 2, gaps))
            return b_cache[dsum_need]

        found_here = 0
        # A 侧为惰性生成器，凑满两份配对即停止枚举。
        for ma in _enumerate_core(sa, chron, ta, e, opt["ka"][e], mstar,
                                  miss_caps[0], false_caps[0], None, None,
                                  gaps):
            for mb in b_for(sstar - ma.dsum):
                pool.append((ma, mb))
                found_here += 1
                if found_here == 2:
                    break
            if found_here == 2:
                break
    pool.sort(key=lambda t: (t[0].sig, t[1].sig))
    return pool[:2]


# --------------------------------------------------------------------------- #
# 响应组装
# --------------------------------------------------------------------------- #

def _mapping_payload(mp: Mapping, start_year: int, sample: list[int],
                     chron: list[int], include_gaps: bool = False) -> dict[str, Any]:
    payload = {
        "matched": [
            {
                "sample_ring": si + 1,
                "chronology_index": cj,
                "year": start_year + cj,
                "sample_value": sample[si],
                "chronology_value": chron[cj],
                "difference": d,
            }
            for si, cj, d in mp.pairs
        ],
        # 木材缺环：首末匹配之间跳过的“可用”年表项；缺测下标不在此列。
        "missing_chronology_years": [start_year + cj for cj in mp.missing],
        "missing_chronology_indices": list(mp.missing),
        "false_sample_rings": [si + 1 for si in mp.false],
        "missing_count": mp.miss,
        "false_count": mp.fls,
        "match_count": len(mp.pairs),
        "end_year": start_year + mp.pairs[-1][1],
        "max_difference": mp.mmax,
        "difference_sum": mp.dsum,
    }
    if include_gaps:
        # 资料缺口：首末匹配之间实际跨过的年表缺测下标，不计缺环 / 跳环。
        payload["crossed_chronology_gap_indices"] = list(mp.crossed)
        payload["crossed_chronology_gap_years"] = [
            start_year + cj for cj in mp.crossed]
        payload["crossed_gap_count"] = len(mp.crossed)
    return payload


def _witness_payload(ma: Mapping, mb: Mapping, start_year: int,
                     sa: list[int], sb: list[int],
                     chron: list[int], include_gaps: bool = False) -> dict[str, Any]:
    e = ma.pairs[-1][1]
    return {
        "common_end_year": start_year + e,
        "mappings": {
            "a": _mapping_payload(ma, start_year, sa, chron, include_gaps),
            "b": _mapping_payload(mb, start_year, sb, chron, include_gaps),
        },
    }


# --------------------------------------------------------------------------- #
# 无解证据
# --------------------------------------------------------------------------- #

def _relaxed_info(sample, chron, tol, miss_cap, false_cap,
                  min_matches, gaps=()) -> tuple[dict[int, list[_Rec]], Optional[int]]:
    front = core_front(sample, chron, tol, miss_cap, false_cap,
                       min_matches, gaps)
    k = None
    if front:
        k = min(min(r.miss + r.fls for r in rs) for rs in front.values())
    return front, k


def _best_by_end(front: dict[int, list[_Rec]]) -> dict[int, _Rec]:
    best: dict[int, _Rec] = {}
    for e, rs in front.items():
        best[e] = min(rs, key=lambda r: (r.miss + r.fls, r.mmax, r.dsum))
    return best


def _crossed_gaps(j0: int, end: int, gaps: tuple[int, ...]) -> list[int]:
    """首匹配 j0 与终年 end 之间被跨过的缺测下标（均为可用匹配时端点不在缺测）。"""
    a = bisect_right(gaps, j0)
    b = bisect_left(gaps, end)
    return list(gaps[a:b])


def _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                         caps, fa, fb, gaps=()) -> dict[str, Any]:
    (miss_a, false_a), (miss_b, false_b) = caps
    n = (len(sa), len(sb))
    m = len(chron)
    gaps_t = tuple(gaps)
    gap_set = frozenset(gaps_t)

    rel8_a, k8_a = _relaxed_info(sa, chron, ta, m, n[0] - 8, 8, gaps_t)
    rel8_b, k8_b = _relaxed_info(sb, chron, tb, m, n[1] - 8, 8, gaps_t)
    rel1_a, k1_a = _relaxed_info(sa, chron, ta, m, n[0] - 1, 1, gaps_t)
    rel1_b, k1_b = _relaxed_info(sb, chron, tb, m, n[1] - 1, 1, gaps_t)

    # 缺测下标只有占位值，末环不允许落在其上。
    last_a = [j for j in range(m)
              if j not in gap_set and abs(sa[-1] - chron[j]) <= ta[-1]]
    last_b = [j for j in range(m)
              if j not in gap_set and abs(sb[-1] - chron[j]) <= tb[-1]]

    shared_relaxed = sorted(set(rel1_a) & set(rel1_b))

    if k1_a is None or k1_b is None:
        reason = "last_ring_unmatchable"
    elif k8_a is None or k8_b is None:
        reason = "insufficient_matchable_rings"
    elif not fa or not fb:
        reason = "skip_limits_exceeded"
    elif not shared_relaxed:
        reason = "no_common_end_year"
    else:
        reason = "common_end_requires_more_skips"

    def core_ev(rel8, rel1, k8, last, strict, miss_cap, false_cap, nn):
        ev = {
            "feasible_within_limits": bool(strict),
            "last_ring_candidate_years": [start_year + j for j in last],
            "reachable_end_years_within_limits": [
                start_year + e for e in sorted(strict)
            ],
            "min_skips_within_limits": (
                min(min(r.miss + r.fls for r in rs)
                    for rs in strict.values()) if strict else None
            ),
            "min_skips_with_eight_matches": k8,
            "min_skips_with_any_match": (
                min(min(r.miss + r.fls for r in rs)
                    for rs in rel1.values()) if rel1 else None
            ),
            "limits": {"missing": miss_cap, "false": false_cap},
        }
        if gaps_t:
            # 资料缺口与木材跳环分列：给出 8 环达标时至少要跨过的缺测项数。
            ev["min_chronology_gaps_crossed_with_eight_matches"] = (
                min(len(_crossed_gaps(r.j0, end, gaps_t))
                    for end, rs in rel8.items() for r in rs)
                if rel8 else None
            )
        return ev

    ev_a = core_ev(rel8_a, rel1_a, k8_a, last_a, fa, miss_a, false_a, n[0])
    ev_b = core_ev(rel8_b, rel1_b, k8_b, last_b, fb, miss_b, false_b, n[1])

    def best_records(rel8, rel1):
        """每个终年一份代表记录：先跳环，再跨过缺测数，再差值指标。"""
        best: dict[int, _Rec] = {}
        for pool in (rel8, rel1):
            for e, rs in pool.items():
                r = min(rs, key=lambda q: (
                    q.miss + q.fls,
                    len(_crossed_gaps(q.j0, e, gaps_t)),
                    q.mmax, q.dsum))
                if e not in best or (
                    r.miss + r.fls,
                    len(_crossed_gaps(r.j0, e, gaps_t)),
                    r.mmax, r.dsum
                ) < (
                    best[e].miss + best[e].fls,
                    len(_crossed_gaps(best[e].j0, e, gaps_t)),
                    best[e].mmax, best[e].dsum,
                ):
                    best[e] = r
        return best

    if gaps_t:
        best8 = (best_records(rel8_a, {}), best_records(rel8_b, {}))
        best1 = (best_records({}, rel1_a), best_records({}, rel1_b))
    else:
        best8 = (_best_by_end(rel8_a), _best_by_end(rel8_b))
        best1 = (_best_by_end(rel1_a), _best_by_end(rel1_b))

    shared_details = []
    for j in shared_relaxed:
        entry: dict[str, Any] = {
            "year": start_year + j, "chronology_index": j,
        }
        for idx, (miss_cap, false_cap) in enumerate(
            ((miss_a, false_a), (miss_b, false_b))
        ):
            r = best8[idx].get(j) or best1[idx].get(j)
            if r is None:
                entry["ab"[idx]] = None
                continue
            detail = {
                "missing": r.miss,
                "false": r.fls,
                "total_skips": r.miss + r.fls,
                "match_count": n[idx] - r.fls,
                "max_difference": r.mmax,
                "difference_sum": r.dsum,
                "within_limits": (
                    r.miss <= miss_cap and r.fls <= false_cap
                    and n[idx] - r.fls >= 8
                ),
            }
            if gaps_t:
                # 年表缺测（资料缺口）与木材缺环 / 伪环分开说明。
                crossed = _crossed_gaps(r.j0, j, gaps_t)
                detail["chronology_gaps_crossed"] = crossed
                detail["crossed_gap_years"] = [start_year + g for g in crossed]
                detail["crossed_gap_count"] = len(crossed)
            entry["ab"[idx]] = detail
        shared_details.append(entry)

    evidence: dict[str, Any] = {
        "reason": reason,
        "shared_candidate_years": [start_year + j for j in shared_relaxed],
        "cores": {"a": ev_a, "b": ev_b},
        "shared_end_details": shared_details,
    }
    if gaps_t:
        evidence["chronology_gaps"] = {
            "indices": list(gaps_t),
            "years": [start_year + g for g in gaps_t],
        }

    return {
        "status": "no_solution",
        "common_end_year": None,
        "mappings": None,
        "objectives": None,
        "witnesses": [],
        "evidence": evidence,
    }


# --------------------------------------------------------------------------- #
# 顶层入口
# --------------------------------------------------------------------------- #

def align_joint(payload: dict[str, Any]) -> dict[str, Any]:
    chron: list[int] = payload["chronology"]
    start_year: int = payload["start_year"]
    sa, sb = payload["sample_a"], payload["sample_b"]
    ta, tb = payload["tolerance_a"], payload["tolerance_b"]
    miss_a, false_a = payload["max_missing_a"], payload["max_false_a"]
    miss_b, false_b = payload["max_missing_b"], payload["max_false_b"]
    caps = ((miss_a, false_a), (miss_b, false_b))
    gaps_t = tuple(payload.get("chronology_gaps") or ())
    include_gaps = bool(gaps_t)

    fa = core_front(sa, chron, ta, miss_a, false_a, 8, gaps_t)
    fb = core_front(sb, chron, tb, miss_b, false_b, 8, gaps_t)

    opt = _joint_objective(fa, fb)
    if opt is None:
        return _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                                    caps, fa, fb, gaps_t)

    witnesses = find_witnesses(
        sa, sb, chron, ta, tb, fa, fb,
        (miss_a, miss_b), (false_a, false_b), opt, gaps_t,
    )
    if not witnesses:
        # 理论不可达（最优值必有映射支撑），保守按无解处理。
        return _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                                    caps, fa, fb, gaps_t)

    ma0, mb0 = witnesses[0]
    status = "unique" if len(witnesses) == 1 else "ambiguous"
    primary = _witness_payload(ma0, mb0, start_year, sa, sb, chron,
                               include_gaps)
    return {
        "status": status,
        "common_end_year": primary["common_end_year"],
        "mappings": primary["mappings"],
        "objectives": {
            "total_skips": opt["skips"],
            "max_difference": opt["mmax"],
            "difference_sum": opt["dsum"],
        },
        "witnesses": [
            _witness_payload(ma, mb, start_year, sa, sb, chron, include_gaps)
            for ma, mb in witnesses
        ],
        "evidence": None,
    }
