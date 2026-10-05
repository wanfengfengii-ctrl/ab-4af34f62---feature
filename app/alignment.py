"""两芯联合树轮对齐核心算法。

同一根梁材的两支树芯必须落入同一地区年表，并保证两支末环对应同一日历年。

单支树芯在“样芯环 × 年表项”网格上做路径搜索：

* 匹配 (i,j) -> (i+1,j+1)：样芯第 i 环对入年表第 j 项，差值不得越过该环容差；
* 伪环 (i,j) -> (i+1,j)   ：跳过一个样芯环（年表中没有对应年）；
* 缺环 (i,j) -> (i,j+1)   ：跳过一个年表项（该年样芯缺轮）。

年表中另有少量“缺测项”（chronology_gaps）：数字化时只有占位值、根本不可
用于比对。缺测项：

* 不得与任一树芯环匹配（连候选位置都不是）；
* 位于首末匹配之间时可以被两支树芯共同跨越，跨越不消耗任一芯的缺环限额，
  也不计入三级目标的总跳环——它是资料缺口，不是木材缺环；
* 不能充当共同终年（末环必须落在真实年表项上）。

日历年仍按原年表下标连续换算（start_year + 下标），缺测年份照常给出。

路径第一步与最后一步都必须是匹配，因此样芯只是放在年表窗口内，
年表前后多余项不会被误记为缺环；所有跳过都严格位于已匹配环之间。

一条完整映射也可等价表述为一串严格双增的匹配对
    (0, j0) < (i1, j1) < ... < (n-1, e)
相邻匹配对之间的样芯环即伪环、年表项即缺环。

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
    start: int = 0            # 首环匹配的年表下标（仅供证据展示，不参与支配比较）


def _v(r: _Rec) -> tuple[int, int, int, int]:
    return (r.miss, r.fls, r.mmax, r.dsum)


def _add(bucket: list[_Rec], r: _Rec) -> None:
    """并入同一网格状态：严格支配者丢弃；指标向量相同只留一条代表。

    指标相同但首匹配下标不同（只影响证据里的窗口缺测列表）时保留更早的
    首匹配，使证据展示的可跨越缺测集合不被剪枝缩小。
    """
    rv = _v(r)
    for k, q in enumerate(bucket):
        qv = _v(q)
        if qv[0] <= rv[0] and qv[1] <= rv[1] and qv[2] <= rv[2] and qv[3] <= rv[3]:
            if qv == rv and r.start < q.start:
                bucket[k] = r
            return
    bucket.append(r)


def _gap_set(gaps: Optional[list[int]]) -> frozenset[int]:
    return frozenset(gaps or ())


def pareto_fronts(
    sample: list[int],
    chron: list[int],
    tol: list[int],
    miss_cap: int,
    false_cap: int,
    gaps: frozenset[int] = frozenset(),
) -> list[list[list[_Rec]]]:
    """返回 dp[i][j]：处理完前 i 个样芯环、年表推进到 j 的帕累托前沿。"""
    n, m = len(sample), len(chron)
    dp: list[list[list[_Rec]]] = [[[] for _ in range(m + 1)] for _ in range(n + 1)]

    # 首环必须匹配，不允许开头跳过；缺测项不具备匹配资格。
    for s in range(m):
        if s in gaps:
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
                         _Rec(r.miss, r.fls + 1, r.mmax, r.dsum, "F", r.start))
                if j < m:
                    if j in gaps:
                        # 缺测项只能被跨越：非缺环、不耗限额、不计跳环。
                        _add(dp[i][j + 1],
                             _Rec(r.miss, r.fls, r.mmax, r.dsum, "X", r.start))
                    else:
                        d = abs(sample[i] - chron[j])
                        # 缺环：跳过真实年表第 j 项
                        if r.miss < miss_cap:
                            _add(dp[i][j + 1],
                                 _Rec(r.miss + 1, r.fls, r.mmax, r.dsum,
                                      "X", r.start))
                        # 匹配
                        if d <= tol[i]:
                            _add(dp[i + 1][j + 1],
                                 _Rec(r.miss, r.fls, max(r.mmax, d), r.dsum + d,
                                      "M", r.start))
    return dp


def end_records(dp: list[list[list[_Rec]]], n: int,
                min_matches: int,
                gaps: frozenset[int] = frozenset()) -> dict[int, list[_Rec]]:
    """终年表下标 -> 末步为匹配且匹配数达标的记录（终年不允许是缺测项）。"""
    out: dict[int, list[_Rec]] = {}
    for j, bucket in enumerate(dp[n]):
        e = j - 1
        if e in gaps:
            continue
        rs = [r for r in bucket if r.last == "M" and n - r.fls >= min_matches]
        if rs:
            out[e] = rs
    return out


def core_front(
    sample: list[int], chron: list[int], tol: list[int],
    miss_cap: int, false_cap: int, min_matches: int = 8,
    gaps: frozenset[int] = frozenset(),
) -> dict[int, list[_Rec]]:
    dp = pareto_fronts(sample, chron, tol, miss_cap, false_cap, gaps)
    return end_records(dp, len(sample), min_matches, gaps)


# --------------------------------------------------------------------------- #
# 阶段 B：固定目标的精确见证枚举
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Mapping:
    pairs: tuple[tuple[int, int, int], ...]   # (样芯下标, 年表下标, 差值)
    missing: tuple[int, ...]
    false: tuple[int, ...]
    gaps_crossed: tuple[int, ...] = ()        # 被跨越的缺测年表下标

    @property
    def miss(self) -> int:
        return len(self.missing)

    @property
    def fls(self) -> int:
        return len(self.false)

    @property
    def skips(self) -> int:
        # 缺测项是资料缺口而非木材缺环，不计入跳环。
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

    年表缺测项不能作为匹配位置；相邻匹配对之间夹着的缺测项属于“被跨越”，
    既不占 missing 也不消耗跳环预算 sk_left。
    """

    def __init__(self, sample, chron, tol, end, skips, mmax_cap,
                 miss_cap, false_cap, gaps=frozenset()):
        self.s = sample
        self.c = chron
        self.t = tol
        self.n = len(sample)
        self.end = end
        self.skips = skips
        self.M = mmax_cap
        self.miss_cap = miss_cap
        self.false_cap = false_cap
        self.gaps = gaps

    def _matchable(self, jp: int) -> bool:
        return jp not in self.gaps

    def _terminal_resolve(self, i: int, j: int, sk_left: int,
                          dsum_left: Optional[int]
                          ) -> Optional[tuple[list[int], list[int],
                                              list[int], int]]:
        """直接收尾：跳过样芯环 i..n-2（伪环），末环 n-1 匹配 end。

        末环与 end 之间允许夹着若干年表项：真实项计入缺环，缺测项计入
        跨越。两类真实跳环数之和必须恰好用完 sk_left（末环之后不再有
        任何跳环）。
        返回 (尾段伪环, 尾段缺环, 尾段跨越缺测, 末环差值)，不可行返回 None。
        """
        if i > self.n - 1 or j > self.end or not self._matchable(self.end):
            return None
        false_tail = list(range(i, self.n - 1))
        tail = list(range(j, self.end))
        miss_tail = [q for q in tail if q not in self.gaps]
        gap_tail = [q for q in tail if q in self.gaps]
        if len(false_tail) + len(miss_tail) != sk_left:
            return None
        d = abs(self.s[self.n - 1] - self.c[self.end])
        if d > self.t[self.n - 1] or d > self.M:
            return None
        if dsum_left is not None and d != dsum_left:
            return None
        return false_tail, miss_tail, gap_tail, d

    @lru_cache(maxsize=None)
    def _feasible(self, i: int, j: int, sk_left: int,
                  dsum_left: Optional[int]) -> bool:
        """从“下一对样芯下标 >= i、年表下标 >= j”出发能否合法收尾。"""
        if sk_left < 0 or i > self.n - 1 or j > self.end:
            return False
        # 区间 [j, end] 内可作缺环 / 匹配的真实年表槽位数（缺测项可自由跨越）。
        span = self.end - j + 1 - sum(1 for g in self.gaps if j <= g <= self.end)
        # 直接收尾（其后样芯环全部伪环、年表真实项全部缺环）。
        if self._terminal_resolve(i, j, sk_left, dsum_left) is not None:
            return True
        if i == self.n - 1:
            return False
        rings_left = self.n - i
        # 剩余 f 个伪环、x 个缺环且 f+x=sk_left：
        # f ≤ rings_left-1（至少末环要匹配）、x ≤ span-1（至少一个年表槽匹配），
        # 有解当且仅当 span-1 ≥ sk_left-(rings_left-1)。
        if span - 1 < sk_left - (rings_left - 1):
            return False
        # 枚举下一对 (ip, jp)：跳过样芯 i..ip-1（伪环）与年表 j..jp-1 中的
        # 真实项（缺环）；区间内的缺测项只算跨越。
        max_ip = min(self.n - 2, i + sk_left)
        for ip in range(i, max_ip + 1):
            f = ip - i
            rem = sk_left - f
            jp = j
            while jp < self.end:
                # 越过 rem 个真实槽后，更远的 jp 只会更差（x 单调不减）。
                x_real = jp - j - sum(
                    1 for g in self.gaps if j <= g < jp)
                if x_real > rem:
                    break
                if self._matchable(jp):
                    used = f + x_real
                    d = abs(self.s[ip] - self.c[jp])
                    if d <= self.t[ip] and d <= self.M:
                        ndl = None if dsum_left is None else dsum_left - d
                        if ndl is None or ndl >= 0:
                            # 伪环与缺环的分顶限额只在整体计数时校验：
                            # 累计伪环 / 缺环数无法仅由 sk_left 表达，保守
                            # 放行，由产出端重建映射后复核（限额 ≤3，额外
                            # 搜索量很小）。
                            if self._feasible(ip + 1, jp + 1,
                                              sk_left - used, ndl):
                                return True
                jp += 1
        return False

    def dfs_yield(self, i, j, sk_left, dsum_left, pairs, missing, false,
                  crossed):
        """按映射序惰性产出；dsum_left 给定时只产出差值和恰为目标值的映射。"""
        if i == self.n - 1:
            # 末环收尾（此时不存在尾段伪环）。
            resolved = self._terminal_resolve(i, j, sk_left, dsum_left)
            if resolved is not None:
                _, miss_tail, gap_tail, d = resolved
                all_missing = tuple(sorted(missing + miss_tail))
                all_crossed = tuple(sorted(crossed + gap_tail))
                if (len(all_missing) <= self.miss_cap
                        and len(false) <= self.false_cap
                        and len(pairs) + 1 >= 8):
                    pairs.append((i, self.end, d))
                    yield Mapping(tuple(pairs), all_missing,
                                  tuple(sorted(false)), all_crossed)
                    pairs.pop()
            return
        max_ip = min(self.n - 2, i + sk_left)
        for ip in range(i, max_ip + 1):
            f_skip = ip - i
            if len(false) + f_skip > self.false_cap:
                continue
            jp = j
            while jp < self.end:
                between = range(j, jp)
                x_skip = sum(1 for q in between if q not in self.gaps)
                # x_skip 随 jp 单调不减；越过跳环预算后更远位置无需再试。
                if x_skip > sk_left - f_skip:
                    break
                if self._matchable(jp):
                    used = f_skip + x_skip
                    if len(missing) + x_skip <= self.miss_cap:
                        d = abs(self.s[ip] - self.c[jp])
                        if d <= self.t[ip] and d <= self.M:
                            ndl = None if dsum_left is None else dsum_left - d
                            if ndl is None or ndl >= 0:
                                if self._feasible(
                                        ip + 1, jp + 1,
                                        sk_left - used, ndl):
                                    pairs.append((ip, jp, d))
                                    yield from self.dfs_yield(
                                        ip + 1, jp + 1, sk_left - used, ndl,
                                        pairs,
                                        missing + [q for q in between
                                                   if q not in self.gaps],
                                        false + list(range(i, ip)),
                                        crossed + [q for q in between
                                                   if q in self.gaps])
                                    pairs.pop()
                jp += 1

        # 直接收尾（映射序上排在所有“下一对仍落在 i..n-2 环”的续接之后）：
        # i..n-2 环全部作伪环，末环匹配 end，尾段真实年表项计缺环。
        resolved = self._terminal_resolve(i, j, sk_left, dsum_left)
        if resolved is not None:
            false_tail, miss_tail, gap_tail, d = resolved
            all_false = tuple(sorted(false + false_tail))
            all_missing = tuple(sorted(missing + miss_tail))
            all_crossed = tuple(sorted(crossed + gap_tail))
            if (len(all_missing) <= self.miss_cap
                    and len(all_false) <= self.false_cap
                    and len(pairs) + 1 >= 8):
                pairs.append((self.n - 1, self.end, d))
                yield Mapping(tuple(pairs), all_missing, all_false,
                              all_crossed)
                pairs.pop()

    def iter_mappings(self, dsum_target, limit=None):
        """从首环对起按映射序惰性产出（至多 limit 份）。"""
        emitted = 0
        for j0 in range(0, self.end + 1):
            if not self._matchable(j0):
                continue
            d = abs(self.s[0] - self.c[j0])
            if d > self.t[0] or d > self.M:
                continue
            ndl = None if dsum_target is None else dsum_target - d
            if ndl is not None and ndl < 0:
                continue
            # 首对之前不允许跳过，故起点之前的缺测项不会被跨越。
            if not self._feasible(1, j0 + 1, self.skips, ndl):
                continue
            for mp in self.dfs_yield(1, j0 + 1, self.skips, ndl,
                                     [(0, j0, d)], [], [], []):
                yield mp
                emitted += 1
                if limit is not None and emitted >= limit:
                    return


def _enumerate_simple(sample, chron, tol, end, skips, mmax_cap, dsum_target,
                      miss_cap, false_cap, limit=None, gaps=frozenset()):
    """枚举入口；首对之前不允许跳过，跳环全部发生在匹配对之间。

    返回按映射序排列的生成器；limit 给定时至多产出 limit 份。
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
                    dsum_target, limit, gaps=frozenset()) -> list[Mapping]:
    return _enumerate_simple(sample, chron, tol, end, skips, mstar,
                             dsum_target, miss_cap, false_cap, limit, gaps)


def find_witnesses(sa, sb, chron, ta, tb, fa, fb,
                   miss_caps, false_caps, opt,
                   gaps=frozenset()) -> list[tuple[Mapping, Mapping]]:
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
                     chron: list[int], gap_fields: bool = False) -> dict[str, Any]:
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
        # 缺环只列真实年表项；缺测项（资料缺口）单列在下方。
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
    if gap_fields:
        # 仅在启用 chronology_gaps 时输出；省略该字段时响应保持旧形状。
        # 实际被本映射跨越的年表缺测下标 / 年份（不耗缺环限额、不计跳环）。
        payload["crossed_chronology_gap_indices"] = list(mp.gaps_crossed)
        payload["crossed_chronology_gap_years"] = [
            start_year + cj for cj in mp.gaps_crossed]
        payload["crossed_gap_count"] = len(mp.gaps_crossed)
    return payload


def _witness_payload(ma: Mapping, mb: Mapping, start_year: int,
                     sa: list[int], sb: list[int],
                     chron: list[int], gap_fields: bool = False) -> dict[str, Any]:
    e = ma.pairs[-1][1]
    return {
        "common_end_year": start_year + e,
        "mappings": {
            "a": _mapping_payload(ma, start_year, sa, chron, gap_fields),
            "b": _mapping_payload(mb, start_year, sb, chron, gap_fields),
        },
    }


# --------------------------------------------------------------------------- #
# 无解证据
# --------------------------------------------------------------------------- #

def _relaxed_info(sample, chron, tol, miss_cap, false_cap,
                  min_matches, gaps=frozenset()
                  ) -> tuple[dict[int, list[_Rec]], Optional[int]]:
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


def _gap_window(r: Optional[_Rec], end: int,
                gaps: frozenset[int]) -> list[int]:
    """某条终年映射的窗口（首末匹配之间）内包含的缺测下标。

    放宽前沿的代表只记录指标向量与首匹配下标；缺测跨越不耗任何代价，
    故窗口内缺测下标即为被跨越集合（缺测项不能充当匹配，也不可能落在
    首匹配之前 / 末匹配之后）。
    """
    if r is None:
        return []
    return sorted(g for g in gaps if r.start <= g <= end)


def _relaxed_starts(sample, chron, tol, false_cap,
                    gaps: frozenset[int]) -> list[list[dict[int, int]]]:
    """放宽缺环限额（真实年表项可全部跳过）后的可达性 DP。

    返回 states[i][j]：{首匹配下标 -> 最少已用伪环数}。只保留同一首匹配下
    伪环更少的代表（伪环越少支配），专供“缺测项能否在 ≥8 环匹配下跨越”
    这类不关心差值的证据查询。
    """
    n, m = len(sample), len(chron)
    states: list[list[dict[int, int]]] = [
        [{} for _ in range(m + 1)] for _ in range(n + 1)]

    def offer(i, j, start, f):
        old = states[i][j].get(start)
        if old is None or f < old:
            states[i][j][start] = f

    for s in range(m):
        if s in gaps:
            continue
        if abs(sample[0] - chron[s]) <= tol[0]:
            offer(1, s + 1, s, 0)

    for i in range(1, n):
        for j in range(m + 1):
            for start, f in tuple(states[i][j].items()):
                # 伪环：放宽后仍受“至少 8 环匹配”折算的伪环上限约束。
                if f + 1 <= false_cap:
                    offer(i + 1, j, start, f + 1)
                if j < m:
                    if j in gaps:
                        offer(i, j + 1, start, f)          # 缺测：自由跨越
                    else:
                        offer(i, j + 1, start, f)          # 真实项：缺环放宽
                        if abs(sample[i] - chron[j]) <= tol[i]:
                            offer(i + 1, j + 1, start, f)  # 匹配
    return states


def _gap_crossable(sample, chron, tol, states, gaps) -> dict[int, bool]:
    """每个缺测下标是否存在窗口覆盖它且至少 8 环匹配的放宽映射。

    states[n-1][e] 为处理完前 n-1 环、年表推进到 e 的状态；末环再匹配真实
    年表项 e 收尾。首匹配 s < g ≤ e 即保证缺测 g 落在窗口内被跨越。
    """
    n = len(sample)
    out = {}
    for g in gaps:
        ok = False
        for e in range(g, len(chron)):
            if e in gaps or abs(sample[n - 1] - chron[e]) > tol[n - 1]:
                continue
            if any(start < g for start in states[n - 1][e]):
                ok = True
                break
        out[g] = ok
    return out


def _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                         caps, fa, fb, gaps=frozenset()) -> dict[str, Any]:
    (miss_a, false_a), (miss_b, false_b) = caps
    n = (len(sa), len(sb))
    m = len(chron)

    rel8_a, k8_a = _relaxed_info(sa, chron, ta, m, n[0] - 8, 8, gaps)
    rel8_b, k8_b = _relaxed_info(sb, chron, tb, m, n[1] - 8, 8, gaps)
    rel1_a, k1_a = _relaxed_info(sa, chron, ta, m, n[0] - 1, 1, gaps)
    rel1_b, k1_b = _relaxed_info(sb, chron, tb, m, n[1] - 1, 1, gaps)

    # 缺测项不具备匹配资格，末环候选只列真实年表项。
    last_a = [j for j in range(m)
              if j not in gaps and abs(sa[-1] - chron[j]) <= ta[-1]]
    last_b = [j for j in range(m)
              if j not in gaps and abs(sb[-1] - chron[j]) <= tb[-1]]

    shared_relaxed = sorted(set(rel1_a) & set(rel1_b))

    # 若取消缺测限制后两芯放宽终年仍能汇合，且汇合终年“全部”落在缺测
    # 下标上，同时启用缺测后再无共享终年，则共同终年只能落在缺测项。
    gap_end_blocked = False
    if gaps and not shared_relaxed:
        ng1_a = _relaxed_info(sa, chron, ta, m, n[0] - 1, 1, frozenset())[0]
        ng1_b = _relaxed_info(sb, chron, tb, m, n[1] - 1, 1, frozenset())[0]
        shared_ng = set(ng1_a) & set(ng1_b)
        gap_end_blocked = bool(shared_ng) and shared_ng <= set(gaps)

    if gap_end_blocked:
        # 两芯放宽后唯一能汇合的终年落在缺测项上：优先于其它原因标出。
        reason = "common_end_blocked_by_gaps"
    elif k1_a is None or k1_b is None:
        reason = "last_ring_unmatchable"
    elif k8_a is None or k8_b is None:
        reason = "insufficient_matchable_rings"
    elif not fa or not fb:
        reason = "skip_limits_exceeded"
    elif not shared_relaxed:
        reason = "no_common_end_year"
    else:
        reason = "common_end_requires_more_skips"

    def core_ev(rel8, rel1, k8, last, strict, miss_cap, false_cap, nn,
                crossable):
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
        if gaps:
            # 每个缺测项：能否在至少 8 环匹配的放宽窗口内被该芯跨越。
            ev["chronology_gaps"] = [
                {
                    "chronology_index": g,
                    "year": start_year + g,
                    "crossable_with_eight_matches": crossable[g],
                }
                for g in sorted(gaps)
            ]
        return ev

    crossable_a = crossable_b = {}
    if gaps:
        states_a = _relaxed_starts(sa, chron, ta, n[0] - 8, gaps)
        states_b = _relaxed_starts(sb, chron, tb, n[1] - 8, gaps)
        crossable_a = _gap_crossable(sa, chron, ta, states_a, gaps)
        crossable_b = _gap_crossable(sb, chron, tb, states_b, gaps)

    ev_a = core_ev(rel8_a, rel1_a, k8_a, last_a, fa, miss_a, false_a, n[0],
                   crossable_a)
    ev_b = core_ev(rel8_b, rel1_b, k8_b, last_b, fb, miss_b, false_b, n[1],
                   crossable_b)

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
            detail: dict[str, Any] = {
                # 木材缺环 / 伪环（真实跳环，受限额约束）
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
            if gaps:
                # 年表缺测（资料缺口，可共同跨越、不耗限额、不计跳环）。
                gap_idx = _gap_window(r, j, gaps)
                detail["chronology_gap_indices"] = gap_idx
                detail["chronology_gap_years"] = [
                    start_year + g for g in gap_idx]
                detail["chronology_gap_count"] = len(gap_idx)
            entry["ab"[idx]] = detail
        shared_details.append(entry)

    evidence: dict[str, Any] = {
        "reason": reason,
        "shared_candidate_years": [start_year + j for j in shared_relaxed],
        "cores": {"a": ev_a, "b": ev_b},
        "shared_end_details": shared_details,
    }
    if gaps:
        evidence["chronology_gaps"] = {
            "indices": sorted(gaps),
            "years": [start_year + g for g in sorted(gaps)],
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
    gaps = _gap_set(payload.get("chronology_gaps"))

    fa = core_front(sa, chron, ta, miss_a, false_a, 8, gaps)
    fb = core_front(sb, chron, tb, miss_b, false_b, 8, gaps)

    opt = _joint_objective(fa, fb)
    if opt is None:
        return _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                                    caps, fa, fb, gaps)

    witnesses = find_witnesses(
        sa, sb, chron, ta, tb, fa, fb,
        (miss_a, miss_b), (false_a, false_b), opt, gaps,
    )
    if not witnesses:
        # 理论不可达（最优值必有映射支撑），保守按无解处理。
        return _no_solution_payload(start_year, chron, sa, sb, ta, tb,
                                    caps, fa, fb, gaps)

    ma0, mb0 = witnesses[0]
    status = "unique" if len(witnesses) == 1 else "ambiguous"
    gf = bool(gaps)
    primary = _witness_payload(ma0, mb0, start_year, sa, sb, chron, gf)
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
            _witness_payload(ma, mb, start_year, sa, sb, chron, gf)
            for ma, mb in witnesses
        ],
        "evidence": None,
    }
