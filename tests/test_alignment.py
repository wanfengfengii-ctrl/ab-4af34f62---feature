"""两芯联合对齐算法测试，含小规模穷举对照。"""

from __future__ import annotations

import itertools
import random

import pytest

from app.alignment import align_joint, Mapping


# --------------------------------------------------------------------------- #
# 穷举参考实现
# --------------------------------------------------------------------------- #

def brute_paths(sample, chron, tol, miss_cap, false_cap, min_matches=8,
                gaps=frozenset()):
    """穷举单芯所有合法映射：{end: [Mapping]}。

    映射即严格双增匹配对 (i_k, j_k)：首对 i_0=0、末对 i=n-1；
    相邻对之间未匹配的样芯环为伪环、真实年表项为缺环；
    窗口内夹着的缺测项记入 gaps_crossed（资料缺口，不算缺环）。
    """
    n, m = len(sample), len(chron)
    out: dict[int, list[Mapping]] = {}

    def rec(i, j_lo, pairs, miss_used, false_used):
        # i：当前待决策的样芯下标；j_lo：年表下一对的最小下标。
        if i >= n:
            return
        # 分支 1：把第 i 环当伪环跳过（首环与末环不允许）。
        if 0 < i < n - 1 and false_used < false_cap:
            rec(i + 1, j_lo, pairs, miss_used, false_used + 1)
        # 分支 2：第 i 环匹配年表第 j 项（缺测项不允许匹配）。
        prev_j = pairs[-1][1] if pairs else None
        for j in range(j_lo, m):
            if j in gaps:
                continue
            if prev_j is None:
                add_miss = 0
            else:
                add_miss = sum(1 for q in range(prev_j + 1, j)
                               if q not in gaps)
            if miss_used + add_miss > miss_cap:
                break
            d = abs(sample[i] - chron[j])
            if d > tol[i]:
                continue
            p = pairs + [(i, j, d)]
            if i == n - 1:
                matched_i = {pp[0] for pp in p}
                matched_j = {pp[1] for pp in p}
                j0 = p[0][1]
                all_false = tuple(
                    q for q in range(n) if q not in matched_i)
                window = range(j0, j + 1)
                all_miss = tuple(
                    q for q in window if q not in matched_j and q not in gaps)
                all_gaps = tuple(q for q in window if q in gaps)
                mp = Mapping(tuple(p), all_miss, all_false, all_gaps)
                if (mp.miss <= miss_cap and mp.fls <= false_cap
                        and len(p) >= min_matches):
                    out.setdefault(j, []).append(mp)
            else:
                rec(i + 1, j + 1, p, miss_used + add_miss, false_used)

    rec(0, 0, [], 0, 0)
    return out


def brute_joint(sa, sb, chron, ta, tb, ma_c, fa_c, mb_c, fb_c,
                gaps=frozenset()):
    pa = brute_paths(sa, chron, ta, ma_c, fa_c, gaps=gaps)
    pb = brute_paths(sb, chron, tb, mb_c, fb_c, gaps=gaps)
    best = None
    pairs_opt = []
    for e in sorted(set(pa) & set(pb)):
        for a, b in itertools.product(pa[e], pb[e]):
            score = (a.skips + b.skips, max(a.mmax, b.mmax),
                     a.dsum + b.dsum)
            if best is None or score < best:
                best = score
                pairs_opt = [(a, b)]
            elif score == best:
                pairs_opt.append((a, b))
    pairs_opt.sort(key=lambda t: (t[0].sig, t[1].sig))
    return best, pairs_opt


def make_payload(chron, sa, sb, ta=None, tb=None, tol=2, start_year=2000,
                 caps=(3, 3, 3, 3), gaps=None):
    p = dict(
        chronology=chron, start_year=start_year,
        sample_a=sa, sample_b=sb,
        tolerance_a=[tol] * len(sa) if ta is None else ta,
        tolerance_b=[tol] * len(sb) if tb is None else tb,
        max_missing_a=caps[0], max_false_a=caps[1],
        max_missing_b=caps[2], max_false_b=caps[3],
    )
    if gaps is not None:
        p["chronology_gaps"] = list(gaps)
    return p


# --------------------------------------------------------------------------- #
# 构造场景
# --------------------------------------------------------------------------- #

def test_basic_unique_alignment():
    chron = [100 + (i * 7) % 50 for i in range(40)]
    sa = [chron[i] for i in range(5, 17)]       # 12 环，终年 16
    sb = [chron[i] for i in [8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
    # A 终年调整为 20 与 B 一致
    sa = [chron[i] for i in range(9, 21)]       # 12 环，终年 20
    # B 在 10 处缺环
    assert sb[1] != sb[2]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 12, tb=[0] * 12))
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2020
    assert res["objectives"] == {
        "total_skips": 1, "max_difference": 0, "difference_sum": 0}
    wa = res["mappings"]["a"]
    wb = res["mappings"]["b"]
    assert wa["missing_chronology_indices"] == []
    assert wa["false_sample_rings"] == []
    assert wb["missing_chronology_indices"] == [10]
    assert wb["false_sample_rings"] == []
    assert wa["match_count"] == 12
    assert wb["match_count"] == 12
    assert len(res["witnesses"]) == 1


def test_false_ring_skipped_in_sample():
    chron = [50 + i * 3 for i in range(40)]
    # A 芯 11 环：第 5 环是伪环（不在年表），其余取年表 3..13，终年 13
    sa = [chron[3], chron[4], chron[5], chron[6], chron[7],
          999,
          chron[8], chron[9], chron[10], chron[11], chron[12 + 1]]
    # 重新设计：年表 3..12 共 10 年，插入 1 个伪环 -> 11 环，终年 12
    sa = [chron[3], chron[4], chron[5], chron[6], chron[7], 999,
          chron[8], chron[9], chron[10], chron[11], chron[12]]
    sb = [chron[i] for i in range(2, 13)]  # 11 环，终年 12
    res = align_joint(make_payload(chron, sa, sb))
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2012
    assert res["objectives"]["total_skips"] == 1
    assert res["mappings"]["a"]["false_sample_rings"] == [6]
    assert res["mappings"]["a"]["missing_chronology_indices"] == []
    assert res["mappings"]["b"]["false_sample_rings"] == []


def test_missing_and_false_caps_zero():
    chron = [100 + (i * 11) % 47 for i in range(35)]
    # 两芯终年同为 18；B 缺年表下标 9。
    sa = [chron[i] for i in range(8, 19)]       # 11 环
    sb = [chron[i] for i in [8, 10, 11, 12, 13, 14, 15, 16, 17, 18]]
    p = make_payload(chron, sa, sb, ta=[0] * 11, tb=[0] * 10,
                     caps=(3, 3, 0, 3))
    res = align_joint(p)
    assert res["status"] == "no_solution"
    assert res["evidence"]["cores"]["b"]["limits"]["missing"] == 0
    # 放开缺环限额后可解
    p["max_missing_b"] = 1
    res = align_joint(p)
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2018
    assert res["mappings"]["b"]["missing_chronology_indices"] == [9]


def test_no_common_end_year():
    # A 精确落在下标 3..14（终年 14，最多再插 3 个缺环 -> 最晚终年 17）；
    # B 精确落在下标 16..25（终年 25）。两者终年集合不相交。
    chron = list(range(100, 130))
    sa = chron[3:15]
    sb = chron[16:26]
    p = make_payload(chron, sa, sb, ta=[0] * 12, tb=[0] * 10)
    res = align_joint(p)
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "no_common_end_year"
    assert res["evidence"]["shared_candidate_years"] == []
    assert res["evidence"]["cores"]["a"][
        "reachable_end_years_within_limits"]
    assert res["evidence"]["cores"]["b"][
        "reachable_end_years_within_limits"] == [2025]


def test_last_ring_unmatchable_evidence():
    chron = list(range(100, 130))
    sa = chron[3:15]
    sb = [101, 102, 103, 104, 105, 106, 107, 108, 109, 999]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 12, tb=[0] * 10))
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "last_ring_unmatchable"
    assert res["evidence"]["cores"]["b"]["last_ring_candidate_years"] == []


def test_skip_limits_exceeded_evidence():
    # 两芯真值都需要 4 个缺环，超过各自限额 3。
    chron = list(range(100, 160))
    sa = [chron[i] for i in list(range(0, 6)) + list(range(10, 14))]
    sb = sa[:]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 10, tb=[0] * 10))
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "skip_limits_exceeded"
    ev = res["evidence"]
    assert ev["cores"]["a"]["min_skips_within_limits"] is None
    # 放宽限额（缺环允许 4 无法通过 API，此处直接核对放宽后信息存在）
    assert ev["cores"]["a"]["min_skips_with_eight_matches"] is not None
    assert "shared_end_details" in ev


def test_insufficient_matchable_rings_evidence():
    # A 芯有 3 个中环值在年表中根本不存在（容差 0），即便全部按伪环跳掉，
    # 伪环限额 3 时仍只剩 7 环可匹配，不足 8 环。
    chron = list(range(100, 130))
    sa = [chron[0], -1, -2, -3, chron[4], chron[5], chron[6],
          chron[7], chron[8], chron[9]]
    sb = chron[0:10]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 10, tb=[0] * 10))
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "insufficient_matchable_rings"
    assert res["evidence"]["cores"]["a"]["min_skips_with_eight_matches"] is None
    assert res["evidence"]["cores"]["a"]["min_skips_with_any_match"] is not None


def test_ambiguity_two_witnesses_sorted_by_mapping_order():
    # 构造两支都存在两种同优落点的情形：年表周期性重复。
    chron = [10] * 25 + [20] * 10
    sa = [10] * 12
    sb = [10] * 10
    p = make_payload(chron, sa, sb, ta=[0] * 12, tb=[0] * 10)
    res = align_joint(p)
    assert res["status"] == "ambiguous"
    assert len(res["witnesses"]) == 2
    # 见证按映射序：第一份终年早于第二份（j0=0 起，最早终年 11）。
    ends = [w["common_end_year"] for w in res["witnesses"]]
    assert ends[0] <= ends[1]
    assert res["common_end_year"] == ends[0]
    sigs = [
        tuple((m["sample_ring"], m["chronology_index"])
              for m in w["mappings"]["a"]["matched"])
        for w in res["witnesses"]
    ]
    assert sigs[0] <= sigs[1]


def test_tiebreak_max_diff_before_sum():
    # 终年同为 26、A 芯同为 4 跳环（缺 2 + 伪 2，8 环匹配）的候选：
    #   胜者 mmax=2、dsum=7（缺年表 22,23，伪环为样芯第 2、4 环）；
    #   败者 mmax=3、dsum=6（差值和更小，但最大差值更大）。
    # 次级目标先压最大差值，故胜者入选，尽管其差值和更大。
    chron = [1, 6, 6, 4, 8, 7, 6, 9, 6, 2, 0, 7, 9, 7, 2, 2, 4,
             7, 1, 2, 9, 7, 0, 2, 4, 7, 4, 5]
    sa = [8, 9, 1, 9, 0, 9, 5, 4, 8, 5]
    ta = [1, 2, 1, 1, 3, 1, 3, 1, 1, 3]
    sb = chron[17:27]                     # 精确落在年表下标 17..26
    tb = [0] * 10
    res = align_joint(make_payload(chron, sa, sb, ta=ta, tb=tb, start_year=2000))
    assert res["common_end_year"] == 2026
    assert res["objectives"]["total_skips"] == 4
    assert res["objectives"]["max_difference"] == 2
    assert res["objectives"]["difference_sum"] == 7
    ma = res["mappings"]["a"]
    assert ma["missing_chronology_indices"] == [22, 23]
    assert ma["false_sample_rings"] == [2, 4]
    assert ma["match_count"] == 8
    assert res["mappings"]["b"]["max_difference"] == 0
    assert res["mappings"]["b"]["match_count"] == 10


# --------------------------------------------------------------------------- #
# 末段紧邻跳环（回归：枚举器必须支持伪环 / 缺环恰好位于末对之前）
# --------------------------------------------------------------------------- #

def test_false_ring_immediately_before_last_pair():
    """倒数第二个样芯环是伪环（尾段伪环），映射仍须被精确枚举。"""
    chron = list(range(100, 140))
    # A：前 8 环对年表 3..10，第 9 环伪环 999，末环对年表 11（共 10 环）。
    sa = [chron[i] for i in range(3, 11)] + [999, chron[11]]
    assert len(sa) == 10
    sb = chron[2:12]                      # 10 环精确落在 2..11，终年 11
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 10, tb=[0] * 10,
                                   caps=(3, 1, 0, 0)))
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2011
    ma = res["mappings"]["a"]
    assert ma["false_sample_rings"] == [9]    # 1 基的倒数第二环
    assert ma["missing_chronology_indices"] == []
    assert ma["match_count"] == 9


def test_missing_immediately_before_common_end():
    """末环对前恰好缺两个真实年表项（尾段缺环），映射仍须被精确枚举。"""
    chron = list(range(100, 140))
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 11]]
    sb = chron[2:12]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 10, tb=[0] * 10,
                                   caps=(2, 0, 0, 0)))
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2011
    ma = res["mappings"]["a"]
    assert ma["missing_chronology_indices"] == [9, 10]
    assert ma["false_sample_rings"] == []
    assert ma["match_count"] == 10
    assert res["objectives"]["total_skips"] == 2


# --------------------------------------------------------------------------- #
# 穷举对照（随机小规模）
# --------------------------------------------------------------------------- #

def _random_case(rng, m=26, na=10, nb=10):
    chron = [rng.randint(0, 9) for _ in range(m)]
    # 以年表子序列为真值，注入噪声与缺/伪环
    def make_core():
        n = na
        start = rng.randint(0, m - n - 3)
        idx = list(range(start, start + n))
        # 随机把 0-2 个位置变缺环（抽掉年表项，顺延取后续年）
        miss = rng.randint(0, 2)
        for _ in range(miss):
            if idx[-1] + 1 < m:
                k = rng.randint(1, len(idx) - 2)
                idx = idx[:k] + [x + 1 for x in idx[k:]]
        vals = [chron[j] + rng.choice([-1, 0, 1]) for j in idx]
        # 随机插伪环（同时补逐环容差）
        false = rng.randint(0, 2)
        for _ in range(false):
            k = rng.randint(1, len(vals) - 1)
            vals.insert(k, rng.randint(0, 9))
        tol = [rng.randint(0, 2) for _ in vals]
        return vals, tol
    sa, ta = make_core()
    sb, tb = make_core()
    return chron, sa, sb, ta, tb


@pytest.mark.parametrize("seed", range(60))
def test_random_cases_match_bruteforce(seed):
    rng = random.Random(seed)
    chron, sa, sb, ta, tb = _random_case(rng)
    if len(sa) < 10 or len(sb) < 10:
        return
    p = dict(
        chronology=chron, start_year=1900,
        sample_a=sa, sample_b=sb, tolerance_a=ta, tolerance_b=tb,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
    )
    best, opt_pairs = brute_joint(sa, sb, chron, ta, tb, 3, 3, 3, 3)
    res = align_joint(p)
    if best is None:
        assert res["status"] == "no_solution"
        return
    assert res["status"] in ("unique", "ambiguous")
    assert res["objectives"] == {
        "total_skips": best[0],
        "max_difference": best[1],
        "difference_sum": best[2],
    }
    assert res["common_end_year"] - 1900 == opt_pairs[0][0].pairs[-1][1]
    # 见证数量与映射序前两份
    expect_n = 1 if len(opt_pairs) == 1 else 2
    assert len(res["witnesses"]) == expect_n
    assert res["status"] == ("unique" if expect_n == 1 else "ambiguous")
    first = opt_pairs[0]
    got_a = tuple((m["sample_ring"] - 1, m["chronology_index"])
                  for m in res["witnesses"][0]["mappings"]["a"]["matched"])
    assert got_a == first[0].sig
    got_b = tuple((m["sample_ring"] - 1, m["chronology_index"])
                  for m in res["witnesses"][0]["mappings"]["b"]["matched"])
    assert got_b == first[1].sig
    if expect_n == 2:
        second = opt_pairs[1]
        got_a2 = tuple((m["sample_ring"] - 1, m["chronology_index"])
                       for m in res["witnesses"][1]["mappings"]["a"]["matched"])
        got_b2 = tuple((m["sample_ring"] - 1, m["chronology_index"])
                       for m in res["witnesses"][1]["mappings"]["b"]["matched"])
        assert (got_a2, got_b2) == (second[0].sig, second[1].sig)
    # 主映射即第一见证
    assert res["mappings"] == res["witnesses"][0]["mappings"]


def test_witnesses_payload_internal_consistency():
    chron = [100 + (i * 7) % 50 for i in range(32)]
    sa = [chron[i] + ((i * 3) % 2) for i in range(4, 18)]   # 14 环
    sb = [chron[i] for i in [6, 7, 9, 10, 11, 13, 14, 15, 16, 17]]
    res = align_joint(dict(
        chronology=chron, start_year=1850, sample_a=sa, sample_b=sb,
        tolerance_a=[2] * 14, tolerance_b=[1] * 10,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3))
    assert res["status"] in ("unique", "ambiguous")
    for w in res["witnesses"]:
        for key, sample, tol in (("a", sa, [2] * 14), ("b", sb, [1] * 10)):
            mp = w["mappings"][key]
            assert mp["match_count"] >= 8
            assert mp["end_year"] == w["common_end_year"]
            years = [m["year"] for m in mp["matched"]]
            assert years == sorted(set(years))
            rings = [m["sample_ring"] for m in mp["matched"]]
            assert rings == sorted(set(rings))
            for m in mp["matched"]:
                assert m["difference"] == abs(
                    m["sample_value"] - m["chronology_value"])
                assert m["difference"] <= tol[m["sample_ring"] - 1]
            assert mp["missing_count"] == len(
                mp["missing_chronology_indices"])
            assert mp["false_count"] == len(mp["false_sample_rings"])
            assert mp["match_count"] == len(sample) - mp["false_count"]
            assert mp["max_difference"] == max(
                m["difference"] for m in mp["matched"])
            assert mp["difference_sum"] == sum(
                m["difference"] for m in mp["matched"])
        # 两芯终年一致
        assert w["mappings"]["a"]["end_year"] == w["mappings"]["b"]["end_year"]


# --------------------------------------------------------------------------- #
# chronology_gaps（年表缺测项）
# --------------------------------------------------------------------------- #

def test_gap_crossed_free_and_reported():
    """缺测项可被两芯共同跨越，不耗缺环限额、不计跳环，并分别列示。"""
    chron = list(range(100, 140))
    # A 连续 12 环覆盖下标 0..11（跨缺测 5）；B 另有真实缺环 8。
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12]]
    sb = [chron[i] for i in [0, 1, 2, 3, 4, 6, 7, 9, 10, 11, 12]]
    res = align_joint(make_payload(chron, sa, sb, ta=[0] * 12, tb=[0] * 11,
                                   gaps=[5]))
    assert res["status"] == "unique"
    assert res["common_end_year"] == 2012
    # 总跳环只计木材缺环（B 的 [8]），缺测跨越不计入。
    assert res["objectives"] == {
        "total_skips": 1, "max_difference": 0, "difference_sum": 0}
    ma, mb = res["mappings"]["a"], res["mappings"]["b"]
    assert ma["missing_chronology_indices"] == []
    assert ma["crossed_chronology_gap_indices"] == [5]
    assert ma["crossed_chronology_gap_years"] == [2005]
    assert ma["crossed_gap_count"] == 1
    assert mb["missing_chronology_indices"] == [8]
    assert mb["crossed_chronology_gap_indices"] == [5]
    assert mb["crossed_gap_count"] == 1
    # 每份同优见证都分别列示实际跨过的缺测下标 / 年份。
    for w in res["witnesses"]:
        assert w["mappings"]["a"]["crossed_chronology_gap_indices"] == [5]
        assert w["mappings"]["b"]["crossed_chronology_gap_years"] == [2005]


def test_gap_never_matched_and_year_continuous():
    """缺测项不得匹配；日历年份按下标连续换算，不因缺测断开。"""
    chron = list(range(100, 140))
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10,
                                   gaps=[5]))
    assert res["status"] == "unique"
    mp = res["mappings"]["a"]
    idxs = [m["chronology_index"] for m in mp["matched"]]
    assert 5 not in idxs
    assert idxs == [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]
    years = [m["year"] for m in mp["matched"]]
    assert years == [2000 + j for j in idxs]
    assert mp["crossed_chronology_gap_indices"] == [5]
    assert mp["crossed_chronology_gap_years"] == [2005]


def test_gap_outside_window_not_counted():
    """首末匹配之外的缺测项既不匹配也不计入跨越。"""
    chron = list(range(100, 130))
    sa = chron[10:20]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10,
                                   gaps=[5, 15, 25]))
    assert res["status"] == "unique"
    mp = res["mappings"]["a"]
    assert mp["crossed_chronology_gap_indices"] == [15]
    assert mp["missing_chronology_indices"] == []


def test_gap_immediately_before_common_end():
    """共同终年前一个位置是缺测项：可跨越，且不再需要任何真实跳环。"""
    chron = list(range(100, 140))
    # 前 10 环落下标 0..9，末环落 11，中间只夹缺测 10（无真实缺环）。
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11]]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 11, tb=[0] * 11,
                                   gaps=[10], caps=(0, 0, 0, 0)))
    assert res["status"] == "unique"
    assert res["objectives"]["total_skips"] == 0
    mp = res["mappings"]["a"]
    assert mp["missing_chronology_indices"] == []
    assert mp["crossed_chronology_gap_indices"] == [10]
    assert mp["end_year"] == 2011


def test_gap_with_real_missing_over_cap_returns_no_solution():
    """跨越缺测不解决真实缺环：真值还缺真实年，超限额仍 no_solution。"""
    chron = list(range(100, 140))
    # 真值窗口跨缺测 5，且两芯各自还缺真实年表项 3（共需 1 个真实缺环）。
    sa = [chron[i] for i in [0, 1, 2, 4, 6, 7, 8, 9, 10, 11]]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10,
                                   gaps=[5], caps=(0, 3, 0, 3)))
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "skip_limits_exceeded"
    # 证据中年表缺测与木材跳环分开说明。
    detail = res["evidence"]["shared_end_details"][0]["a"]
    assert detail["total_skips"] == 1          # 只含木材缺环
    assert detail["chronology_gap_indices"] == [5]
    assert detail["chronology_gap_count"] == 1
    assert detail["within_limits"] is False
    # 缺测 5 在放宽窗口内可跨越（失败原因是真实缺环限额，不是缺测本身）。
    gap_ev = res["evidence"]["cores"]["a"]["chronology_gaps"]
    assert gap_ev == [{
        "chronology_index": 5, "year": 2005,
        "crossable_with_eight_matches": True}]

    # 放宽缺环限额到 1 后可解：缺测跨越仍不计入跳环。
    res2 = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10,
                                    gaps=[5], caps=(1, 3, 1, 3)))
    assert res2["status"] == "unique"
    assert res2["objectives"]["total_skips"] == 2
    mp = res2["mappings"]["a"]
    assert mp["missing_chronology_indices"] == [3]
    assert mp["crossed_chronology_gap_indices"] == [5]


def test_gap_at_common_end_blocks_and_reports_gaps_separately():
    """共同终年只能落在缺测项：no_solution + 缺口与跳环分列的两芯证据。"""
    chron = list(range(100, 130))
    sa = chron[0:10]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10,
                                   gaps=[9]))
    assert res["status"] == "no_solution"
    assert res["evidence"]["reason"] == "common_end_blocked_by_gaps"
    assert res["evidence"]["chronology_gaps"] == {
        "indices": [9], "years": [2009]}
    # 末环候选不含缺测年。
    assert 2009 not in res["evidence"]["cores"]["a"][
        "last_ring_candidate_years"]
    for side in ("a", "b"):
        assert res["evidence"]["cores"][side]["chronology_gaps"] == [{
            "chronology_index": 9, "year": 2009,
            "crossable_with_eight_matches": False}]


def test_gap_not_in_response_when_omitted():
    """省略 chronology_gaps 时响应形状保持旧版（无任何 gap 字段）。"""
    chron = list(range(100, 130))
    sa = chron[0:10]
    res = align_joint(make_payload(chron, sa, sa[:], ta=[0] * 10, tb=[0] * 10))
    mp = res["mappings"]["a"]
    assert "crossed_chronology_gap_indices" not in mp
    assert "crossed_chronology_gap_years" not in mp
    assert "crossed_gap_count" not in mp
    assert res["evidence"] is None


def test_gap_omitted_same_inputs_identical_objectives():
    """同一输入无缺测时，加与不加不存在的缺口结果等价。"""
    chron = [100 + (i * 7) % 50 for i in range(30)]
    sa = [chron[i] for i in range(9, 20)]
    sb = [chron[i] for i in [9, 10, 12, 13, 14, 15, 16, 17, 18, 19]]
    p = make_payload(chron, sa, sb, ta=[0] * 11, tb=[0] * 10)
    r1 = align_joint(p)
    p2 = dict(p); p2["chronology_gaps"] = [28]   # 落在窗口外
    r2 = align_joint(p2)
    assert r1["objectives"] == r2["objectives"]
    assert r1["common_end_year"] == r2["common_end_year"]
    assert r2["mappings"]["a"]["crossed_chronology_gap_indices"] == []
    assert r2["mappings"]["b"]["crossed_chronology_gap_indices"] == []


@pytest.mark.parametrize("seed", range(40))
def test_random_gap_cases_match_bruteforce(seed):
    """含 1–3 个缺测项的随机算例：与穷举参考逐目标值、逐见证对照。"""
    rng = random.Random(seed)
    m = 26
    chron = [rng.randint(0, 9) for _ in range(m)]
    gaps = frozenset(rng.sample(range(2, m - 2), rng.randint(1, 3)))

    def make_core():
        n = 10
        start = rng.randint(0, m - n - 4)
        idx = list(range(start, start + n))
        miss = rng.randint(0, 2)
        for _ in range(miss):
            if idx[-1] + 1 < m:
                k = rng.randint(1, len(idx) - 2)
                idx = idx[:k] + [x + 1 for x in idx[k:]]
        vals = [chron[j] + rng.choice([-1, 0, 1]) for j in idx]
        false = rng.randint(0, 2)
        for _ in range(false):
            k = rng.randint(1, len(vals) - 1)
            vals.insert(k, rng.randint(0, 9))
        tol = [rng.randint(0, 2) for _ in vals]
        return vals, tol

    sa, ta = make_core()
    sb, tb = make_core()
    if len(sa) < 10 or len(sb) < 10:
        return
    p = dict(
        chronology=chron, start_year=1900,
        sample_a=sa, sample_b=sb, tolerance_a=ta, tolerance_b=tb,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=sorted(gaps),
    )
    best, opt_pairs = brute_joint(sa, sb, chron, ta, tb, 3, 3, 3, 3, gaps)
    res = align_joint(p)
    if best is None:
        assert res["status"] == "no_solution"
        return
    assert res["status"] in ("unique", "ambiguous")
    assert res["objectives"] == {
        "total_skips": best[0],
        "max_difference": best[1],
        "difference_sum": best[2],
    }
    assert len(res["witnesses"]) == (1 if len(opt_pairs) == 1 else 2)
    first = opt_pairs[0]
    ga = tuple((mm["sample_ring"] - 1, mm["chronology_index"])
               for mm in res["witnesses"][0]["mappings"]["a"]["matched"])
    gb = tuple((mm["sample_ring"] - 1, mm["chronology_index"])
               for mm in res["witnesses"][0]["mappings"]["b"]["matched"])
    assert (ga, gb) == (first[0].sig, first[1].sig)
    # 跨越缺测与穷举一致
    ca = res["witnesses"][0]["mappings"]["a"]["crossed_chronology_gap_indices"]
    cb = res["witnesses"][0]["mappings"]["b"]["crossed_chronology_gap_indices"]
    assert ca == list(first[0].gaps_crossed)
    assert cb == list(first[1].gaps_crossed)
    # 缺环只含真实年表项、见证间终年一致
    for w in res["witnesses"]:
        for key in ("a", "b"):
            mp = w["mappings"][key]
            assert set(mp["missing_chronology_indices"]).isdisjoint(gaps)
            assert set(mp["crossed_chronology_gap_indices"]) <= gaps
        assert w["mappings"]["a"]["end_year"] == w["mappings"]["b"]["end_year"]
