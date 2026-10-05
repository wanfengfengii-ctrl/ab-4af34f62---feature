"""chronology_gaps（年表缺测项）的算法与接口测试。

含一个缺测感知的小规模穷举参考实现，与主算法逐目标值、逐见证对照。
"""

from __future__ import annotations

import itertools
import random

import pytest
from fastapi.testclient import TestClient

from app.alignment import align_joint
from app.main import app

client = TestClient(app)


# --------------------------------------------------------------------------- #
# 缺测感知的穷举参考实现
# --------------------------------------------------------------------------- #

def brute_paths_gaps(sample, chron, tol, miss_cap, false_cap,
                     gaps, min_matches=8):
    """穷举单芯所有合法映射：{end: [dict]}。

    规则：缺测下标不可匹配；首匹配前不允许推进年表；首末匹配之间的
    缺测下标免费“跨过”，不计缺环；其余年表项计缺环。
    """
    n, m = len(sample), len(chron)
    gap_set = set(gaps)
    out: dict[int, list[dict]] = {}

    def rec(i, j_lo, pairs, miss_used, false_used):
        if i >= n:
            return
        if 0 < i < n - 1 and false_used < false_cap:
            rec(i + 1, j_lo, pairs, miss_used, false_used + 1)
        prev_j = pairs[-1][1] if pairs else None
        for j in range(j_lo, m):
            if j in gap_set:
                continue
            add_miss = 0
            if prev_j is not None:
                add_miss = sum(
                    1 for q in range(prev_j + 1, j) if q not in gap_set)
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
                all_false = tuple(q for q in range(n) if q not in matched_i)
                all_miss = tuple(
                    q for q in range(j0, j + 1)
                    if q not in matched_j and q not in gap_set)
                all_cross = tuple(
                    g for g in gaps if j0 < g < j)
                mp = {
                    "pairs": tuple((a, b) for a, b, _ in p),
                    "missing": all_miss,
                    "false": all_false,
                    "crossed": all_cross,
                    "mmax": max(dd for _, _, dd in p),
                    "dsum": sum(dd for _, _, dd in p),
                }
                if (len(mp["missing"]) <= miss_cap
                        and len(mp["false"]) <= false_cap
                        and len(p) >= min_matches):
                    out.setdefault(j, []).append(mp)
            else:
                rec(i + 1, j + 1, p, miss_used + add_miss, false_used)

    rec(0, 0, [], 0, 0)
    return out


def brute_joint_gaps(sa, sb, chron, ta, tb, ma_c, fa_c, mb_c, fb_c, gaps):
    pa = brute_paths_gaps(sa, chron, ta, ma_c, fa_c, gaps)
    pb = brute_paths_gaps(sb, chron, tb, mb_c, fb_c, gaps)
    best = None
    opt_pairs = []
    for e in sorted(set(pa) & set(pb)):
        for a, b in itertools.product(pa[e], pb[e]):
            sk = len(a["missing"]) + len(a["false"]) \
                + len(b["missing"]) + len(b["false"])
            score = (sk, max(a["mmax"], b["mmax"]), a["dsum"] + b["dsum"])
            if best is None or score < best:
                best = score
                opt_pairs = [(a, b)]
            elif score == best:
                opt_pairs.append((a, b))
    opt_pairs.sort(key=lambda t: (t[0]["pairs"], t[1]["pairs"]))
    return best, opt_pairs


def _payload(chron, sa, sb, ta, tb, gaps):
    return dict(
        chronology=chron, start_year=1900,
        sample_a=sa, sample_b=sb, tolerance_a=ta, tolerance_b=tb,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=list(gaps),
    )


# --------------------------------------------------------------------------- #
# 构造场景
# --------------------------------------------------------------------------- #

def test_gap_crossed_for_free_between_matches():
    # 真值落在年表 0..8、10..17（缺测下标 9）；两芯取同一 12 环窗口
    # 年表 3..8、10..15，共同终年 15，中间跨过缺测 9。
    chron = list(range(100, 125))
    idx = list(range(0, 9)) + list(range(10, 18))
    window = idx[3:15]
    assert len(window) == 12 and window[-1] == 15
    sa = [chron[i] for i in window]
    sb = sa[:]
    p = _payload(chron, sa, sb, [0] * 12, [0] * 12, gaps=(9,))
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["objectives"] == {
        "total_skips": 0, "max_difference": 0, "difference_sum": 0}
    ma, mb = r["mappings"]["a"], r["mappings"]["b"]
    assert ma["crossed_chronology_gap_indices"] == [9]
    assert ma["crossed_chronology_gap_years"] == [1909]
    assert ma["crossed_gap_count"] == 1
    assert ma["missing_chronology_indices"] == []
    assert mb["crossed_chronology_gap_indices"] == [9]
    assert mb["missing_chronology_indices"] == []
    assert r["common_end_year"] == 1915
    assert ma["end_year"] == mb["end_year"] == 1915
    # 日历年按下标连续换算（缺测也占年）。
    matched_years = [m["year"] for m in ma["matched"]]
    assert matched_years == [1900 + i for i in window]
    # 见证里同样列出。
    w = r["witnesses"][0]["mappings"]["a"]
    assert w["crossed_chronology_gap_indices"] == [9]


def test_gap_item_cannot_be_matched():
    # 缺测下标上的占位值恰好与样芯相等：仍不得被匹配。
    chron = list(range(100, 125))
    chron[5] = 7
    sa = [7] + [chron[i] for i in range(0, 11) if i != 5]
    assert len(sa) == 11
    sb = sa[:]
    p = _payload(chron, sa, sb, [0] * 11, [0] * 11, gaps=(5,))
    r = align_joint(p)
    # 不允许把首环对到 5；任何解也不能在 matched 中出现下标 5。
    for side in ("a", "b"):
        for w in r["witnesses"]:
            idxs = [m["chronology_index"]
                    for m in w["mappings"][side]["matched"]]
            assert 5 not in idxs


def test_gap_before_first_match_not_crossed():
    # 样芯真值落在年表 10..21；缺测下标 3 在首匹配之前，不得计入跨过。
    chron = list(range(100, 130))
    sa = [chron[i] for i in range(10, 22)]
    sb = sa[:]
    p = _payload(chron, sa, sb, [0] * 12, [0] * 12, gaps=(3,))
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["mappings"]["a"]["crossed_chronology_gap_indices"] == []
    assert r["mappings"]["a"]["crossed_gap_count"] == 0


def test_end_on_gap_is_no_solution():
    # 年表值互不相同、容差 0：两芯末环唯一可落之处就是缺测下标 15，
    # 末环不得匹配缺测项 -> no_solution，且 1915 不进候选年。
    chron = list(range(100, 130))
    sa = [chron[i] for i in range(6, 16)]   # 10 环，唯一终年 15
    sb = sa[:]
    p = _payload(chron, sa, sb, [0] * 10, [0] * 10, gaps=(15,))
    r = align_joint(p)
    assert r["status"] == "no_solution", r
    for side in ("a", "b"):
        cev = r["evidence"]["cores"][side]
        assert 1915 not in cev["last_ring_candidate_years"]
    assert r["evidence"]["shared_candidate_years"] == []


def test_gap_does_not_consume_missing_cap():
    # 真值需要跨过 4 个年表位置：3 个木材缺环（达限额 3）+ 1 个缺测。
    chron = list(range(100, 140))
    wood_miss = [4, 6, 8]
    gap = [5]
    used = set(wood_miss) | set(gap)
    idx = [i for i in range(0, 14) if i not in used]
    assert len(idx) == 10
    sa = [chron[i] for i in idx]
    sb = sa[:]
    p = dict(
        chronology=chron, start_year=1900, sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 10, tolerance_b=[0] * 10,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=gap,
    )
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["objectives"]["total_skips"] == 6
    ma = r["mappings"]["a"]
    assert ma["missing_chronology_indices"] == wood_miss
    assert ma["crossed_chronology_gap_indices"] == gap
    # 若缺测误算缺环（共需 4 > 3），则应无解——确保不是这种情况。


def test_wood_missing_immediately_before_final_pair():
    # 回归：木材缺环落在最后一对匹配之前（年表 0..8 后跳过 9，末环对 10）。
    # 旧枚举器在叶子要求 j==end 且 sk_left==0，会漏掉这类 DP 可达映射。
    chron = list(range(100, 130))
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]]
    p = dict(
        chronology=chron, start_year=1900, sample_a=sa, sample_b=sa[:],
        tolerance_a=[0] * 10, tolerance_b=[0] * 10,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
    )
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["objectives"]["total_skips"] == 2
    assert r["mappings"]["a"]["missing_chronology_indices"] == [9]
    sig = [m["chronology_index"] for m in r["mappings"]["a"]["matched"]]
    assert sig[-2:] == [8, 10]


def test_without_gaps_same_case_needs_missing_slot():
    # 对照：省略 chronology_gaps 时，下标 5 是可用项，必须计为木材缺环。
    chron = list(range(100, 140))
    wood_miss = [4, 6, 8]
    used = set(wood_miss) | {5}
    idx = [i for i in range(0, 14) if i not in used]
    sa = [chron[i] for i in idx]
    sb = sa[:]
    # 对齐核心直接放宽缺环限额到 4（API 仍限 0–3），以核对缺环计数。
    p = dict(
        chronology=chron, start_year=1900, sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 10, tolerance_b=[0] * 10,
        max_missing_a=4, max_false_a=4, max_missing_b=4, max_false_b=4,
    )
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["objectives"]["total_skips"] == 8
    ma = r["mappings"]["a"]
    assert ma["missing_chronology_indices"] == [4, 5, 6, 8]
    assert "crossed_chronology_gap_indices" not in ma


def test_no_solution_evidence_separates_gaps_and_wood_skips():
    # A 严格终年 {9}（0 跳环），B 严格终年 {17}（0 跳环）；A 末环值在年表
    # 下标 17 处重复，故放宽后共同终年可为 17：A 需越过 9..16 共 8 个槽，
    # 其中下标 13 是缺测（免费跨过），其余 7 个仍是木材缺环（> 限额 3）。
    chron = list(range(100, 117)) + [109] + list(range(118, 140))
    assert chron[9] == chron[17] == 109 and len(chron) == 40
    sa = [chron[i] for i in range(0, 10)]
    sb = [chron[i] for i in range(8, 18)]
    p = dict(
        chronology=chron, start_year=1900, sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 10, tolerance_b=[0] * 10,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=[13],
    )
    r = align_joint(p)
    assert r["status"] == "no_solution", r
    assert r["evidence"]["reason"] == "common_end_requires_more_skips"
    ev = r["evidence"]
    assert ev["chronology_gaps"]["indices"] == [13]
    first = next(d for d in ev["shared_end_details"]
                 if d["chronology_index"] == 17)
    da = first["a"]
    # 年表缺测（资料缺口）与木材跳环分列：跨过缺测 13 不抵缺环数。
    assert da["chronology_gaps_crossed"] == [13]
    assert da["crossed_gap_years"] == [1913]
    assert da["missing"] == 7 and da["within_limits"] is False
    assert da["total_skips"] == 7
    # B 在该终年 0 跳环即可（窗口 8..17 同样免费跨过缺测 13）。
    assert first["b"]["missing"] == 0
    assert first["b"]["chronology_gaps_crossed"] == [13]
    assert first["b"]["within_limits"] is True
    # 两芯证据各自的限额字段保持原样。
    assert ev["cores"]["a"]["limits"] == {"missing": 3, "false": 3}


# --------------------------------------------------------------------------- #
# 随机穷举对照
# --------------------------------------------------------------------------- #

def _random_gap_case(rng, m=26, na=10, nb=10):
    chron = [rng.randint(0, 9) for _ in range(m)]
    gaps = tuple(sorted(rng.sample(range(2, m - 2), rng.randint(1, 3))))
    usable = [j for j in range(m) if j not in gaps]

    def make_core():
        start_pos = rng.randint(0, len(usable) - na - 2)
        idx = list(usable[start_pos:start_pos + na])
        # 随机把 0-2 个可用位置变木材缺环（顺延取后续可用年表项）
        for _ in range(rng.randint(0, 2)):
            k = rng.randint(1, len(idx) - 2)
            later = [j for j in usable if j > idx[k]
                     and j not in idx[k:]]
            if later:
                idx = idx[:k] + later[:1] + idx[k + 1:]
        vals = [chron[j] + rng.choice([-1, 0, 1]) for j in idx]
        for _ in range(rng.randint(0, 2)):
            k = rng.randint(1, len(vals) - 1)
            vals.insert(k, rng.randint(0, 9))
        tol = [rng.randint(0, 2) for _ in vals]
        return vals, tol

    sa, ta = make_core()
    sb, tb = make_core()
    return chron, sa, sb, ta, tb, gaps


@pytest.mark.parametrize("seed", range(80))
def test_random_gap_cases_match_bruteforce(seed):
    rng = random.Random(1000 + seed)
    chron, sa, sb, ta, tb, gaps = _random_gap_case(rng)
    if len(sa) < 10 or len(sb) < 10:
        return
    p = _payload(chron, sa, sb, ta, tb, gaps)
    best, opt_pairs = brute_joint_gaps(
        sa, sb, chron, ta, tb, 3, 3, 3, 3, gaps)
    r = align_joint(p)
    if best is None:
        assert r["status"] == "no_solution", r
        return
    assert r["status"] in ("unique", "ambiguous")
    assert r["objectives"] == {
        "total_skips": best[0],
        "max_difference": best[1],
        "difference_sum": best[2],
    }
    expect_n = 1 if len(opt_pairs) == 1 else 2
    assert len(r["witnesses"]) == expect_n
    assert r["status"] == ("unique" if expect_n == 1 else "ambiguous")
    first = opt_pairs[0]
    for side, mp_ref in (("a", first[0]), ("b", first[1])):
        got = r["witnesses"][0]["mappings"][side]
        got_sig = tuple((m["sample_ring"] - 1, m["chronology_index"])
                        for m in got["matched"])
        assert got_sig == mp_ref["pairs"]
        assert got["missing_chronology_indices"] == list(mp_ref["missing"])
        assert got["crossed_chronology_gap_indices"] == list(
            mp_ref["crossed"])
        assert got["false_sample_rings"] == [i + 1 for i in mp_ref["false"]]
        assert got["crossed_gap_count"] == len(mp_ref["crossed"])
    if expect_n == 2:
        second = opt_pairs[1]
        for side, mp_ref in (("a", second[0]), ("b", second[1])):
            got = r["witnesses"][1]["mappings"][side]
            got_sig = tuple((m["sample_ring"] - 1, m["chronology_index"])
                            for m in got["matched"])
            assert got_sig == mp_ref["pairs"]


# --------------------------------------------------------------------------- #
# HTTP 422 校验
# --------------------------------------------------------------------------- #

def _ok_body():
    chron = list(range(100, 130))
    sa = [chron[i] for i in range(0, 12)]
    return dict(
        chronology=chron, start_year=1990, sample_a=sa, sample_b=sa[:],
        tolerance_a=[0] * 12, tolerance_b=[0] * 12,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
    )


def _post(body):
    return client.post("/api/dendro/align", json=body)


def test_api_duplicate_gap_index_422():
    body = _ok_body()
    body["chronology_gaps"] = [5, 5]
    r = _post(body)
    assert r.status_code == 422
    errs = r.json()["errors"]
    e = next(e for e in errs if e["loc"][-1] == "chronology_gaps")
    assert "严格递增" in e["msg"]


def test_api_unordered_gap_index_422():
    body = _ok_body()
    body["chronology_gaps"] = [7, 5]
    r = _post(body)
    assert r.status_code == 422
    e = next(e for e in r.json()["errors"]
             if e["loc"][-1] == "chronology_gaps")
    assert "严格递增" in e["msg"]


def test_api_gap_index_out_of_range_422():
    body = _ok_body()
    body["chronology_gaps"] = [30]       # 年表长度 30，合法范围 0..29
    r = _post(body)
    assert r.status_code == 422
    e = next(e for e in r.json()["errors"]
             if e["loc"][-1] == "chronology_gaps")
    assert "越界" in e["msg"]


def test_api_gap_count_bounds():
    body = _ok_body()
    body["chronology_gaps"] = []
    r = _post(body)
    assert r.status_code == 422
    body = _ok_body()
    body["chronology_gaps"] = [1, 2, 3, 4, 5, 6]
    r = _post(body)
    assert r.status_code == 422


def test_api_gap_negative_422():
    body = _ok_body()
    body["chronology_gaps"] = [-1]
    r = _post(body)
    assert r.status_code == 422
    assert any(e["loc"][-1] == "chronology_gaps"
               for e in r.json()["errors"])


def test_gap_immediately_before_final_pair_is_crossed_free():
    # 缺测下标紧贴最后一对匹配之前：叶子路径必须能免费跨过。
    chron = list(range(100, 130))
    sa = [chron[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]]
    p = dict(
        chronology=chron, start_year=1900, sample_a=sa, sample_b=sa[:],
        tolerance_a=[0] * 10, tolerance_b=[0] * 10,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=[9],
    )
    r = align_joint(p)
    assert r["status"] == "unique", r
    assert r["objectives"]["total_skips"] == 0
    ma = r["mappings"]["a"]
    assert ma["crossed_chronology_gap_indices"] == [9]
    assert ma["missing_chronology_indices"] == []


def test_api_omitted_gaps_stays_compatible():
    r = _post(_ok_body())
    assert r.status_code == 200, r.text
    ma = r.json()["mappings"]["a"]
    assert "crossed_chronology_gap_indices" not in ma


def test_api_gaps_enabled_end_to_end():
    chron = list(range(100, 130))
    idx = list(range(0, 9)) + list(range(10, 13))
    sa = [chron[i] for i in idx]
    sb = sa[:]
    body = dict(
        chronology=chron, start_year=1990, sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 12, tolerance_b=[0] * 12,
        max_missing_a=3, max_false_a=3, max_missing_b=3, max_false_b=3,
        chronology_gaps=[9],
    )
    r = _post(body)
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["status"] == "unique"
    assert data["objectives"]["total_skips"] == 0
    assert data["mappings"]["a"]["crossed_chronology_gap_indices"] == [9]
    assert data["mappings"]["a"]["crossed_chronology_gap_years"] == [1999]
