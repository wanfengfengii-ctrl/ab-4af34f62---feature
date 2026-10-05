"""/api/dendro/align 与 /health 的接口测试。"""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


BASE = dict(
    chronology=[100 + (i * 7) % 50 for i in range(30)],
    start_year=1990,
    sample_a=None,
    sample_b=None,
    tolerance_a=None,
    tolerance_b=None,
    max_missing_a=3, max_false_a=3,
    max_missing_b=3, max_false_b=3,
)


def _req():
    chron = BASE["chronology"]
    # A 芯精确落在下标 9..19（终年 19）；B 芯 10 环落在同区间、缺下标 11。
    sa = [chron[i] for i in range(9, 20)]
    sb = [chron[i] for i in [9, 10, 12, 13, 14, 15, 16, 17, 18, 19]]
    p = dict(BASE)
    p.update(
        sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 11, tolerance_b=[0] * 10,
    )
    return p


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_align_ok():
    r = client.post("/api/dendro/align", json=_req())
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "unique"
    assert body["common_end_year"] == 2009
    assert body["objectives"]["total_skips"] == 1
    assert body["mappings"]["b"]["missing_chronology_indices"] == [11]
    assert len(body["witnesses"]) == 1


def test_chronology_too_short():
    p = _req()
    p["chronology"] = [1] * 24
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422
    body = r.json()
    assert body["status"] == "validation_error"
    assert any("chronology" in str(e.get("loc", "")) for e in body["errors"])


def test_chronology_too_long():
    p = _req()
    p["chronology"] = [1] * 61
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_sample_length_bounds():
    p = _req()
    p["sample_a"] = [1] * 9
    p["tolerance_a"] = [1] * 9
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422
    assert body_status(r) == "validation_error"


def test_tolerance_length_mismatch():
    p = _req()
    p["tolerance_a"] = [0] * 10      # sample_a 有 11 环
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422
    assert "tolerance_a" in r.text


def test_negative_tolerance_rejected():
    p = _req()
    p["tolerance_b"] = [-1] * 10
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_limit_out_of_range():
    p = _req()
    p["max_missing_a"] = 4
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422
    p = _req()
    p["max_false_b"] = -1
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_extra_field_rejected():
    p = _req()
    p["unexpected"] = 1
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_missing_field_rejected():
    p = _req()
    del p["sample_b"]
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_non_integer_rejected():
    p = _req()
    p["sample_a"] = [1.5] * 12
    p["tolerance_a"] = [0] * 12
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 422


def test_no_solution_shape():
    p = _req()
    p["sample_b"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 99999]
    p["tolerance_b"] = [0] * 10
    r = client.post("/api/dendro/align", json=p)
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "no_solution"
    assert body["common_end_year"] is None
    assert body["mappings"] is None
    ev = body["evidence"]
    assert ev["reason"] == "last_ring_unmatchable"
    assert "a" in ev["cores"] and "b" in ev["cores"]
    assert ev["cores"]["b"]["last_ring_candidate_years"] == []
    assert "shared_end_details" in ev


# --------------------------------------------------------------------------- #
# chronology_gaps
# --------------------------------------------------------------------------- #

def _gap_req(gaps=None, caps=(3, 3, 3, 3)):
    chron = list(range(100, 130))
    # 两芯真值共同终年 20、窗口跨缺测 12；B 另缺真实年表项 14。
    # 样芯取值本身不含缺测年（缺测年只被跨越，不产生伪环）。
    sa = [chron[i]
          for i in [9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20]]
    sb = [chron[i]
          for i in [9, 10, 11, 13, 15, 16, 17, 18, 19, 20]]
    p = dict(
        chronology=chron, start_year=1990,
        sample_a=sa, sample_b=sb,
        tolerance_a=[0] * 11, tolerance_b=[0] * 10,
        max_missing_a=caps[0], max_false_a=caps[1],
        max_missing_b=caps[2], max_false_b=caps[3],
    )
    if gaps is not None:
        p["chronology_gaps"] = gaps
    return p


def test_gap_alignment_ok():
    # 缺测下标 12 落在两芯窗口内（9..20），可共同跨越；B 另有真实缺环 14。
    r = client.post("/api/dendro/align", json=_gap_req([12]))
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["status"] == "unique"
    assert body["common_end_year"] == 2010
    # 总跳环只计 B 的木材缺环，缺测跨越不计入。
    assert body["objectives"]["total_skips"] == 1
    ma, mb = body["mappings"]["a"], body["mappings"]["b"]
    assert ma["missing_chronology_indices"] == []
    assert ma["crossed_chronology_gap_indices"] == [12]
    assert ma["crossed_chronology_gap_years"] == [2002]
    assert mb["missing_chronology_indices"] == [14]
    assert mb["crossed_chronology_gap_indices"] == [12]
    # 缺测占位值绝不与任何环匹配。
    for mp in (ma, mb):
        assert 12 not in [m["chronology_index"] for m in mp["matched"]]
    # 每份见证都分别列示实际跨过的缺测下标 / 年份。
    for w in body["witnesses"]:
        for key in ("a", "b"):
            wm = w["mappings"][key]
            assert wm["crossed_chronology_gap_indices"] == [12]
            assert wm["crossed_chronology_gap_years"] == [2002]


def test_gap_omitted_response_has_no_gap_fields():
    r = client.post("/api/dendro/align", json=_gap_req())
    assert r.status_code == 200
    mp = r.json()["mappings"]["a"]
    assert "crossed_chronology_gap_indices" not in mp
    assert "crossed_chronology_gap_years" not in mp


def test_gap_year_conversion_continuous():
    r = client.post("/api/dendro/align", json=_gap_req([12]))
    matched = r.json()["mappings"]["a"]["matched"]
    for m in matched:
        assert m["year"] == 1990 + m["chronology_index"]


def test_gap_duplicate_index_422():
    r = client.post("/api/dendro/align", json=_gap_req([12, 12]))
    assert r.status_code == 422
    body = r.json()
    assert body["status"] == "validation_error"
    locs = [tuple(e["loc"]) for e in body["errors"]]
    assert ("body", "chronology_gaps") in locs


def test_gap_unsorted_index_422():
    r = client.post("/api/dendro/align", json=_gap_req([15, 12]))
    assert r.status_code == 422
    assert ("body", "chronology_gaps") in {
        tuple(e["loc"]) for e in r.json()["errors"]}


def test_gap_out_of_range_422():
    r = client.post("/api/dendro/align", json=_gap_req([30]))
    assert r.status_code == 422
    assert "chronology_gaps" in r.text


def test_gap_negative_index_422():
    r = client.post("/api/dendro/align", json=_gap_req([-1]))
    assert r.status_code == 422


def test_gap_empty_and_too_many_422():
    r = client.post("/api/dendro/align", json=_gap_req([]))
    assert r.status_code == 422
    r = client.post("/api/dendro/align", json=_gap_req([1, 2, 3, 4, 5, 6]))
    assert r.status_code == 422


def test_gap_at_common_end_no_solution_evidence():
    # 共同终年只能落在缺测项：no_solution，证据区分缺测与跳环。
    r = client.post("/api/dendro/align", json=_gap_req([20]))
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "no_solution"
    ev = body["evidence"]
    assert ev["reason"] == "common_end_blocked_by_gaps"
    assert ev["chronology_gaps"]["indices"] == [20]
    for side in ("a", "b"):
        gev = ev["cores"][side]["chronology_gaps"]
        assert gev[0]["chronology_index"] == 20
        assert gev[0]["crossable_with_eight_matches"] is False


def body_status(r):
    return r.json()["status"]
