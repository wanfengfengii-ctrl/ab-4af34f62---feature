"""两芯联合树轮对齐 Web 服务。"""

from __future__ import annotations

from typing import List, Optional

from fastapi import FastAPI
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, ConfigDict, field_validator, model_validator

from .alignment import align_joint


# 模型级校验中可定位到具体字段的错误前缀；异常处理器据此修正 loc。
_FIELD_PREFIX = "chronology_gaps:"


def _sanitize_errors(errors: list[dict]) -> list[dict]:
    """Pydantic v2 的 ctx 可能携带不可 JSON 序列化的异常对象，统一清洗。"""
    clean = []
    for e in errors:
        item = {k: v for k, v in e.items() if k != "ctx"}
        msg = str(item.get("msg", ""))
        if _FIELD_PREFIX in msg:
            # 跨字段校验失败时把错误定位到 chronology_gaps；
            # Pydantic 会给 ValueError 加 "Value error, " 前缀，故用子串定位。
            item["loc"] = ["body", "chronology_gaps"]
            item["msg"] = msg.split(_FIELD_PREFIX, 1)[1].strip()
        ctx = e.get("ctx")
        if ctx:
            item["ctx"] = {
                k: str(v) if not isinstance(
                    v, (str, int, float, bool, type(None), list, dict)
                ) else v
                for k, v in ctx.items()
            }
        clean.append(item)
    return clean


class AlignRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    chronology: List[int] = Field(
        ..., description="地区年表生长指数（整数），长度 25–60。")
    start_year: int = Field(
        ..., description="年表首项对应的日历年。")
    sample_a: List[int] = Field(..., min_length=10, max_length=24,
                                description="A 芯生长指数序列，长度 10–24。")
    sample_b: List[int] = Field(..., min_length=10, max_length=24,
                                description="B 芯生长指数序列，长度 10–24。")
    tolerance_a: List[int] = Field(
        ..., description="A 芯逐环容差（非负整数），与 sample_a 等长。")
    tolerance_b: List[int] = Field(
        ..., description="B 芯逐环容差（非负整数），与 sample_b 等长。")
    max_missing_a: int = Field(..., ge=0, le=3, description="A 芯缺环限额 0–3。")
    max_false_a: int = Field(..., ge=0, le=3, description="A 芯伪环限额 0–3。")
    max_missing_b: int = Field(..., ge=0, le=3, description="B 芯缺环限额 0–3。")
    max_false_b: int = Field(..., ge=0, le=3, description="B 芯伪环限额 0–3。")
    chronology_gaps: Optional[List[int]] = Field(
        None, min_length=1, max_length=5,
        description="年表缺测项零基下标：1–5 个、严格递增、位于年表范围内；"
                    "缺测项不可匹配，可在首末匹配之间免费跨过。")

    @field_validator("chronology")
    @classmethod
    def _check_chronology(cls, v: List[int]) -> List[int]:
        if not 25 <= len(v) <= 60:
            raise ValueError("chronology 长度必须在 25 至 60 之间")
        return v

    @field_validator("sample_a", "sample_b")
    @classmethod
    def _check_sample(cls, v: List[int]) -> List[int]:
        if not 10 <= len(v) <= 24:
            raise ValueError("样芯序列长度必须在 10 至 24 之间")
        return v

    @field_validator("tolerance_a", "tolerance_b")
    @classmethod
    def _check_tol_nonneg(cls, v: List[int]) -> List[int]:
        if any(x < 0 for x in v):
            raise ValueError("逐环容差必须为非负整数")
        return v

    @field_validator("chronology_gaps")
    @classmethod
    def _check_gap_values(cls, v: Optional[List[int]]) -> Optional[List[int]]:
        if v is not None and any(x < 0 for x in v):
            raise ValueError("chronology_gaps 下标必须为非负整数")
        return v

    @model_validator(mode="after")
    def _check_cross_fields(self) -> "AlignRequest":
        if len(self.tolerance_a) != len(self.sample_a):
            raise ValueError("tolerance_a 长度必须与 sample_a 相同")
        if len(self.tolerance_b) != len(self.sample_b):
            raise ValueError("tolerance_b 长度必须与 sample_b 相同")
        gaps = self.chronology_gaps
        if gaps is not None:
            m = len(self.chronology)
            if any(b <= a for a, b in zip(gaps, gaps[1:])):
                # 同时覆盖重复（b == a）与乱序（b < a）。
                raise ValueError(
                    f"{_FIELD_PREFIX} 下标必须严格递增，不得重复或乱序")
            if gaps[-1] >= m:
                raise ValueError(
                    f"{_FIELD_PREFIX} 下标越界：必须位于 [0, {m - 1}] 内")
        return self


app = FastAPI(
    title="木构文物两芯联合树轮对齐服务",
    version="1.0.0",
    description=(
        "将同一根梁材的两支树芯共同对入地区年表：先最小化两芯缺环与伪环"
        "总数，再依次最小化最大差值与差值总和，并对同优映射判定唯一 / 歧义。"
    ),
)


@app.exception_handler(RequestValidationError)
async def _validation_handler(request, exc: RequestValidationError):
    return JSONResponse(
        status_code=422,
        content={
            "status": "validation_error",
            "message": "请求参数未通过校验",
            "errors": _sanitize_errors(exc.errors()),
        },
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/api/dendro/align")
def align(req: AlignRequest) -> dict:
    return align_joint(req.model_dump())
