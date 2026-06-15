"""
FastAPI 应用入口 —— 生物学年龄预测 API + Agent Chat。
"""

import io
import time
import uuid
from pathlib import Path
from typing import Dict, Optional

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, StreamingResponse
from jinja2 import Environment, FileSystemLoader

from app.model import load_model, predict_bio_age
from app.utils import generate_demo, read_input_file
from agent.memory import ConversationMemory
from agent.orchestrator import AgentOrchestrator, _classify_prediction_error

# --------------------------------------------------------------------------
# FastAPI 应用
# --------------------------------------------------------------------------

app = FastAPI(
    title="Biological Age Predictor",
    description="基于甲基化数据的生物学年龄预测 API",
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

TEMPLATES_DIR = Path(__file__).resolve().parent / "templates"
jinja_env = Environment(
    loader=FileSystemLoader(str(TEMPLATES_DIR)),
    auto_reload=False,
)


def render_template(name: str, context: dict = None) -> str:
    """渲染 Jinja2 模板并返回 HTML 字符串。"""
    template = jinja_env.get_template(name)
    return template.render(context or {})

# --------------------------------------------------------------------------
# 全局状态
# --------------------------------------------------------------------------

orchestrator = AgentOrchestrator()
sessions: Dict[str, ConversationMemory] = {}


def _get_session(session_id: Optional[str] = None) -> ConversationMemory:
    """获取或创建 session 记忆。"""
    if session_id and session_id in sessions:
        return sessions[session_id]
    new_id = uuid.uuid4().hex[:8]
    sessions[new_id] = ConversationMemory()
    return sessions[new_id]


# --------------------------------------------------------------------------
# 事件
# --------------------------------------------------------------------------

@app.on_event("startup")
async def startup():
    load_model()


# --------------------------------------------------------------------------
# 路由：首页
# --------------------------------------------------------------------------

@app.get("/", response_class=HTMLResponse)
async def home():
    return render_template("index.html")


# --------------------------------------------------------------------------
# 路由：下载示例
# --------------------------------------------------------------------------

@app.get("/download-demo")
async def download_demo():
    from app.model import train_cpg_names
    demo_file = generate_demo(train_cpg_names)
    return StreamingResponse(
        demo_file,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": "attachment; filename=bio_age_demo_input_file.xlsx"},
    )


# --------------------------------------------------------------------------
# 路由：经典版预测
# --------------------------------------------------------------------------

@app.post("/predict")
async def predict(
    file: UploadFile = File(..., description="甲基化数据文件"),
    fill_missing: bool = Query(default=True),
):
    t0 = time.time()
    try:
        content = await file.read()
        df = read_input_file(content, file.filename)
        result_df = predict_bio_age(df, fill_missing=fill_missing)
        output = io.BytesIO()
        result_df.to_excel(output, index=False, engine="openpyxl")
        output.seek(0)
        METRICS["total_requests"] += 1
        METRICS["total_duration_ms"] += (time.time() - t0) * 1000
        return StreamingResponse(
            output,
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={"Content-Disposition": "attachment; filename=bio_age_prediction_result.xlsx"},
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=_classify_prediction_error(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=_classify_prediction_error(e))


# --------------------------------------------------------------------------
# 路由：Agent Chat
# --------------------------------------------------------------------------

@app.post("/chat")
async def chat(
    message: str = Form(...),
    file: Optional[UploadFile] = File(None),
    session_id: Optional[str] = Form(None),
):
    t0 = time.time()
    try:
        memory = _get_session(session_id)
        current_sid = next(k for k, v in sessions.items() if v is memory)

        methylation_df = None
        existing_results = None
        if file and file.filename:
            content = await file.read()
            if len(content) > 0:
                df = read_input_file(content, file.filename)
                # 判断是否为预测结果文件（有 Predicted_Age 列）
                if "Predicted_Age" in df.columns:
                    existing_results = df.to_dict(orient="records")
                else:
                    methylation_df = df

        result = await orchestrator.run(
            user_input=message,
            methylation_data=methylation_df,
            memory=memory,
            existing_results=existing_results,
        )

        reply = _format_agent_reply(result, methylation_df is not None)
        memory.add_turn(message, reply, analysis_type=result.get("mode"))

        duration_ms = (time.time() - t0) * 1000
        METRICS["total_requests"] += 1
        METRICS["total_duration_ms"] += duration_ms

        return {
            "reply": reply,
            "session_id": current_sid,
            "trace": result,
            "duration_ms": round(duration_ms, 1),
        }

    except ValueError as e:
        return {"reply": _classify_prediction_error(e), "session_id": session_id or "", "trace": {}}
    except Exception as e:
        return {"reply": f"❌ 系统异常：{str(e)}", "session_id": session_id or "", "trace": {}}


def _format_agent_reply(result: dict, has_file: bool) -> str:
    """将 Agent 结果格式化为可读文本。

    优先级：LLM 报告 > 结构化模板。
    """
    mode = result.get("mode", "knowledge")

    # ★ 追问优先（必须在 guidance 之前检查）
    if result.get("type") == "follow_up":
        llm_answer = result.get("llm_answer")
        if llm_answer:
            return llm_answer
        return result.get("message", "请具体描述你想了解的内容。")

    # ★ 预测错误 → 友好修正指引
    if result.get("type") == "prediction_error":
        return result.get("message", "⚠️ 分析异常，请检查数据格式后重试。")

    # 无数据 → 引导
    if mode in ("guidance", "knowledge"):
        return (
            "🧬 **生物学年龄分析助手**\n\n"
            "请上传你的甲基化数据文件（支持 xlsx/csv/tsv/txt），"
            "我会为你预测生物学年龄并生成分析报告。"
        )

    # LLM 报告优先
    llm_report = result.get("llm_report")
    if llm_report and result.get("report_source") == "llm":
        ktext = result.get("knowledge_text", "")
        if ktext:
            return llm_report + "\n\n" + ktext
        return llm_report

    # 模板兜底
    tool_data = result.get("tool_results", {})
    meth_data = tool_data.get("methylation_clock")
    if not meth_data:
        return "分析完成，但未获取到甲基化预测结果。"

    if isinstance(meth_data, dict):
        rows = [meth_data]
    elif isinstance(meth_data, list):
        rows = meth_data
    else:
        return f"分析完成，但结果格式异常: {type(meth_data)}"

    if len(rows) == 0:
        return "分析完成，但结果为空。"

    import datetime
    today = datetime.date.today().strftime("%Y-%m-%d")
    lines = [f"📊 **分析结果（{today}）**\n"]

    for i, row in enumerate(rows):
        sid = (
            row.get("Sample_ID")
            or row.get("SampleID")
            or row.get("sample_id")
            or f"样本{i + 1}"
        )
        age = row.get("Predicted_Age", "N/A")
        risk = row.get("Risk_Level", "N/A")
        ca = row.get("Chronological_Age")
        gap = row.get("Age_Gap")
        suggestion = row.get("Clinical_Suggestion", "")

        lines.append(f"**{sid}**")
        lines.append(f"- 生物学年龄：**{age} 岁**")
        if ca and ca != "N/A" and gap and gap != "N/A":
            lines.append(f"- 实际年龄：{ca} 岁（年龄差：{gap} 岁）")
        lines.append(f"- 风险等级：{risk}")
        if suggestion:
            lines.append(f"- 建议：{suggestion}")
        lines.append("")

    # 知识引用
    ktext = result.get("knowledge_text", "")
    if ktext:
        lines.append(ktext)

    return "\n".join(lines).strip()


# --------------------------------------------------------------------------
# 路由：历史记录
# --------------------------------------------------------------------------

@app.get("/history/{session_id}")
async def get_history(session_id: str):
    if session_id not in sessions:
        raise HTTPException(status_code=404, detail="Session not found")
    memory = sessions[session_id]
    return {
        "session_id": session_id,
        "turn_count": memory.turn_count,
        "is_compressed": memory.is_compressed,
        "context": memory.get_context(),
    }


# --------------------------------------------------------------------------
# 路由：用户反馈
# --------------------------------------------------------------------------

@app.post("/feedback")
async def feedback(
    session_id: str = Form(...),
    rating: str = Form(...),  # "positive" | "negative"
):
    if session_id not in sessions:
        raise HTTPException(status_code=404, detail="Session not found")
    memory = sessions[session_id]
    ok = memory.set_feedback(rating)
    if not ok:
        raise HTTPException(status_code=400, detail="No turn to rate")
    # 收集统计
    turn = memory.get_last_turn()
    if rating == "positive":
        METRICS["feedback_positive"] += 1
    elif rating == "negative":
        METRICS["feedback_negative"] += 1
    return {
        "ok": True,
        "rating": rating,
        "turn_type": turn.analysis_type if turn else None,
    }


# --------------------------------------------------------------------------
# 路由：监控指标
# --------------------------------------------------------------------------

METRICS = {
    "total_requests": 0,
    "total_duration_ms": 0.0,
    "feedback_positive": 0,
    "feedback_negative": 0,
}


@app.get("/metrics")
async def get_metrics():
    avg_ms = METRICS["total_duration_ms"] / max(METRICS["total_requests"], 1)
    return {
        **METRICS,
        "avg_duration_ms": round(avg_ms, 1),
        "feedback_rate": (
            round(METRICS["feedback_positive"] / max(METRICS["feedback_negative"] + METRICS["feedback_positive"], 1), 2)
            if METRICS["feedback_positive"] + METRICS["feedback_negative"] > 0
            else 0
        ),
    }
