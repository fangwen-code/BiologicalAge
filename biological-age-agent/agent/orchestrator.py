"""
Agent Orchestrator —— 智能调度核心。

职责：理解用户意图 → 决定调用哪些工具 → 调度执行 → 输出结果。
采用单 Agent 架构，仅在最必要时调用 LLM。
"""

import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional

import pandas as pd

from app.model import predict_bio_age
from agent.knowledge_tool import format_knowledge_reply, search_knowledge_by_risk
from agent.llm_client import get_llm_client


# ──────────────────────────────────────────────────────────────────────
# 错误分类 —— 将原始异常翻译为用户可操作的中文指引
# ──────────────────────────────────────────────────────────────────────

def _classify_prediction_error(exc: Exception) -> str:
    """分析预测异常，返回面向用户的中文修正建议。"""
    msg = str(exc)

    # 模型加载失败
    if "Model file not found" in msg:
        return "⚠️ 模型文件缺失，请联系管理员恢复 model/bio_age_prediction_model.pkl。"

    if "Failed to load model" in msg:
        return "⚠️ 模型初始化失败，依赖版本可能不兼容。尝试：pip install scikit-learn==1.7.2"

    # 文件格式问题
    if "Cannot detect valid CpG columns" in msg:
        return (
            "⚠️ **未检测到有效的 CpG 位点列**\n\n"
            "文件需要包含 hg38 chr_pos 格式的位点（如 chr1_123456）。\n\n"
            "请先下载 Demo 文件参考格式：首页 → 第一步 → 下载 Demo 文件"
        )

    if "cg" in msg and "prefix" in msg.lower():
        return (
            "⚠️ **检测到 cg 开头的探针**\n\n"
            "需要先转换为 hg38 chr_pos 格式再上传。\n"
            "可使用 UCSC Genome Browser 或 R 包（ChAMP / minfi）进行转换。"
        )

    # 性别列问题
    if "Gender" in msg:
        return (
            "⚠️ **性别列格式有误**\n\n"
            "支持的格式：女/男、0/1、Female/Male、F/M。\n"
            "请检查文件中 Gender 列的值。"
        )

    # 重复样本
    if "Duplicate SampleID" in msg:
        return (
            "⚠️ **发现重复的样本编号**\n\n"
            "请确保每个 SampleID 唯一，删除重复行后重试。"
        )

    # 缺失位点（fill_missing=False 时）
    if "Missing" in msg and "CpG" in msg:
        return "⚠️ 数据缺失大量 CpG 位点，请补充后重试。"

    # sklearn / numpy 版本不兼容
    if "scalar arrays" in msg or "InconsistentVersion" in msg:
        return (
            "⚠️ **模型版本兼容性问题**\n\n"
            "模型由旧版 sklearn 训练，当前环境可能不兼容。\n"
            "尝试：pip install scikit-learn==1.7.2"
        )

    # 通用兜底
    return f"⚠️ 分析异常：{msg}"


class AnalysisMode(Enum):
    """分析模式，由输入数据决定。"""
    METHYLATION = "meth"        # 有甲基化数据
    KNOWLEDGE = "knowledge"     # 无数据，纯知识问答


@dataclass
class ToolResult:
    """单个工具的执行结果。"""
    tool_name: str
    success: bool
    data: Any = None
    error: Optional[str] = None
    duration_ms: float = 0.0


@dataclass
class AgentContext:
    """Agent 执行上下文。"""
    mode: AnalysisMode
    user_input: str
    has_methylation: bool = False
    tool_results: Dict[str, ToolResult] = field(default_factory=dict)
    memory: Optional[Any] = None


class AgentOrchestrator:
    """Agent 主调度器。"""

    def __init__(self):
        self.llm = None  # 仅在需要时初始化

    async def run(
        self,
        user_input: str,
        methylation_data: Optional[pd.DataFrame] = None,
        memory: Optional[Any] = None,
        existing_results: Optional[list] = None,
    ) -> Dict[str, Any]:
        """执行一次完整的 Agent 分析流程。

        Args:
            existing_results: 经典版预测结果（带 Predicted_Age 列），
                              跳过模型预测直接进入解读。
        """
        if existing_results:
            return await self._interpret_existing_results(user_input, existing_results)

        # 1. 确定模式
        ctx = AgentContext(
            mode=self._determine_mode(user_input, methylation_data),
            user_input=user_input,
            has_methylation=methylation_data is not None,
            memory=memory,
        )

        # 2. 根据模式执行
        try:
            if ctx.mode == AnalysisMode.METHYLATION:
                return await self._run_methylation_pipeline(ctx, methylation_data)
            else:
                return await self._run_knowledge_pipeline(ctx)
        except Exception as e:
            friendly = _classify_prediction_error(e)
            return {
                "type": "prediction_error",
                "mode": ctx.mode.value,
                "error": str(e),
                "message": friendly,
            }

    async def _interpret_existing_results(
        self,
        user_input: str,
        results: list,
    ) -> Dict[str, Any]:
        """解读经典版预测结果，跳过模型预测直接到知识检索+报告生成。"""
        analysis_data = {
            "mode": "result_interpretation",
            "tool_results": {"methylation_clock": results},
            "has_methylation": True,
        }
        return await self._generate_report_from_data(analysis_data, user_input)

    async def _generate_report_from_data(
        self,
        analysis_data: dict,
        user_input: str,
    ) -> Dict[str, Any]:
        """基于已有的分析数据生成报告（供解读模式使用）。"""
        from agent.knowledge_tool import format_knowledge_reply, search_knowledge_by_risk
        from agent.llm_client import get_llm_client

        meth_data = analysis_data.get("tool_results", {}).get("methylation_clock")
        if meth_data and isinstance(meth_data, list) and len(meth_data) > 0:
            row = meth_data[0]
            risk = row.get("Risk_Level", "")
            gap_str = row.get("Age_Gap", "")
            gap = float(gap_str) if isinstance(gap_str, (int, float)) else 0.0
            entries = search_knowledge_by_risk(gap, risk)
            if entries:
                analysis_data["knowledge_text"] = format_knowledge_reply(entries)

        llm = get_llm_client()
        if llm.is_available():
            try:
                llm_report = await llm.generate_report(
                    analysis_data=analysis_data, user_input=user_input,
                )
                analysis_data["llm_report"] = llm_report
                analysis_data["report_source"] = "llm"
            except Exception:
                analysis_data["report_source"] = "template"
        else:
            analysis_data["report_source"] = "template"

        return analysis_data

    def _determine_mode(
        self,
        user_input: str,
        methylation_data: Optional[pd.DataFrame],
    ) -> AnalysisMode:
        """根据输入数据确定分析模式。"""
        has_meth = methylation_data is not None and len(methylation_data) > 0
        return AnalysisMode.METHYLATION if has_meth else AnalysisMode.KNOWLEDGE

    async def _run_methylation_pipeline(
        self,
        ctx: AgentContext,
        methylation_data: pd.DataFrame,
    ) -> Dict[str, Any]:
        """甲基化数据的分析流程。"""
        try:
            meth_result = predict_bio_age(methylation_data)
            ctx.tool_results["methylation_clock"] = ToolResult(
                tool_name="methylation_clock", success=True, data=meth_result
            )
            return await self._generate_report(ctx)
        except Exception as e:
            friendly = _classify_prediction_error(e)
            return {
                "type": "prediction_error",
                "mode": "meth",
                "error": str(e),
                "message": friendly,
            }

    async def _run_knowledge_pipeline(
        self,
        ctx: AgentContext,
    ) -> Dict[str, Any]:
        """无数据时的处理：有历史则追问解读，无历史则引导。"""
        # 检查 memory 中是否有历史预测结果
        if ctx.memory and (ctx.memory.turn_count > 0 or ctx.memory.is_compressed):
            # 有历史对话，是追问场景
            context = ctx.memory.get_context()
            llm = get_llm_client()
            if llm.is_available():
                try:
                    memory_text = "\n".join(
                        f"{m['role']}: {m['content'][:200]}"
                        for m in context if isinstance(m, dict)
                    )
                    answer = await llm.chat([
                        {
                            "role": "system",
                            "content": "你是生物学年龄分析助手。用户正在就之前的分析结果提问。"
                            "基于以下对话历史回答用户的问题。回答要简洁准确。",
                        },
                        *[m for m in context if isinstance(m, dict)],
                        {"role": "user", "content": ctx.user_input},
                    ])
                    return {"type": "follow_up", "llm_answer": answer, "mode": "knowledge"}
                except Exception:
                    pass

            # LLM 不可用或失败时，基于记忆做简单回复
            return {
                "type": "follow_up",
                "message": "请具体描述你想了解的内容，比如："
                           "'解释一下样本A的风险' 或 '有什么干预建议'。",
                "mode": "knowledge",
            }

        # 无历史，首次咨询 → 引导上传
        return {
            "type": "guidance",
            "message": "请上传甲基化数据以开始分析。",
            "options": [
                {
                    "type": "methylation",
                    "description": "上传甲基化数据 → 精确的生物学年龄预测",
                },
            ],
            "mode": "knowledge",
        }

    def _build_analysis_data(self, ctx: AgentContext) -> dict:
        """构建供 LLM 使用的结构化分析数据。"""
        data = {
            "mode": ctx.mode.value,
            "has_methylation": ctx.has_methylation,
            "tool_results": {},
        }

        for name, tr in ctx.tool_results.items():
            if isinstance(tr.data, pd.DataFrame):
                data["tool_results"][name] = tr.data.to_dict(orient="records")
            else:
                data["tool_results"][name] = tr.data

        return data

    async def _generate_report(self, ctx: AgentContext) -> Dict[str, Any]:
        """生成最终报告。

        流程：构建数据 → 检索知识库 → LLM报告（可选）→ 返回。
        """
        analysis_data = self._build_analysis_data(ctx)

        # 1. 从预测结果中提取年龄加速和风险等级，检索知识库
        knowledge_text = ""
        tool_data = analysis_data.get("tool_results", {})
        meth_data = tool_data.get("methylation_clock")
        if meth_data and isinstance(meth_data, list) and len(meth_data) > 0:
            row = meth_data[0]
            risk = row.get("Risk_Level", "")
            gap_str = row.get("Age_Gap", "")
            gap = float(gap_str) if isinstance(gap_str, (int, float)) else 0.0
            entries = search_knowledge_by_risk(gap, risk)
            if entries:
                knowledge_text = format_knowledge_reply(entries)
                analysis_data["knowledge_entries"] = entries
                analysis_data["knowledge_text"] = knowledge_text

        # 2. 尝试 LLM 报告
        llm = get_llm_client()
        if llm.is_available():
            try:
                llm_report = await llm.generate_report(
                    analysis_data=analysis_data,
                    user_input=ctx.user_input,
                )
                analysis_data["llm_report"] = llm_report
                analysis_data["report_source"] = "llm"
            except Exception as e:
                print(f"LLM report failed, using template: {e}")
                analysis_data["report_source"] = "template"
        else:
            analysis_data["report_source"] = "template"

        return analysis_data
