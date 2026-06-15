"""
单元测试 —— 覆盖核心模块，不依赖模型文件。
"""

import io
import sys
from pathlib import Path

import pandas as pd
import pytest

# 确保项目根目录在 path 中
_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_root))

from app.utils import read_input_file, generate_demo
from agent.knowledge_tool import search_knowledge, search_knowledge_by_risk, format_knowledge_reply
from agent.memory import ConversationMemory


# ======================================================================
# Utils — 文件解析
# ======================================================================

class TestReadInputFile:
    def test_xlsx(self):
        """解析 xlsx 文件。"""
        df = pd.DataFrame({"SampleID": ["S1"], "Gender": ["女"]})
        buf = io.BytesIO()
        df.to_excel(buf, index=False, engine="openpyxl")
        result = read_input_file(buf.getvalue(), "test.xlsx")
        assert result["SampleID"].iloc[0] == "S1"

    def test_csv(self):
        """解析 csv 文件。"""
        content = "SampleID,Gender\nS1,女\n".encode("utf-8")
        result = read_input_file(content, "test.csv")
        assert result["SampleID"].iloc[0] == "S1"

    def test_tsv(self):
        """解析 tsv 文件。"""
        content = "SampleID\tGender\nS1\t女\n".encode("utf-8")
        result = read_input_file(content, "test.tsv")
        assert result["SampleID"].iloc[0] == "S1"

    def test_txt(self):
        """解析 txt（制表符分隔）文件。"""
        content = "SampleID\tGender\nS1\t女\n".encode("utf-8")
        result = read_input_file(content, "test.txt")
        assert result["SampleID"].iloc[0] == "S1"

    def test_unsupported_format(self):
        """不支持的格式应抛 ValueError。"""
        with pytest.raises(ValueError, match="Unsupported"):
            read_input_file(b"", "test.pdf")


class TestGenerateDemo:
    def test_generates_valid_excel(self):
        """Demo 文件应生成有效的 Excel 字节流。"""
        cpgs = ["chr1_100", "chr2_200", "chr3_300"]
        buf = generate_demo(cpgs, n_samples=2)
        df = pd.read_excel(buf, engine="openpyxl")
        assert len(df) == 2
        assert "SampleID" in df.columns
        assert "Gender" in df.columns
        assert "Note" in df.columns

    def test_includes_cpg_columns(self):
        """Demo 文件应包含传入的 CpG 位点列。"""
        cpgs = ["chr1_100", "chr2_200"]
        buf = generate_demo(cpgs, n_samples=1)
        df = pd.read_excel(buf, engine="openpyxl")
        for cpg in cpgs:
            assert cpg in df.columns


# ======================================================================
# Knowledge Tool — 知识检索
# ======================================================================

class TestKnowledgeSearch:
    def test_search_by_keyword(self):
        """关键词检索应返回匹配条目。"""
        results = search_knowledge("吸烟", max_results=3)
        assert len(results) > 0
        assert any("吸烟" in str(r.content) for r in results)

    def test_search_empty_query(self):
        """空查询应返回空列表。"""
        results = search_knowledge("")
        assert len(results) == 0

    def test_search_by_risk_lower(self):
        """低风险应返回生活方式相关条目。"""
        results = search_knowledge_by_risk(-2.0, "lower risk", max_results=2)
        assert len(results) > 0

    def test_search_by_risk_higher(self):
        """高风险应返回风险因素条目。"""
        results = search_knowledge_by_risk(6.0, "higher risk", max_results=2)
        assert len(results) > 0

    def test_format_reply(self):
        """格式化函数应返回带引用的文本。"""
        entries = search_knowledge("吸烟", max_results=1)
        formatted = format_knowledge_reply(entries)
        assert "吸烟" in formatted or "Smoking" in formatted or formatted == ""

    def test_format_empty(self):
        """空条目应返回空字符串。"""
        assert format_knowledge_reply([]) == ""


# ======================================================================
# Memory — 对话记忆
# ======================================================================

class TestConversationMemory:
    def test_add_turn(self):
        """添加一轮对话后 turn_count 应增加。"""
        mem = ConversationMemory()
        mem.add_turn("你好", "你好！")
        assert mem.turn_count == 1

    def test_compression(self):
        """超过 max_raw_rounds 后应触发压缩。"""
        mem = ConversationMemory(max_raw_rounds=3)
        for i in range(4):
            mem.add_turn(f"问题{i}", f"回答{i}")
        assert mem.is_compressed

    def test_context_includes_summary_after_compression(self):
        """压缩后 context 应包含摘要。"""
        mem = ConversationMemory(max_raw_rounds=2)
        mem.add_turn("问题1", "回答1")
        mem.add_turn("问题2", "回答2")
        mem.add_turn("问题3", "回答3")  # 触发压缩
        ctx = mem.get_context()
        assert any("History summary" in str(m.get("content", "")) for m in ctx)

    def test_feedback(self):
        """反馈标记应更新最后一轮。"""
        mem = ConversationMemory()
        mem.add_turn("你好", "你好！")
        assert mem.set_feedback("positive") is True
        assert mem.get_last_turn().feedback == "positive"

    def test_feedback_no_turns(self):
        """无对话时设置反馈应失败。"""
        mem = ConversationMemory()
        assert mem.set_feedback("positive") is False

    def test_feedback_invalid_rating(self):
        """无效反馈值应失败。"""
        mem = ConversationMemory()
        mem.add_turn("你好", "你好！")
        assert mem.set_feedback("invalid") is False
