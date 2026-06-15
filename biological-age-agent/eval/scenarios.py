"""
Eval 测试场景定义。

覆盖 4 类场景：正常流程、边界情况、异常恢复、记忆边界。
全部使用合成数据，不依赖真实数据。
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class TestScenario:
    """单个测试场景。"""
    name: str
    category: str  # happy_path / boundary / error / consistency / memory
    input: Dict[str, Any]
    expected: Dict[str, Any]
    description: Optional[str] = None
    mock: Optional[Dict[str, Any]] = None  # 模拟外部服务


@dataclass
class EvalResult:
    """单个场景的评估结果。"""
    scenario_name: str
    passed: bool
    details: Dict[str, Any]
    error: Optional[str] = None


# ======================================================================
# 场景定义
# ======================================================================

SCENARIOS: List[TestScenario] = []


def _register(scenario: TestScenario):
    SCENARIOS.append(scenario)
    return scenario


# --- 1. Happy Path ---

_register(TestScenario(
    name="methylation_prediction",
    category="happy_path",
    description="正常预测：上传甲基化数据，Agent 正确路由到预测流程",
    input={"has_methylation": True},
    expected={
        "mode": "meth",
        "graceful_on_error": True,
    },
))

_register(TestScenario(
    name="no_data_guidance",
    category="happy_path",
    description="无数据时输出引导信息",
    input={"has_methylation": False},
    expected={
        "mode": "knowledge",
        "output_contains_options": True,
    },
))

# --- 2. Boundary ---

_register(TestScenario(
    name="empty_methylation_file",
    category="boundary",
    description="甲基化数据为空文件，系统不崩溃",
    input={"has_methylation": True, "empty_file": True},
    expected={
        "behavior": "graceful_error",
    },
))

_register(TestScenario(
    name="extreme_age_prediction",
    category="boundary",
    description="极端高龄预测（90+），系统不崩溃",
    input={"has_methylation": True, "target_age": 95},
    expected={
        "behavior": "graceful_error",
    },
))

# --- 3. Error ---

_register(TestScenario(
    name="malformed_data_format",
    category="error",
    description="上传文件格式错误",
    input={"has_methylation": True, "corrupt_data": True},
    expected={
        "behavior": "graceful_error",
        "should_not_crash": True,
    },
))

_register(TestScenario(
    name="missing_model_file",
    category="error",
    description="模型文件丢失时的错误处理",
    input={"has_methylation": True},
    mock={"model_missing": True},
    expected={
        "behavior": "graceful_error",
        "should_not_crash": True,
    },
))

# --- 4. Memory ---

_register(TestScenario(
    name="memory_compression",
    category="memory",
    description="6 轮对话后自动触发压缩",
    input={"has_methylation": False, "turns": 6},
    expected={
        "memory_compressed": True,
    },
))

_register(TestScenario(
    name="memory_across_sessions",
    category="memory",
    description="跨 session 记忆保留用户档案",
    input={"has_methylation": False, "session_id": "repeat_user"},
    expected={
        "profile_loaded": True,
    },
))
