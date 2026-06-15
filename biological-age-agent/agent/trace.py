"""
Trace 记录 —— 追踪 Agent 的每一步决策和执行。

用于：
- 调试 Agent 行为
- 评估（Eval）数据源
- 面试时展示（Trace 可视化）
"""

import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class StepStatus(Enum):
    SUCCESS = "success"
    ERROR = "error"
    TIMEOUT = "timeout"
    DEGRADED = "degraded"


@dataclass
class TraceStep:
    """Agent 执行轨迹中的单个步骤。"""
    step_id: str
    parent_id: Optional[str]
    tool_name: Optional[str]
    input_snapshot: Dict[str, Any]
    output_snapshot: Optional[Any]
    start_time: float
    end_time: float
    status: StepStatus
    thinking: Optional[str] = None
    error_message: Optional[str] = None


class TraceRecorder:
    """Agent 执行轨迹记录器。"""

    def __init__(self):
        self._steps: List[TraceStep] = []
        self._root_id = str(uuid.uuid4())
        self._start_time = time.time()

    @contextmanager
    def step(
        self,
        tool_name: Optional[str] = None,
        thinking: Optional[str] = None,
        parent_id: Optional[str] = None,
    ):
        """记录一个执行步骤。"""
        step_id = str(uuid.uuid4())
        start = time.time()
        input_snapshot = {"tool": tool_name, "thinking": thinking}

        try:
            yield
            self._steps.append(TraceStep(
                step_id=step_id,
                parent_id=parent_id or self._root_id,
                tool_name=tool_name,
                input_snapshot=input_snapshot,
                output_snapshot=None,  # 外部更新
                start_time=start,
                end_time=time.time(),
                status=StepStatus.SUCCESS,
                thinking=thinking,
            ))
        except Exception as e:
            self._steps.append(TraceStep(
                step_id=step_id,
                parent_id=parent_id or self._root_id,
                tool_name=tool_name,
                input_snapshot=input_snapshot,
                output_snapshot=None,
                start_time=start,
                end_time=time.time(),
                status=StepStatus.ERROR,
                thinking=thinking,
                error_message=str(e),
            ))
            raise

    def to_dict(self) -> Dict[str, Any]:
        """导出为字典（用于可视化或 Eval）。"""
        return {
            "total_duration_ms": (time.time() - self._start_time) * 1000,
            "total_steps": len(self._steps),
            "success_count": sum(
                1 for s in self._steps if s.status == StepStatus.SUCCESS
            ),
            "error_count": sum(
                1 for s in self._steps if s.status == StepStatus.ERROR
            ),
            "tool_usage": self._tool_usage(),
            "steps": [
                {
                    "tool": s.tool_name,
                    "status": s.status.value,
                    "duration_ms": (s.end_time - s.start_time) * 1000,
                    "error": s.error_message,
                }
                for s in self._steps
            ],
        }

    def _tool_usage(self) -> Dict[str, int]:
        """统计各工具调用次数。"""
        counts: Dict[str, int] = {}
        for s in self._steps:
            if s.tool_name:
                counts[s.tool_name] = counts.get(s.tool_name, 0) + 1
        return counts
