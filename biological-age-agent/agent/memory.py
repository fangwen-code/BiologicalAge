"""
分层记忆系统。

- 短期记忆：保留最近 5 轮原始对话
- 压缩记忆：超过 5 轮后压缩为摘要
- 用户档案：持久化存储，不因压缩丢失
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class UserProfile:
    """用户持久化档案。"""
    user_id: str
    baseline_age: Optional[float] = None
    history_analyses: List[Dict[str, Any]] = field(default_factory=list)
    preferences: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ConversationTurn:
    """单轮对话记录。"""
    user_input: str
    agent_output: str
    timestamp: float
    analysis_type: Optional[str] = None
    feedback: Optional[str] = None  # "positive" / "negative" / None


class ConversationMemory:
    """分层对话记忆。

    使用方式:
        memory = ConversationMemory(max_raw_rounds=5)
        memory.add_turn(user_msg, agent_msg)

        # 获取提供给 LLM 的上下文
        context = memory.get_context()
    """

    def __init__(self, max_raw_rounds: int = 5):
        self.max_raw_rounds = max_raw_rounds
        self._raw_turns: List[ConversationTurn] = []
        self._summary: str = ""
        self.profile: Optional[UserProfile] = None

    def add_turn(
        self,
        user_input: str,
        agent_output: str,
        analysis_type: Optional[str] = None,
    ) -> None:
        """添加一轮对话，超限时自动压缩。"""
        import time

        self._raw_turns.append(ConversationTurn(
            user_input=user_input,
            agent_output=agent_output,
            timestamp=time.time(),
            analysis_type=analysis_type,
        ))

        if len(self._raw_turns) >= self.max_raw_rounds:
            self._compress()

    def _compress(self) -> None:
        """将原始对话压缩为摘要。"""
        # 当前用简单截断策略，后续可接入 LLM 压缩
        compressed = []
        for turn in self._raw_turns:
            compressed.append(
                f"User: {turn.user_input[:100]} | "
                f"Agent: {turn.agent_output[:200]}"
            )

        new_summary = "\n".join(compressed)

        if self._summary:
            self._summary = f"{self._summary}\n{new_summary}"
        else:
            self._summary = new_summary

        self._raw_turns.clear()

    def get_context(self) -> List[Dict[str, str]]:
        """获取供 LLM 使用的上下文消息列表。"""
        context = []

        if self._summary:
            context.append({
                "role": "system",
                "content": f"History summary:\n{self._summary}",
            })

        for turn in self._raw_turns:
            context.append({"role": "user", "content": turn.user_input})
            context.append({"role": "assistant", "content": turn.agent_output})

        return context

    @property
    def turn_count(self) -> int:
        """当前原始对话轮数。"""
        return len(self._raw_turns)

    @property
    def is_compressed(self) -> bool:
        """是否已发生压缩。"""
        return bool(self._summary)

    def get_last_turn(self) -> Optional[ConversationTurn]:
        """获取最后一轮原始对话（用于反馈绑定）。"""
        if self._raw_turns:
            return self._raw_turns[-1]
        return None

    def set_feedback(self, rating: str) -> bool:
        """给最后一轮对话打反馈标记。"""
        turn = self.get_last_turn()
        if turn is None:
            return False
        if rating not in ("positive", "negative"):
            return False
        turn.feedback = rating
        return True
