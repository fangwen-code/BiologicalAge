"""
自定义 LLM 客户端 —— 支持任意 OpenAI 兼容 API。

通过环境变量配置：
    LLM_API_KEY      你的 API Key（默认从环境变量读取）
    LLM_BASE_URL     API 地址（默认 https://api.deepseek.com/v1）
    LLM_MODEL        模型名（默认 deepseek-chat）

用法:
    client = LLMClient()
    report = await client.generate_report(analysis_data)
"""

import os
import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import httpx


@dataclass
class LLMConfig:
    """LLM 客户端配置。

    优先级：构造函数参数 > 环境变量 > 默认值
    """
    api_key: str = ""
    base_url: str = "https://api.deepseek.com/v1"
    model: str = "deepseek-chat"
    temperature: float = 0.0
    max_tokens: int = 2048
    timeout_sec: int = 30

    def __post_init__(self):
        if not self.api_key:
            self.api_key = os.getenv("LLM_API_KEY", "")
        self.base_url = os.getenv("LLM_BASE_URL", self.base_url)
        self.model = os.getenv("LLM_MODEL", self.model)


class LLMClient:
    """通用 LLM 客户端（OpenAI 兼容接口）。"""

    def __init__(self, config: Optional[LLMConfig] = None):
        self.config = config or LLMConfig()
        self._client: Optional[httpx.AsyncClient] = None

    @property
    def client(self) -> httpx.AsyncClient:
        if self._client is None:
            self._client = httpx.AsyncClient(
                base_url=self.config.base_url,
                timeout=self.config.timeout_sec,
                headers={
                    "Authorization": f"Bearer {self.config.api_key}",
                    "Content-Type": "application/json",
                },
            )
        return self._client

    def is_available(self) -> bool:
        """检查是否配置了 API Key。"""
        return bool(self.config.api_key)

    async def chat(
        self,
        messages: List[Dict[str, str]],
        **kwargs,
    ) -> str:
        """通用聊天补全。

        Args:
            messages: 消息列表，如 [{"role": "user", "content": "..."}]
            **kwargs: 覆盖 config 的参数（temperature, max_tokens 等）

        Returns:
            模型回复文本。

        Raises:
            ConnectionError: API 调用失败。
        """
        payload = {
            "model": kwargs.get("model", self.config.model),
            "messages": messages,
            "temperature": kwargs.get("temperature", self.config.temperature),
            "max_tokens": kwargs.get("max_tokens", self.config.max_tokens),
            "stream": False,
        }

        try:
            resp = await self.client.post("/chat/completions", json=payload)
            resp.raise_for_status()
            data = resp.json()
            return data["choices"][0]["message"]["content"]
        except httpx.HTTPStatusError as e:
            raise ConnectionError(
                f"LLM API error ({e.response.status_code}): "
                f"{e.response.text[:200]}"
            )
        except Exception as e:
            raise ConnectionError(f"LLM request failed: {str(e)}")

    async def generate_report(
        self,
        analysis_data: Dict[str, Any],
        user_input: str = "",
    ) -> str:
        """根据分析数据生成可读报告。

        Args:
            analysis_data: Agent 分析结果（含预测数据、模式、上下文）。
            user_input: 用户的原始输入。

        Returns:
            自然语言报告。
        """
        system_prompt = """你是专业的衰老科学分析助手。你的任务是根据用户的甲基化数据分
析结果，生成一份可读、专业、个性化的生物学年龄评估报告。

要求：
1. 报告结构清晰，包含：核心结论、数据说明、风险解读、建议
2. 语言自然友好，但保持科学严谨
3. 不要编造数据，严格基于给定的分析结果
4. 如果数据不完整，明确告知局限性
5. 报告中提及样本时，必须使用数据中提供的原始 Sample_ID，不能简化为"样本1、样本2"等
6. 如果下方提供了「参考资料」，请选择性引用其结论来支撑报告内容（引用标注 [1][2] 等）
7. 不要编造参考文献，只引用下方提供的内容。报告长度控制在 300-500 字"""

        # 构造数据摘要
        data_summary = json.dumps(
            {k: v for k, v in analysis_data.items() if k != "knowledge_entries"},
            ensure_ascii=False, indent=2,
        )

        import datetime
        today = datetime.date.today().strftime("%Y-%m-%d")

        knowledge_text = analysis_data.get("knowledge_text", "")
        ref_section = f"\n\n参考资料：\n{knowledge_text}" if knowledge_text else ""

        user_prompt = f"""用户输入：{user_input or '（无）'}
报告日期：{today}

分析数据：
{data_summary}{ref_section}

请生成一份完整的生物学年龄评估报告。
注意：
1. 报告日期必须使用上面提供的「报告日期」
2. 每个样本必须使用 Sample_ID 原始名称，不得简化为"样本1、样本2"等
3. 报告开头注明「报告生成于 {today}」"""

        return await self.chat([
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ])

    async def close(self):
        if self._client:
            await self._client.aclose()
            self._client = None


# --------------------------------------------------------------------------
# 全局单例（懒加载）
# --------------------------------------------------------------------------

_global_client: Optional[LLMClient] = None


def get_llm_client() -> LLMClient:
    """获取全局 LLM 客户端（懒加载）。"""
    global _global_client
    if _global_client is None:
        _global_client = LLMClient()
    return _global_client
