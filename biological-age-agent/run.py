"""
启动入口。

用法:
    python run.py

环境变量配置（可选，用于 LLM 报告生成）:
    LLM_API_KEY    你的 API Key
    LLM_BASE_URL   API 地址（默认 https://api.deepseek.com/v1）
    LLM_MODEL      模型名（默认 deepseek-chat）

也可以创建 .env 文件（参考 .env.example）。
"""

import os
from pathlib import Path

# 加载 .env 文件（如果存在）
dotenv_path = Path(__file__).parent / ".env"
if dotenv_path.exists():
    try:
        from dotenv import load_dotenv
        load_dotenv(dotenv_path)
        print(f"Loaded environment from {dotenv_path}")
    except ImportError:
        pass

import uvicorn

if __name__ == "__main__":
    # 启动提示
    api_key = os.getenv("LLM_API_KEY", "")
    if api_key:
        print("LLM report generation: enabled")
    else:
        print("LLM report generation: disabled (set LLM_API_KEY to enable)")

    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        log_level="info",
    )
