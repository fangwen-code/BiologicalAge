"""
MCP Server —— 将生物学年龄模型封装为标准 MCP Tool。

任何 MCP Client（Claude Desktop、Cursor 等）均可发现和调用。
"""

from pathlib import Path
from typing import Optional

import pandas as pd

from app.model import load_model, predict_bio_age

# 尝试导入 fastmcp（如果环境中安装了的话）
try:
    from fastmcp import FastMCP

    mcp = FastMCP(name="biological-age-predictor")

    @mcp.tool()
    def predict_methylation_age(
        file_path: str,
        fill_missing: bool = True,
    ) -> str:
        """基于甲基化数据预测生物学年龄。

        Args:
            file_path: 甲基化数据文件路径（xlsx/csv/tsv/txt）。
            fill_missing: 缺失位点时是否用均值填充。

        Returns:
            预测结果的 JSON 字符串。
        """
        ext = Path(file_path).suffix.lower()
        if ext in (".xlsx", ".xls"):
            df = pd.read_excel(file_path, engine="openpyxl")
        elif ext == ".csv":
            df = pd.read_csv(file_path)
        elif ext in (".tsv", ".txt"):
            df = pd.read_csv(file_path, sep="\t")
        else:
            return f"Unsupported format: {ext}"

        result = predict_bio_age(df, fill_missing=fill_missing)
        return result.to_json(orient="records", force_ascii=False)

    def run():
        """启动 MCP Server。"""
        load_model()
        mcp.run()

except ImportError:
    import warnings

    warnings.warn("fastmcp not installed. MCP server unavailable.")

    def run():
        print("fastmcp is required to run MCP server.")
        print("Install: pip install fastmcp")
