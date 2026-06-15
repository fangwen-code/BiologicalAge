"""
工具函数模块。
"""

import io
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

SUPPORTED_FORMATS = {".xlsx", ".xls", ".csv", ".tsv", ".txt"}


def read_input_file(content: bytes, filename: str) -> pd.DataFrame:
    """根据文件扩展名自动选择解析方式。

    Args:
        content: 上传文件的二进制内容。
        filename: 原始文件名（用于判断扩展名）。

    Returns:
        解析后的 DataFrame。

    Raises:
        ValueError: 不支持的格式。
    """
    ext = Path(filename).suffix.lower()

    if ext in (".xlsx", ".xls"):
        return pd.read_excel(io.BytesIO(content), engine="openpyxl")
    elif ext == ".csv":
        return pd.read_csv(io.BytesIO(content))
    elif ext in (".tsv", ".txt"):
        return pd.read_csv(io.BytesIO(content), sep="\t")
    else:
        raise ValueError(
            f"Unsupported file format: {ext}. "
            f"Supported: {', '.join(sorted(SUPPORTED_FORMATS))}"
        )


def generate_demo(train_cpg_names: List[str], n_samples: int = 3) -> io.BytesIO:
    """生成示例数据文件，供用户下载参考格式。

    Args:
        train_cpg_names: 模型使用的所有 CpG 位点名称。
        n_samples: 生成的样本数量。

    Returns:
        包含示例数据的 Excel 文件的 BytesIO 对象。
    """
    np.random.seed(42)
    n_demo = min(5, len(train_cpg_names))
    demo_cpgs = np.random.choice(train_cpg_names, size=n_demo, replace=False).tolist()

    demo_data = {
        "SampleID": [f"Demo_Sample_{i+1:02d}" for i in range(n_samples)],
        **{
            cpg: [
                round(np.random.uniform(0.1, 0.9), 3)
                for _ in range(n_samples)
            ]
            for cpg in demo_cpgs
        },
        "Gender": ["女", 1, "男"][:n_samples] + ["女"] * max(0, n_samples - 3),
        "Chronological_Age": [32.5, 48.7, 65.2][:n_samples] + [40.0] * max(0, n_samples - 3),
    }

    column_names = ["SampleID"] + demo_cpgs + ["Gender", "Chronological_Age"]
    demo_df = pd.DataFrame(demo_data)[column_names]

    note = (
        "format: hg38 chr_pos, such as chr1_1234. "
        "CpGs starting with cg need to convert to hg38 chr_pos format; "
        "beta value: 0-1; "
        "Gender: 女/男, 0/1, Female/Male, F/M"
    )
    demo_df["Note"] = [note] * n_samples

    output = io.BytesIO()
    demo_df.to_excel(output, index=False, engine="openpyxl")
    output.seek(0)
    return output
