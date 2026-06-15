"""
合成数据生成器 —— 用于 Agent 行为评估，不需要真实数据。
"""

from typing import List

import numpy as np
import pandas as pd


def generate_synthetic_methylation(
    seed: int = 42,
    target_age: float = 45.0,
    n_sites: int = 100,
) -> pd.DataFrame:
    """生成模拟甲基化 Beta 值数据。

    Args:
        seed: 随机种子。
        target_age: 期望的预测年龄（通过偏移模拟）。
        n_sites: 生成的 CpG 位点数量。

    Returns:
        包含 SampleID、CpG 位点、Gender 的 DataFrame。
    """
    rng = np.random.RandomState(seed)

    sites = [f"chr{r}_{pos}" for r in range(1, 6) for pos in range(rng.randint(1000, 9999, 20))]
    sites = sites[:n_sites]

    data = {"SampleID": [f"synth_{seed}_1"], "Gender": ["女"]}
    for site in sites:
        # 基线 beta 值 + 微小的年龄相关偏移
        base = rng.beta(2, 2)
        age_shift = (target_age - 40) * 0.002 * rng.randn()
        data[site] = [round(np.clip(base + age_shift, 0, 1), 4)]

    df = pd.DataFrame(data)
    df["Chronological_Age"] = target_age
    return df


def generate_synthetic_batch(
    n_samples: int = 5,
    seed: int = 0,
) -> pd.DataFrame:
    """批量生成合成甲基化数据。"""
    dfs = []
    for i in range(n_samples):
        age = 30 + i * 10
        df = generate_synthetic_methylation(seed=seed + i, target_age=age)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)
