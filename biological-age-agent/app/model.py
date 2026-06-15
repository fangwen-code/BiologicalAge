"""
模型加载与预测逻辑模块。

封装甲基化数据的验证、预处理和生物学年龄预测流程。
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import joblib
import numpy as np
import pandas as pd
from fastapi import HTTPException

# --------------------------------------------------------------------------
# 路径
# --------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
MODEL_PATH = BASE_DIR / "model" / "bio_age_prediction_model.pkl"

# --------------------------------------------------------------------------
# 全局模型状态
# --------------------------------------------------------------------------

scaler: Any = None
model: Any = None
train_cpg_names: List[str] = []
train_features_mean: pd.Series = pd.Series(dtype=float)

# --------------------------------------------------------------------------
# 常量
# --------------------------------------------------------------------------

CHR_POS_PATTERN: re.Pattern = re.compile(r"^chr[0-9XYM]+_\d+$")

GENDER_MAP: Dict[Any, int] = {
    "女": 0, "男": 1,
    0: 0, 1: 1,
    "Female": 0, "Male": 1,
    "F": 0, "M": 1,
}


# --------------------------------------------------------------------------
# 模型加载
# --------------------------------------------------------------------------

def load_model() -> None:
    """加载模型包（含 model、scaler、CpG 位点名称、训练集均值）。"""
    global scaler, model, train_cpg_names, train_features_mean

    if scaler is not None:
        return

    if not MODEL_PATH.exists():
        raise HTTPException(
            status_code=500,
            detail=f"Model file not found: {MODEL_PATH}",
        )

    try:
        pkg = joblib.load(MODEL_PATH)
        model = pkg["model"]
        scaler = pkg["scaler"]
        train_cpg_names = list(pkg["train_cpg_names"])
        raw_mean = pkg["train_features_mean"]
        train_features_mean = raw_mean.iloc[:len(train_cpg_names)]
        train_features_mean.index = train_cpg_names
        print(f"Model loaded. Features: {len(train_cpg_names)} CpGs")
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to load model: {str(e)}",
        )


# --------------------------------------------------------------------------
# 数据校验
# --------------------------------------------------------------------------

def _detect_cpg_columns(df: pd.DataFrame) -> List[str]:
    """检测 DataFrame 中的 CpG 位点列（hg38 chr_pos 格式）。

    Raises:
        ValueError: 未检测到有效 CpG 列，或检测到 cg 开头的列。
    """
    chr_pos_cols = [
        col for col in df.columns
        if bool(CHR_POS_PATTERN.match(str(col)))
    ]
    cg_cols = [
        col for col in df.columns
        if isinstance(col, str) and col.startswith("cg")
    ]

    if len(chr_pos_cols) == 0 and len(cg_cols) == 0:
        raise ValueError(
            "Cannot detect valid CpG columns! "
            "Format must be hg38 chr_pos (e.g., chr1_1234). "
            "Download the demo file for reference."
        )
    if len(cg_cols) > 0:
        raise ValueError(
            f"Found {len(cg_cols)} CpGs with 'cg' prefix. "
            "Must convert to hg38 chr_pos format before prediction."
        )
    return chr_pos_cols


def _fill_missing_cpgs(
    df: pd.DataFrame,
    chr_pos_cols: List[str],
    fill_missing: bool,
) -> pd.DataFrame:
    """处理缺失的 CpG 位点。

    Raises:
        ValueError: fill_missing=False 且存在缺失时。
    """
    missing = [cpg for cpg in train_cpg_names if cpg not in chr_pos_cols]

    if len(missing) == 0:
        return df

    if fill_missing:
        print(f"Missing {len(missing)} CpG sites, filling with means.")
        fill_data = {cpg: train_features_mean[cpg] for cpg in missing}
        fill_df = pd.DataFrame(fill_data, index=df.index)
        return pd.concat([df, fill_df], axis=1)

    raise ValueError(
        f"Missing {len(missing)} required CpG sites "
        f"(fill_missing=False):\n{','.join(missing)}"
    )


def _validate_gender(df: pd.DataFrame) -> pd.DataFrame:
    """校验并映射性别列。

    Raises:
        ValueError: 存在无法识别的性别值。
    """
    unknown = df.loc[~df["Gender"].isin(GENDER_MAP.keys()), "Gender"]
    if len(unknown) > 0:
        raise ValueError(
            f"Unrecognized Gender values: {unknown.unique().tolist()}. "
            "Support: 女/男, 0/1, Female/Male, F/M"
        )
    df = df.copy()
    df["Gender"] = df["Gender"].map(GENDER_MAP)
    return df


# --------------------------------------------------------------------------
# 风险分级
# --------------------------------------------------------------------------

def _classify_risk(age_acc: float) -> Tuple[str, str]:
    """根据年龄加速值进行风险分级。"""
    if pd.isna(age_acc):
        return "N/A", "Cannot evaluate without chronological age."

    if age_acc < 0:
        return (
            "lower risk",
            "Aging pace is slower than peers. "
            "Keep current healthy lifestyle.",
        )
    elif age_acc <= 5:
        return (
            "medium risk",
            "Aging pace slightly higher than peers. "
            "Regular checkups recommended.",
        )
    else:
        return (
            "higher risk",
            "Aging pace significantly higher than peers. "
            "Consult a doctor for targeted intervention.",
        )


# --------------------------------------------------------------------------
# 主预测函数
# --------------------------------------------------------------------------

def predict_bio_age(
    df: pd.DataFrame,
    fill_missing: bool = True,
) -> pd.DataFrame:
    """预测生物学年龄。

    完整流程：数据校验 → CpG 缺失填充 → 性别映射 →
    标准化 → 预测 → 风险分级。

    Args:
        df: 甲基化数据（含 SampleID、CpG 位点、Gender）。
        fill_missing: 缺失时是否用均值填充。

    Returns:
        包含预测结果的 DataFrame。
    """
    load_model()

    chr_pos_cols = _detect_cpg_columns(df)
    df = _fill_missing_cpgs(df, chr_pos_cols, fill_missing)

    if df["SampleID"].duplicated().any():
        raise ValueError("Duplicate SampleID found.")

    df = _validate_gender(df)

    input_cols = train_cpg_names + ["Gender"]
    data = df[input_cols].copy()
    data = data.fillna(train_features_mean)

    if data.isnull().any().any():
        null_cols = data.columns[data.isnull().any()].tolist()
        raise ValueError(f"Cannot fill missing values for: {null_cols}")

    data_scaled = scaler.transform(data.values)
    predicted = model.predict(data_scaled)
    chronological = df.get(
        "Chronological_Age", pd.Series([pd.NA] * len(df))
    )

    results = []
    for i in range(len(df)):
        age_acc = (
            predicted[i] - chronological.iloc[i]
            if pd.notna(chronological.iloc[i])
            else float("nan")
        )
        risk, suggestion = _classify_risk(age_acc)
        results.append({
            "Sample_ID": df["SampleID"].iloc[i],
            "Chronological_Age": (
                chronological.iloc[i]
                if pd.notna(chronological.iloc[i])
                else "N/A"
            ),
            "Predicted_Age": round(predicted[i], 1),
            "Age_Gap": round(age_acc, 1) if pd.notna(age_acc) else "N/A",
            "Risk_Level": risk,
            "Clinical_Suggestion": suggestion,
        })

    return pd.DataFrame(results)
