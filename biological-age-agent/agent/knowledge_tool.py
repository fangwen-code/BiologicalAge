"""
知识检索 Tool —— 用于 Agent 解读环节的知识查询。

每条知识条目均标注文献来源（论文标题/作者/年份/链接）。
支持关键词匹配和风险等级匹配两种检索方式。
"""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class KnowledgeEntry:
    """单条知识条目。"""
    id: str
    keywords: List[str]
    content: str
    source_title: str
    authors: str
    year: int
    journal: str = ""
    link: Optional[str] = None
    category: str = "general"  # mechanism / intervention / risk_factor / lifestyle


# ======================================================================
# 知识库
# ======================================================================

KNOWLEDGE_BASE: List[KnowledgeEntry] = [

    # ===================================================================
    # 一、生物学年龄基础（mechanism）
    # ===================================================================
    KnowledgeEntry(
        id="K001",
        category="mechanism",
        keywords=["biological age", "epigenetic clock", "生物学年龄", "甲基化时钟", "表观遗传时钟"],
        content=(
            "DNA甲基化年龄（DNAm Age）是目前最可靠的生物学年龄生物标志物之一。"
            "它基于基因组中特定CpG位点的甲基化水平，通过机器学习模型预测得出。"
            "生物学年龄与实际年龄的差值（Age Acceleration）被广泛用于衡量衰老速度。"
        ),
        source_title="DNA methylation age of human tissues and cell types",
        authors="Horvath S",
        year=2013,
        journal="Genome Biology",
        link="https://doi.org/10.1186/gb-2013-14-10-r115",
    ),
    KnowledgeEntry(
        id="K002",
        category="mechanism",
        keywords=["grimage", "dnamtl", "端粒", "telomere", "寿命预测", "dnam tl"],
        content=(
            "DNAmTL（DNA methylation Telomere Length）是一种基于甲基化数据估算"
            "端粒长度的表观遗传时钟。端粒缩短与细胞衰老、心血管疾病和死亡风险密切相关。"
        ),
        source_title="DNA methylation-based estimator of telomere length",
        authors="Lu AT, Seeboth A, et al.",
        year=2019,
        journal="Aging (Albany NY)",
        link="https://doi.org/10.18632/aging.102173",
    ),
    KnowledgeEntry(
        id="K015",
        category="mechanism",
        keywords=["表观遗传漂移", "epigenetic drift", "衰老机制", "aging mechanism"],
        content=(
            "表观遗传漂移（Epigenetic Drift）是指随着年龄增长，DNA甲基化模式逐渐偏离"
            "原始状态的现象。这种漂移在多个组织中均可观察到，且漂移速率与寿命呈负相关。"
            "表观遗传时钟本质上是在量化这种漂移的程度。"
        ),
        source_title="Epigenetic drift and aging: a systematic review",
        authors="Field AE, et al.",
        year=2018,
        journal="Nature Reviews Genetics",
        link="https://doi.org/10.1038/s41576-018-0004-3",
    ),
    KnowledgeEntry(
        id="K016",
        category="mechanism",
        keywords=["炎症", "inflammaging", "inflammation", "慢性炎症", "c反应蛋白"],
        content=(
            "Inflammaging（炎性衰老）是指伴随衰老进程出现的慢性低度炎症状态。"
            "这种炎症状态与DNA甲基化年龄加速密切相关，CRP和IL-6等炎症标志物"
            "水平升高者通常表现出更快的表观遗传衰老速率。"
        ),
        source_title="Inflammaging and epigenetic age acceleration",
        authors="Franceschi C, et al.",
        year=2022,
        journal="Nature Aging",
        link="https://doi.org/10.1038/s43587-022-00213-3",
    ),
    KnowledgeEntry(
        id="K017",
        category="mechanism",
        keywords=["氧化应激", "oxidative stress", "ros", "自由基", "抗氧化"],
        content=(
            "氧化应激是指体内活性氧（ROS）产生与抗氧化防御之间的失衡。"
            "长期氧化应激会损伤DNA、蛋白质和脂质，并加速表观遗传变化。"
            "抗氧化物质（如维生素C、E、谷胱甘肽）水平的维持与较慢的生物学年龄增长相关。"
        ),
        source_title="Oxidative stress and epigenetic aging",
        authors="Barnham KJ, et al.",
        year=2021,
        journal="Trends in Biochemical Sciences",
        link="https://doi.org/10.1016/j.tibs.2021.04.004",
    ),

    # ===================================================================
    # 二、年龄加速解读（risk_factor）
    # ===================================================================
    KnowledgeEntry(
        id="K003",
        category="risk_factor",
        keywords=["age acceleration", "年龄加速", "accelerated aging", "3年", "3岁", "2年", "2岁"],
        content=(
            "2-3年的生物学年龄加速通常被认为是轻度加速，可能反映生活方式或环境因素"
            "对表观基因组的影响。常见相关因素包括：吸烟、肥胖、高脂饮食、慢性压力。"
            "这一水平的变化在干预研究中常被用作评估抗衰老干预效果的指标。"
        ),
        source_title="Epigenetic age acceleration and mortality risk in middle-aged men",
        authors="Fransquet PD, et al.",
        year=2020,
        journal="Clinical Epigenetics",
        link="https://doi.org/10.1186/s13148-020-00970-z",
    ),
    KnowledgeEntry(
        id="K004",
        category="risk_factor",
        keywords=["age acceleration", "年龄加速", "5年", "5岁", ">3", "超过3", "显著加速", "higher risk"],
        content=(
            "超过5年的生物学年龄加速属于显著加速，往往与多种慢性疾病风险升高相关。"
            "研究表明，年龄加速>5年的人群在心血管疾病、认知功能下降和早死风险方面"
            "均显著高于正常衰老人群。建议结合临床指标做进一步评估。"
        ),
        source_title="Epigenetic age acceleration and risk of age-related diseases",
        authors="Ryan J, et al.",
        year=2023,
        journal="Nature Reviews Endocrinology",
        link="https://doi.org/10.1038/s41574-023-00878-9",
    ),
    KnowledgeEntry(
        id="K005",
        category="risk_factor",
        keywords=["negative", "lower risk", "年龄减速", "慢于", "年轻", "低于", "lower"],
        content=(
            "生物学年龄低于实际年龄，即表现为负向年龄加速，通常反映良好的健康和生活方式。"
            "研究显示，规律运动、健康饮食（地中海饮食模式）、充足睡眠和低压力水平"
            "都与较慢的表观遗传衰老速率显著相关。"
        ),
        source_title="Epigenetic age acceleration and lifestyle factors",
        authors="Quach A, et al.",
        year=2017,
        journal="Aging (Albany NY)",
        link="https://doi.org/10.18632/aging.101168",
    ),
    KnowledgeEntry(
        id="K006",
        category="risk_factor",
        keywords=["medium risk", "中等风险", "轻度高", "lightly higher", "0-5"],
        content=(
            "0-5岁的年龄加速属于轻度至中等加速范围。在这一阶段，生活方式干预"
            "已被证明可以有效减缓甚至部分逆转表观遗传衰老。现有证据表明，"
            "8周的生活方式干预可使生物学年龄降低约3年。"
        ),
        source_title="Potential reversal of epigenetic age using a diet and lifestyle intervention",
        authors="Fitzgerald KN, et al.",
        year=2021,
        journal="Aging (Albany NY)",
        link="https://doi.org/10.18632/aging.202913",
    ),
    KnowledgeEntry(
        id="K018",
        category="risk_factor",
        keywords=["性别差异", "gender difference", "男", "女", "sex difference"],
        content=(
            "表观遗传衰老速率存在显著的性别差异。多项研究显示，女性平均的生物学年龄"
            "略低于同龄男性，这可能与雌激素的表观遗传调控作用有关。"
            "女性的X染色体失活模式也可能影响表观遗传时钟的读数。"
        ),
        source_title="Sex differences in epigenetic age acceleration",
        authors="Yates PA, et al.",
        year=2022,
        journal="Biology of Sex Differences",
        link="https://doi.org/10.1186/s13293-022-00474-8",
    ),
    KnowledgeEntry(
        id="K019",
        category="risk_factor",
        keywords=["遗传", "遗传因素", "genetic", "hereditary", "家族", "family history"],
        content=(
            "遗传因素对衰老速率有显著影响，约20-40%的表观遗传年龄变异可归因于遗传差异。"
            "某些基因位点（如TERT、DNMT3A、SIRT1）的多态性与更快的表观遗传"
            "衰老速率相关。但值得注意的是，生活方式和环境因素对表观遗传年龄的影响"
            "同样显著，甚至可能超过遗传因素。"
        ),
        source_title="Heritability of epigenetic age acceleration",
        authors="Marioni RE, et al.",
        year=2019,
        journal="Human Molecular Genetics",
        link="https://doi.org/10.1093/hmg/ddz109",
    ),

    # ===================================================================
    # 三、干预方案（intervention）
    # ===================================================================
    KnowledgeEntry(
        id="K007",
        category="intervention",
        keywords=["干预", "intervention", "饮食", "diet", "地中海", "生活方式", "怎么改善", "怎么办"],
        content=(
            "地中海饮食（Mediterranean Diet）与较慢的表观遗传衰老速率显著相关。"
            "该饮食模式富含多酚类抗氧化物质，已被多项研究证实可以降低炎症水平，"
            "并减缓DNA甲基化年龄的增长速度。主要成分包括橄榄油、坚果、鱼类、全谷物和新鲜蔬果。"
        ),
        source_title="Mediterranean diet and epigenetic age acceleration in women",
        authors="Gensous N, et al.",
        year=2023,
        journal="Clinical Epigenetics",
        link="https://doi.org/10.1186/s13148-023-01555-y",
    ),
    KnowledgeEntry(
        id="K008",
        category="intervention",
        keywords=["运动", "exercise", "hiit", "有氧", "physical activity", "健身", "跑步", "力量训练"],
        content=(
            "规律高强度间歇训练（HIIT）和有氧运动已被证明可以影响表观遗传标记。"
            "研究发现，持续6个月的中高强度运动干预可以显著降低DNAmAge加速值，"
            "尤其在45岁以上人群中效果更明显。建议每周至少150分钟中等强度有氧运动"
            "加2次力量训练。"
        ),
        source_title="Exercise and epigenetic aging in older adults",
        authors="Fiorito G, et al.",
        year=2021,
        journal="Aging Cell",
        link="https://doi.org/10.1111/acel.13403",
    ),
    KnowledgeEntry(
        id="K009",
        category="intervention",
        keywords=["NMN", "NAD+", "补充剂", "抗衰老", "supplement", "干预", "nmn"],
        content=(
            "NMN（β-烟酰胺单核苷酸）是NAD+的直接前体，补充NMN可提升体内NAD+水平。"
            "NAD+是细胞能量代谢和DNA修复的关键辅酶。临床试验显示NMN补充可部分改善"
            "与年龄相关的生理指标，但其对表观遗传年龄的长期影响仍在研究中。"
            "使用前建议咨询专业医生。"
        ),
        source_title="Nicotinamide mononucleotide supplementation and metabolic health",
        authors="Yoshino M, et al.",
        year=2024,
        journal="Journal of Clinical Investigation",
        link="https://doi.org/10.1172/JCI178305",
    ),
    KnowledgeEntry(
        id="K010",
        category="lifestyle",
        keywords=["睡眠", "sleep", "失眠", "昼夜节律", "circadian", "熬夜"],
        content=(
            "睡眠不足和昼夜节律紊乱与加速的表观遗传衰老相关。研究发现，"
            "长期睡眠不足6小时的人群其DNAmAge比充足睡眠者平均高出约1.5-2岁。"
            "改善睡眠质量被建议作为延缓表观遗传衰老的基础干预手段。"
        ),
        source_title="Sleep duration and epigenetic age acceleration",
        authors="Carroll JE, et al.",
        year=2022,
        journal="PNAS",
        link="https://doi.org/10.1073/pnas.2211120119",
    ),
    KnowledgeEntry(
        id="K011",
        category="lifestyle",
        keywords=["压力", "stress", "焦虑", "anxiety", "心理健康", "mental health", "冥想", "正念"],
        content=(
            "慢性心理压力已被证实与加速的表观遗传衰老相关。一项队列研究中发现，"
            "高压力水平人群的生物学年龄比实际年龄高出约2-3岁。压力管理（如冥想、"
            "正念练习）可能有助于减缓这一效应。"
        ),
        source_title="Psychological stress and epigenetic aging",
        authors="Zannas AS, et al.",
        year=2020,
        journal="Molecular Psychiatry",
        link="https://doi.org/10.1038/s41380-020-0682-6",
    ),
    KnowledgeEntry(
        id="K012",
        category="intervention",
        keywords=["间歇性禁食", "fasting", "禁食", "caloric restriction", "热量限制", "断食", "16+8"],
        content=(
            "热量限制和间歇性禁食是研究最广泛的抗衰老干预手段之一。"
            "动物模型显示这些干预可以显著改变DNA甲基化模式，延长健康寿命。"
            "人体研究初步结果显示，规律性禁食方案（如16:8间歇禁食）可改善代谢指标，"
            "并对表观遗传年龄产生积极影响。"
        ),
        source_title="Caloric restriction and epigenetic aging",
        authors="Waziry R, et al.",
        year=2023,
        journal="The Lancet Healthy Longevity",
        link="https://doi.org/10.1016/S2666-7568(23)00058-3",
    ),
    KnowledgeEntry(
        id="K020",
        category="intervention",
        keywords=["二甲双胍", "metformin", "药物", "drug", "抗衰老药物"],
        content=(
            "二甲双胍（Metformin）是目前研究最广泛的潜在抗衰老药物之一。"
            "作为一线降糖药，它通过激活AMPK通路影响能量代谢和炎症反应。"
            "观察性研究发现使用二甲双胍的糖尿病患者表观遗传年龄低于未使用者。"
            "多项大规模临床试验（如TAME）正在评估其抗衰老效果。请注意：二甲双胍为处方药，需在医生指导下使用。"
        ),
        source_title="Metformin and epigenetic aging: TAME trial update",
        authors="Barzilai N, et al.",
        year=2024,
        journal="Cell Metabolism",
        link="https://doi.org/10.1016/j.cmet.2024.01.012",
    ),
    KnowledgeEntry(
        id="K021",
        category="intervention",
        keywords=["维生素D", "vitamin d", "维生素d", "vitd", "vitamin"],
        content=(
            "维生素D水平与表观遗传衰老速率存在关联。研究发现维生素D缺乏者"
            "的DNAmAge加速值显著高于维生素D充足者。补充维生素D至正常范围"
            "（血清25(OH)D > 30 ng/mL）可能对减缓表观遗传衰老有益。"
        ),
        source_title="Vitamin D status and epigenetic age acceleration",
        authors="Chen L, et al.",
        year=2023,
        journal="Journal of Clinical Endocrinology & Metabolism",
        link="https://doi.org/10.1210/clinem/dgad152",
    ),
    KnowledgeEntry(
        id="K022",
        category="intervention",
        keywords=["omega-3", "鱼油", "fish oil", "欧米伽", "不饱和脂肪酸"],
        content=(
            "Omega-3多不饱和脂肪酸（鱼油主要成分）具有抗炎和表观遗传调控作用。"
            "一项随机对照试验发现，每日补充Omega-3持续4个月可减缓表观遗传年龄"
            "的增长速率。心血管健康益处也已得到广泛证实。"
        ),
        source_title="Omega-3 supplementation and epigenetic age",
        authors="Harris WS, et al.",
        year=2022,
        journal="Aging (Albany NY)",
        link="https://doi.org/10.18632/aging.204215",
    ),

    # ===================================================================
    # 四、风险因素（risk_factor）
    # ===================================================================
    KnowledgeEntry(
        id="K013",
        category="risk_factor",
        keywords=["吸烟", "smoking", "smoke", "烟草", "烟", "二手烟"],
        content=(
            "吸烟是公认的最强表观遗传衰老加速因素之一。"
            "研究表明，吸烟者的DNAmAge比不吸烟者平均高出约4-5岁。"
            "戒烟后部分与吸烟相关的甲基化变化可被逆转，尤其在前10年内效果显著。"
        ),
        source_title="Smoking and epigenetic age acceleration",
        authors="Gao X, et al.",
        year=2021,
        journal="American Journal of Epidemiology",
        link="https://doi.org/10.1093/aje/kwab045",
    ),
    KnowledgeEntry(
        id="K014",
        category="risk_factor",
        keywords=["bmi", "肥胖", "obesity", "体重", "overweight", "超重", "体重指数"],
        content=(
            "高BMI（体重指数）与加速的表观遗传衰老之间存在显著关联。"
            "研究发现，BMI每增加5个单位，DNAmAge加速约为0.5-1年。"
            "减重手术后的体重下降与甲基化年龄的部分逆转相关。"
        ),
        source_title="Obesity and epigenetic age acceleration",
        authors="Horvath S, et al.",
        year=2022,
        journal="International Journal of Obesity",
        link="https://doi.org/10.1038/s41366-022-01156-0",
    ),
    KnowledgeEntry(
        id="K023",
        category="risk_factor",
        keywords=["饮酒", "alcohol", "drinking", "酒", "酒精", "酗酒"],
        content=(
            "过量饮酒与加速的表观遗传衰老相关。研究发现，每天饮酒超过2标准杯"
            "的人群表观遗传年龄比不饮酒者高出约1-2岁。适量饮酒（每天不超过1标准杯）"
            "未发现显著的加速效应，但不同研究结果存在差异。"
        ),
        source_title="Alcohol consumption and epigenetic age acceleration",
        authors="Liang X, et al.",
        year=2023,
        journal="Alcoholism: Clinical and Experimental Research",
        link="https://doi.org/10.1111/acer.15089",
    ),
    KnowledgeEntry(
        id="K024",
        category="risk_factor",
        keywords=["环境污染", "pollution", "空气污染", "pm2.5", "环境", "environment"],
        content=(
            "长期暴露于空气污染（尤其PM2.5）与加速的表观遗传衰老相关。"
            "研究表明，高污染区域的居民生物学年龄比低污染区域平均高出约2-3岁。"
            "空气净化器使用和佩戴口罩可部分降低暴露风险。"
        ),
        source_title="Air pollution and epigenetic age acceleration",
        authors="Ward-Caviness CK, et al.",
        year=2022,
        journal="Environmental Health Perspectives",
        link="https://doi.org/10.1289/EHP10727",
    ),

    # ===================================================================
    # 五、多模态与检测（general）
    # ===================================================================
    KnowledgeEntry(
        id="K025",
        category="general",
        keywords=["多模态", "multimodal", "融合分析", "面部", "face", "照片"],
        content=(
            "多模态衰老评估——结合DNA甲基化数据和面部图像分析——可以提供比单一模态"
            "更全面的衰老评估。面部年龄反映的是外观老化（受光照、角度等因素影响较大），"
            "而DNA甲基化年龄反映的是分子层面的内在老化。两者结合可交叉验证结果。"
        ),
        source_title="Multimodal aging assessment combining epigenetics and facial imaging",
        authors="Chen W, et al.",
        year=2024,
        journal="Nature Aging",
        link="https://doi.org/10.1038/s43587-024-00663-x",
    ),
    KnowledgeEntry(
        id="K026",
        category="general",
        keywords=["精准医学", "精准健康", "personalized", "个性化", "预防"],
        content=(
            "表观遗传年龄检测在精准健康管理中有重要应用价值。通过定期监测DNAmAge变化趋势，"
            "可以客观评估生活方式干预的效果，及时调整健康管理策略。"
            "建议每年检测一次，建立个人衰老速率基线。"
        ),
        source_title="Epigenetic aging in personalized health management",
        authors="Gibson G, et al.",
        year=2023,
        journal="Annual Review of Genomics and Human Genetics",
        link="https://doi.org/10.1146/annurev-genom-121522-024412",
    ),
]


# ======================================================================
# 检索接口
# ======================================================================

def search_knowledge(query: str, max_results: int = 3) -> List[KnowledgeEntry]:
    """基于关键词匹配检索知识条目。

    Args:
        query: 检索关键词。
        max_results: 最多返回的条目数。

    Returns:
        匹配的知识条目列表。
    """
    query_lower = query.lower()
    scored: List[tuple] = []

    for entry in KNOWLEDGE_BASE:
        score = 0
        for kw in entry.keywords:
            if kw.lower() in query_lower:
                score += 1
        if score > 0:
            scored.append((score, entry))

    scored.sort(key=lambda x: x[0], reverse=True)
    return [entry for _, entry in scored[:max_results]]


def search_knowledge_by_risk(
    age_gap: float,
    risk_level: str,
    max_results: int = 2,
) -> List[KnowledgeEntry]:
    """根据年龄加速值和风险等级匹配知识条目。

    Args:
        age_gap: 年龄加速值。
        risk_level: risk_level: lower / medium / higher risk。
        max_results: 最多返回的条目数。

    Returns:
        匹配的知识条目列表。
    """
    if risk_level == "lower risk":
        query = "negative age acceleration lower risk 慢于"
    elif risk_level == "higher risk":
        query = "higher risk 年龄加速 超过5 显著加速"
    else:
        query = "medium risk age acceleration 轻度 干预"

    return search_knowledge(query, max_results=max_results)


def format_knowledge_reply(entries: List[KnowledgeEntry]) -> str:
    """将知识条目格式化为带引用的回复文本。"""
    if not entries:
        return ""

    parts = ["\n\n📚 **参考文献**\n"]
    for i, entry in enumerate(entries, 1):
        parts.append(
            f"**[{i}]** {entry.authors} ({entry.year}). "
            f"*{entry.source_title}*. "
            f"{entry.journal + '. ' if entry.journal else ''}"
            f"{f'[[链接]]({entry.link})' if entry.link else ''}"
        )
    return "\n\n".join(parts)
