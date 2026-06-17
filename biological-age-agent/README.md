<p align="center">
  <img src="https://img.shields.io/badge/status-production-ready-brightgreen" alt="Status">
  <img src="https://img.shields.io/badge/python-3.12%2B-blue" alt="Python">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
  <img src="https://img.shields.io/badge/tests-19%2F19-passing-brightgreen" alt="Tests">
  <img src="https://img.shields.io/badge/eval-8%2F8-passing-brightgreen" alt="Eval">
</p>

<h1 align="center">🧬 Epigenetic Age — 生物学年龄预测 Agent</h1>

<p align="center">
  <strong>基于 DNA 甲基化数据的生物学年龄预测 + 智能分析 Agent</strong><br>
  从分子层面评估衰老速率，提供个性化干预建议
</p>

---

## 📋 目录

- [项目简介](#项目简介)
- [核心能力](#核心能力)
- [系统架构](#系统架构)
- [快速开始](#快速开始)
- [API 文档](#api-文档)
- [测试与评估](#测试与评估)
- [技术栈](#技术栈)
- [面试亮点](#面试亮点)
- [后续规划](#后续规划)

---

## 项目简介

**Epigenetic Age** 是一个端到端的生物学年龄预测系统，利用 DNA 甲基化数据（表观遗传时钟）评估个体的生物学年龄，并通过 Agent 系统生成个性化分析报告。

> 表观遗传时钟（Epigenetic Clock）是目前最准确的生物学年龄生物标志物之一。它基于基因组中特定 CpG 位点的甲基化水平，通过机器学习模型预测得出。生物学年龄与实际年龄的差值（Age gap）被广泛用于衡量衰老速率。

### 为什么是生物学年龄？

| 概念 | 说明 |
|------|------|
| **实际年龄** | 出生年份决定，不可改变 |
| **生物学年龄** | 由分子层面的衰老程度决定，**可被生活方式干预** |
| **年龄加速** | 生物学年龄 − 实际年龄。正值表示衰老快于同龄人 |

---

## 核心技术能力

### 🧪 甲基化年龄预测
- **预测模型** — 可替换成自定义的模型
- **年龄加速计算** — 生物学年龄与实足年龄的差值
- **三级风险分层** — 低风险 / 中风险 / 高风险，每级附带临床建议

### 🤖 AI Agent 智能分析
- **双模式调度** — 有数据时预测分析，无数据时知识问答
- **LLM 报告生成** — 基于 DeepSeek / OpenAI 兼容 API 生成个性化报告
- **分层记忆系统** — 短期记忆（最近 5 轮）+ 长时压缩记忆
- **知识库检索** — 内置 24 条衰老科学知识条目，覆盖机制 / 干预 / 风险因素
- **工具执行追踪** — 每一步工具调用可审计、可回溯

### 📊 工程质量保障
- **19 个单元测试** — 覆盖文件解析、知识检索、记忆系统
- **8 个 Eval 场景** — 覆盖正常流程 / 边界 / 异常恢复 / 记忆边界
- **用户反馈闭环** — 👍👎 评分机制 + 满意度统计
- **延迟监控** — 每次请求耗时追踪 + 聚合指标

---

## 系统架构

```
┌─────────────────────────────────────────────────────────┐
│                    用户界面 (Web)                       │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐    │
│  │  经典预测    │  │  Agent 对话  │  │  历史记录   │    │
│  └──────┬──────┘  └──────┬──────┘  └─────────────┘    │
└─────────┼─────────────────┼────────────────────────────┘
          │                 │
          ▼                 ▼
┌─────────────────────────────────────────────────────────┐
│                  FastAPI 服务层                         │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌────────┐ │
│  │ /predict │  │  /chat   │  │/feedback│  │/metrics│ │
│  └────┬─────┘  └────┬─────┘  └──────────┘  └────────┘ │
└───────┼──────────────┼──────────────────────────────────┘
        │              │
        ▼              ▼
┌─────────────────────────────────────────────────────────┐
│                  Agent Orchestrator                     │
│  ┌─────────────────────────────────────────────────┐   │
│  │  AgentContext → 模式判定 → 工具调度 → 报告生成  │   │
│  └──────────┬──────────────────────────┬───────────┘   │
│             │                          │               │
│       ┌─────▼─────┐            ┌───────▼───────┐      │
│       │  预测工具  │            │  LLM 报告     │      │
│       │ (本地模型) │            │ (DeepSeek等)  │      │
│       └─────┬─────┘            └───────┬───────┘      │
└─────────────┼──────────────────────────┼───────────────┘
              │                          │
              ▼                          ▼
┌─────────────────┐  ┌──────────────────────────────┐
│  自定义 Model   │  │     知识库 (24 entries)       │
│                 │  │  ┌─────┐ ┌──────┐ ┌──────┐  │
│                 │  │  │机制  │ │干预  │ │风险  │  │
└─────────────────┘  │  └─────┘ └──────┘ └──────┘  │
                     └──────────────────────────────┘
```

### 数据流

```
用户上传 .xlsx/.csv/.tsv ──→ 文件解析 ──→ CpG 位点检测
                                              │
                                              ▼
                                   缺失值填充（训练集均值）
                                              │
                                              ▼
                                   标准化 + model预测
                                              │
                                              ▼
                                   风险分级 + 知识检索
                                              │
                           ┌──────────────────┐
                           ▼                  ▼
                     LLM 报告生成        结构化模板输出
```

---

## 快速开始

### 前置要求

| 依赖 | 版本 |
|------|------|
| Python | 3.12+ |
| pip | 最新版 |
| Docker（可选）| 20.10+ |

### 本地运行

```bash
# 1. 克隆仓库
git clone <repo-url>
cd biological-age-agent

# 2. 安装依赖
pip install -r requirements.txt

# 3. 配置 LLM API（可选，不配也能运行预测）
#     创建 .env 文件：
echo 'LLM_API_KEY=sk-your-key-here' > .env
echo 'LLM_BASE_URL=https://api.deepseek.com/v1' >> .env
echo 'LLM_MODEL=deepseek-chat' >> .env

# 4. 启动服务
python run.py

# 5. 打开浏览器访问
open http://localhost:8000
```

### Docker 部署

```bash
docker compose up -d
# 访问 http://localhost:8000
```

---

## API 文档

| 方法 | 路径 | 说明 |
|------|------|------|
| **GET** | `/` | Web 界面 |
| **GET** | `/download-demo` | 下载示例数据文件 |
| **POST** | `/predict` | 甲基化数据预测（返回 Excel） |
| **POST** | `/chat` | Agent 智能对话分析 |
| **POST** | `/feedback` | 用户反馈评分 |
| **GET** | `/history/{session_id}` | 对话历史 |
| **GET** | `/metrics` | 服务监控指标 |

---

## 测试与评估

### 单元测试

```bash
python -m pytest tests/ -v
```

覆盖范围：
- ✅ 文件解析（xlsx / csv / tsv / txt）
- ✅ Demo 文件生成
- ✅ 知识库检索（关键词 / 风险等级）
- ✅ 对话记忆（压缩 / 反馈）

### Agent 行为评估

```bash
python -m eval.evaluator
```

8 个场景覆盖：

| 类别 | 场景 | 说明 |
|------|------|------|
| ✅ Happy Path | `methylation_prediction` | 正常预测流程 |
| ✅ Happy Path | `no_data_guidance` | 无数据时引导上传 |
| ✅ Boundary | `empty_methylation_file` | 空文件处理 |
| ✅ Boundary | `extreme_age_prediction` | 极端年龄健壮性 |
| ✅ Error | `malformed_data_format` | 格式异常恢复 |
| ✅ Error | `missing_model_file` | 模型缺失处理 |
| ✅ Memory | `memory_compression` | 长会话压缩 |
| ✅ Memory | `memory_across_sessions` | 跨会话记忆 |

---

## 技术栈

| 层 | 技术 |
|----|------|
| **后端框架** | FastAPI + Uvicorn |
| **机器学习** | model + scikit-learn |
| **Agent 框架** | Orchestrator + Tool 模式 |
| **LLM 集成** | OpenAI 兼容 API（DeepSeek / 任意） |
| **前端** | 原生 HTML + CSS + JavaScript |
| **数据** | pandas + numpy + openpyxl |
| **测试** | pytest |
| **部署** | Docker + docker-compose |
| **CI** | GitHub Actions |

---

## 面试亮点

这个项目在面试中可以从以下几个角度展开：

### 1. 技术选型能力
> "我选择了 Elastic Net做回归预测，适合高维数据和降维。模型用 joblib 序列化，API 用 FastAPI 部署——这是目前 ML 服务端最主流的组合。"

### 2. 工程化思维
> "整个系统分层清晰：数据层（模型预测）→ Agent 层（智能调度）→ API 层（路由分发）→ 展示层（Web UI）。每层职责单一，可独立测试。"

### 3. 质量意识
> "我写了 19 个单元测试覆盖核心函数，8 个 Eval 场景验证 Agent 行为，还接了 CI 自动运行。每次提交都保证测试全绿。"

### 4. 产品意识
> "报告底下加了👍👎反馈按钮，用户可以对结果评分。系统会追踪满意度趋势，低分输出会被自动标记——闭环优化。"

### 5. 持续学习
> "后续计划加入 FAHR-Face（哈佛最新 PyTorch 模型）做面部生物学年龄估算，与甲基化时钟做多模态交叉验证，进一步提升评估准确性。"

---

## 后续规划

- [ ] **FAHR-Face 集成** — 面部图像生物学年龄估算（哈佛 PyTorch 模型，权重待发布）
- [ ] **多模态融合** — 甲基化 + 面部交叉验证报告
- [ ] **用户系统** — 多用户支持、历史记录持久化
- [ ] **模型在线更新** — 支持增量训练，随着数据积累提升精度
- [ ] **LLM 输出漂移监控** — 统计学检测报告质量变化

---

<p align="center">
  Built with ❤️ for biological age research
</p>
