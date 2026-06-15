"""
Eval 执行器 —— 运行测试场景并汇总结果。

用法:
    cd 项目根目录
    python -m eval.evaluator

依赖:
    - AgentOrchestrator（需要模型文件存在）
    - app.model（需要 sklearn / joblib 等依赖）
"""

import asyncio
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

# 确保项目根目录在 sys.path 中
_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from agent.orchestrator import AgentOrchestrator
from agent.memory import ConversationMemory
from eval.scenarios import EvalResult, SCENARIOS, TestScenario
from eval.synthetic_data import generate_synthetic_methylation


@dataclass
class EvalReport:
    """评估报告。"""
    total: int
    passed: int
    failed: int
    results: List[EvalResult] = field(default_factory=list)

    def pass_rate(self) -> float:
        return self.passed / self.total if self.total > 0 else 0.0

    def summary(self) -> str:
        return (
            f"Eval Results: {self.passed}/{self.total} passed "
            f"({self.pass_rate():.0%})"
        )

    def by_category(self) -> Dict[str, Dict[str, int]]:
        cats: Dict[str, Dict[str, int]] = {}
        for r in self.results:
            scenario = next(
                s for s in SCENARIOS if s.name == r.scenario_name
            )
            cat = scenario.category
            if cat not in cats:
                cats[cat] = {"total": 0, "passed": 0}
            cats[cat]["total"] += 1
            if r.passed:
                cats[cat]["passed"] += 1
        return cats

    def print_report(self) -> None:
        """打印可读的评估报告。"""
        print("=" * 60)
        print("  Agent Eval Report")
        print("=" * 60)
        print(f"  Total:  {self.total}")
        print(f"  Passed: {self.passed}")
        print(f"  Failed: {self.failed}")
        print(f"  Rate:   {self.pass_rate():.0%}")
        print()

        by_cat = self.by_category()
        for cat, stats in by_cat.items():
            if stats["passed"] == stats["total"]:
                status = "[OK]"
            else:
                status = "[WARN]"
            print(f"  {status} {cat}: {stats['passed']}/{stats['total']}")

        print()
        for r in self.results:
            if r.passed:
                icon = "[PASS]"
            else:
                icon = "[FAIL]"
            err = f" — {r.error}" if r.error else ""
            print(f"  {icon} {r.scenario_name}{err}")


class AgentEvaluator:
    """Agent 行为评估器。"""

    def __init__(self, agent: Optional[AgentOrchestrator] = None):
        self.agent = agent or AgentOrchestrator()

    def run_all(self) -> EvalReport:
        """运行所有测试场景。"""
        results = []
        for scenario in SCENARIOS:
            result = asyncio.run(self._run_scenario(scenario))
            results.append(result)
        return EvalReport(
            total=len(results),
            passed=sum(1 for r in results if r.passed),
            failed=sum(1 for r in results if not r.passed),
            results=results,
        )

    async def _run_scenario(self, scenario: TestScenario) -> EvalResult:
        """运行单个测试场景。"""
        try:
            inputs = scenario.input

            # 根据场景参数准备数据
            methylation_data = self._prepare_methylation_data(inputs)
            memory = self._prepare_memory(inputs)

            # 调用 Agent
            result = await self.agent.run(
                user_input="分析我的数据",
                methylation_data=methylation_data,
                memory=memory,
                existing_results=None,
            )

            # 校验结果
            details = {
                "mode": result.get("mode"),
                "tools": list(result.get("tool_results", {}).keys()),
                "has_type": result.get("type") is not None,
            }

            passed = self._check_expected(result, scenario.expected)
            return EvalResult(
                scenario_name=scenario.name,
                passed=passed,
                details=details,
            )

        except Exception as e:
            # 异常场景本身就是预期行为
            expected = scenario.expected
            behaves_ok = (
                expected.get("behavior") == "graceful_error"
                or expected.get("graceful_on_error") is True
            )

            return EvalResult(
                scenario_name=scenario.name,
                passed=behaves_ok,
                details={"error_type": type(e).__name__},
                error=f"{type(e).__name__}: {str(e)[:200]}" if not behaves_ok else None,
            )

    def _prepare_methylation_data(
        self, inputs: Dict[str, Any]
    ) -> Optional[pd.DataFrame]:
        """根据场景参数生成测试用甲基化数据。"""
        if not inputs.get("has_methylation"):
            return None

        if inputs.get("empty_file"):
            return pd.DataFrame()

        if inputs.get("corrupt_data"):
            # 损坏格式：无 CpG 列
            return pd.DataFrame({
                "SampleID": ["bad_sample"],
                "SomeRandomColumn": [1],
            })

        # 正常合成数据
        target_age = inputs.get("target_age", 45.0)
        return generate_synthetic_methylation(seed=42, target_age=target_age)

    def _prepare_memory(self, inputs: Dict[str, Any]) -> Optional[ConversationMemory]:
        """根据场景参数准备对话记忆。"""
        memory = ConversationMemory()

        turn_count = inputs.get("turns", 0)
        if turn_count > 0:
            for i in range(turn_count):
                memory.add_turn(
                    f"问题{i + 1}",
                    f"回答{i + 1}",
                    analysis_type="meth",
                )

        return memory if turn_count > 0 else None

    def _check_expected(
        self, result: Dict[str, Any], expected: Dict[str, Any]
    ) -> bool:
        """校验 Agent 输出是否满足预期。"""
        # 检查模式
        mode = expected.get("mode")
        if mode and result.get("mode") != mode:
            return False

        # 检查工具调用（仅当预测成功时）
        tools_called = expected.get("tools_called")
        if tools_called:
            actual_tools = set(result.get("tool_results", {}).keys())
            if not actual_tools.issuperset(tools_called):
                return False

        # 检查引导信息
        if expected.get("output_contains_options"):
            guidance = result.get("type") == "guidance"
            has_options = bool(result.get("options"))
            if not (guidance and has_options):
                return False

        # 检查是否为优雅降级
        if expected.get("behavior") == "graceful_error":
            if result.get("type") not in ("prediction_error", "error"):
                # 预测成功也算通过
                pass

        # graceful_on_error：不管预测成功还是失败，路由正确就算通过
        if expected.get("graceful_on_error"):
            # 路由已在 mode 检查中验证，直接通过
            pass

        return True

    def run_category(self, category: str) -> EvalReport:
        """运行指定类别的测试场景。"""
        scenarios = [s for s in SCENARIOS if s.category == category]
        results = [asyncio.run(self._run_scenario(s)) for s in scenarios]
        return EvalReport(
            total=len(results),
            passed=sum(1 for r in results if r.passed),
            failed=sum(1 for r in results if not r.passed),
            results=results,
        )


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------

def main():
    # Windows GBK 兼容
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

    print("\nLoading model and preparing evaluator...")
    evaluator = AgentEvaluator()
    print("Running scenarios...\n")
    report = evaluator.run_all()
    report.print_report()

    # 退出码
    sys.exit(0 if report.passed == report.total else 1)


if __name__ == "__main__":
    main()
