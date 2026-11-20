"""发布包的命令行入口. main 解析原始输入准备、端到端推理及分阶段审阅命令, 阶段交接与落盘见 pipeline.run_sample."""

import argparse
import os
from pathlib import Path

from .pipeline import run_sample


def main(argv: list[str] | None = None) -> None:
    """读取一个样本的命令行参数. argv 为参数字符串列表或 None, 后者表示读取当前进程参数.

    prepare 写出准备目录. run 将准备得到的 sample.json 交给 run_sample, 在同一次输出下写出 predictions/.
    infer、plan、build 可从已有样本或候选文件继续, 参数示例见 Infer/README.md.
    """
    parser = argparse.ArgumentParser(prog="ligandseek", description="配体位点检测、身份匹配与构象重建")
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (("run", "从原始输入准备并推理"), ("prepare", "从原始输入生成 sample.json"), ("infer", "从实验输入运行"), ("plan", "按已有候选制定对接计划"), ("build", "按已有候选生成构象")):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("sample", type=Path, help="本样本的输入说明 JSON")
        command.add_argument("--output", type=Path, required=True, help="本次运行的新目录")
        command.add_argument("--models", type=Path, default=os.getenv("LIGANDSEEK_MODELS"), required=name != "plan" and not os.getenv("LIGANDSEEK_MODELS"), help="模型资产目录, 可由安装环境 LIGANDSEEK_MODELS 指定")
        if name == "prepare":
            continue
        if name in ("run", "infer"):
            command.add_argument("--until", choices=("find", "match", "build"), default="build")
            command.add_argument("--find-seed", type=int, default=0)
            command.add_argument("--match-seed", type=int, default=0)
        else:
            command.add_argument("--candidates", type=Path, required=True, help="同样本的 find.json 或 match.json")
        command.add_argument("--k", type=int, help="本样本最多提交的 docking 请求数, 省略则不限")
        command.add_argument("--scope", choices=("wide", "narrow", "all"), default="wide")
        command.add_argument("--build-seed", type=int, default=0)
        command.add_argument("--foreground-threshold", type=float, help="可选前景概率阈值")
        command.add_argument("--decisions", type=Path, help="按 blob 指定操作、用于对接的 SMILES 或单请求种子")
        command.add_argument("--reference", type=Path, help="窄范围可用的实例归属与资格 JSON")
    args = parser.parse_args(argv)
    if args.command in ("run", "prepare"):
        from ligandseek_prepare import prepare_sample

        directory = args.output if args.command == "prepare" else args.output / "prepared"
        args.sample = prepare_sample(args.sample, args.models, directory)
        if args.command == "prepare":
            return
        args.command = "infer"
        args.output = args.output / "predictions"
    run_sample(**vars(args))


if __name__ == "__main__":
    main()
