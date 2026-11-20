"""在各阶段的解释器中执行正式推理函数. pipeline 通过 execute 传递本次样本, 子进程回传阶段结果.

模型目录下 libraries 保存当前科学库. Find 和 Match 使用 Pocket Plus 环境, Build 可使用独立的 PocketXMol 环境.
"""

import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


def execute(stage: str, sample: dict, candidates: list[dict], models: Path, directory: Path, seed: int) -> list[dict]:
    """在对应解释器调用 find、match 或 build. 每个子进程仅持有本阶段模型, 退出后释放显存.

    - stage: 阶段名, 如 "build". sample 含样本身份、实验图、受体及候选身份.
    - candidates: Find 使用空列表, Match 接收 Find 区域, Build 接收实际计划中的逐 blob 身份、种子和预测中心.
    - models: 固定模型资产根, 其中 libraries 下的科学源码随资产一同迁移.
    - directory: Find 和 Match 使用本次 find 目录, Build 使用 build 目录. seed 是 Find 或 Match 种子, Build 使用请求自身的 sampling_seed.
    - 返回值: 阶段函数返回的预测区域或实际提交记录, 逐字段定义见对应 find、match、build 包.
    """
    directory = directory.resolve()
    models = models.resolve()
    directory.mkdir(exist_ok=True, parents=True)
    libs = models / "libraries"
    pp, matcher, px = libs / "Pocket_Plus", libs / "Matcher", libs / "PocketXMol"
    ada = libs / "AdaLigand/Data_Preprocessing/Ori_Data"
    # PX 和 Pocket Plus 可各自含有同名的顶层模块. 当前引擎优先, 化学输入的子进程继承包搜索路径.
    roots = [px, matcher, ada] if stage == "build" else [pp, matcher, ada]
    env = os.environ.copy()
    search = list(map(str, roots))
    package_root = Path(__file__).resolve().parents[1]
    # 源码布局下只加入 Infer/. 安装 wheel 时, 各解释器从自身环境导入 ligandseek 及二进制依赖.
    if (package_root / "pyproject.toml").is_file() and package_root.name not in ("site-packages", "dist-packages"):
        search.insert(0, str(package_root))
    if env.get("PYTHONPATH"):
        search.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(search)
    python = env.get("LIGANDSEEK_BUILD_PYTHON", sys.executable) if stage == "build" else sys.executable
    # 两个 JSON 只负责跨解释器传递数据. 地图、置信度与姿态均由科学函数直接写到 directory.
    with TemporaryDirectory(dir=directory, prefix=".call-") as scratch:
        args = Path(scratch) / "arguments.json"
        result = Path(scratch) / "predictions.json"
        args.write_text(json.dumps([stage, sample, candidates, str(models), str(directory), seed]), encoding="utf-8")
        subprocess.run([python, "-B", "-m", "ligandseek.runtime", str(args), str(result)],
                       cwd=px if stage == "build" and px.is_dir() else directory, env=env, check=True)
        return json.loads(result.read_text(encoding="utf-8"))


if __name__ == "__main__":
    stage, sample, candidates, models, directory, seed = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    if stage == "find":
        from .find import find

        predictions = find(sample, Path(models), Path(directory), seed)
    elif stage == "match":
        from .match import match

        predictions = match(sample, candidates, Path(models), Path(directory), seed)
    else:
        from .build import build

        predictions = build(sample, candidates, Path(models), Path(directory))
    Path(sys.argv[2]).write_text(json.dumps(predictions, ensure_ascii=False, allow_nan=False), encoding="utf-8")
