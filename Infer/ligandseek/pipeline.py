"""按单样本交接 Find、Match、对接计划与 Build 结果. run_sample 是发布包的 Python 入口.

run.json 记录解析后的输入与选项. find.json、match.json 和 plan.json 保存同一样本的全部可考察 blob.
build.json 保存实际提交的请求及姿态索引. 各文件通过 sample_id、blob_id 和用户的精确 SMILES 轴连接.
"""

import json
import os
from pathlib import Path
import subprocess
import sys

from .identity import normalize_identities
from .runtime import execute
from .selection import matched_smiles, plan_docking


def read_json(path: Path) -> dict | list:
    """读取 path 指定的 UTF-8 JSON. 返回单样本、候选或参考字典, 人工决定则返回按 blob 编号寻址的列表, 如 [{"blob_id": 17, "action": "skip"}]."""
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict) -> None:
    """把单个阶段的 payload 字典写入 path, 如 match.json 保存 sample_id、身份轴、seed 与逐 blob 预测. UTF-8 JSON 保留空值, 临时文件只在写入完成后改为正式文件名."""
    data = json.dumps(payload, ensure_ascii=False, allow_nan=False, indent=2)
    partial = path.with_suffix(".partial")
    partial.write_text(data + "\n", encoding="utf-8")
    partial.replace(path)


# ================================================================================================
def run_sample(
    sample: Path, output: Path, *, command: str = "infer", models: Path | None = None,
    candidates: Path | None = None, until: str = "build", scope: str = "wide",
    k: int | None = None, foreground_threshold: float | None = None, decisions: Path | None = None,
    reference: Path | None = None, find_seed: int = 0, match_seed: int = 0, build_seed: int = 0,
) -> None:
    """在 output 新目录保存一个样本的阶段文件. 成功写出所选阶段后返回 None.

    输入:
        - sample: Path, 单样本 JSON. sample_id 标识样本与受体条件, 如 "demo-ca2". 样本字段见 SELECTION.md.
            - experimental_map: 实验密度图路径, 体素按 ZYX 编号.
            - receptor_cif: 含所选受体的 mmCIF, 世界坐标单位 Å. prepare_inputs 提取受体原子, 也可读取含配体的完整复合物.
            - simulated_map: 相同受体在实验网格上的模拟密度路径.
            - identities: 候选身份列表. 同一次规范化得到相同 SMILES 的输入共享身份下标. 每种身份可带 type_tag 和 SMI-TED embedding 文件路径.
        - output: Path, 本次调用的新目录, 如 Path("reviewed").
        - command: infer 从实验图开始; plan 和 build 从 candidates 指定的逐区域预测文件继续.
        - models: Path 或 None, 冻结模型的本地资产根. 制定计划时可省略.
        - candidates: Path 或 None, 同一样本的 find.json 或 match.json. 读取旧得分时按文件内 identity_smiles 解释原下标, 再按本样本的身份顺序排列.
        - until: infer 的最后阶段, 取 find、match 或 build.
        - scope: 自动选择范围, 如 wide 按已知小分子类别筛选.
        - k: 每样本的非负整数请求上限, None 为不限.
        - foreground_threshold: 可选前景概率下限, None 关闭.
        - decisions: Path 或 None, 逐 blob 的人工选择, 例如为 17 指定闭集身份 CCN.
        - reference: Path 或 None, narrow 可选的参考身份类别、实例资格与归属.
        - find_seed: Find 种子, 如 7.
        - match_seed: Match 种子, 如 11.
        - build_seed: Build 基础种子, 如 23. decisions 可覆盖单个区域的种子.

    输出文件:
        - run.json: 字典, sample 保存解析过的样本输入, command 保存本次命令、六个路径和科学选项.
        - find.json: 字典, sample_id、identity_smiles、seed 和 candidates. 每个候选的定位字段见 stages.find.
        - match.json: 与 Find 共用样本和区域编号, candidates 增加前景概率、身份得分和 predicted_smiles.
        - plan.json: 样本标识、身份轴、scope、k、seed、foreground_threshold 和全部 candidates, 选择字段见 plan_docking.
        - build.json: 样本标识、身份轴及实际提交的 results, 单请求的姿态文件与排名见 stages.build.
    """
    sample_data = read_json(sample)
    # 实验图、受体和语言向量的相对路径均以 sample.json 所在目录为基准.
    for field in ("experimental_map", "receptor_cif", "simulated_map"):
        sample_data[field] = str((sample.parent / sample_data[field]).resolve())
    for identity in sample_data["identities"]:
        if identity.get("embedding") is not None:
            identity["embedding"] = str((sample.parent / identity["embedding"]).resolve())
    overrides = read_json(decisions) if decisions else []
    previous = read_json(candidates) if candidates is not None else None
    truth = read_json(reference) if scope == "narrow" and reference else None

    original_inputs = normalize_identities(sample_data, overrides, previous, truth)
    if models is not None and any(i.get("type_tag") is None for i in sample_data["identities"]):
        # 只传候选身份. 化学解释器从公共组分与 SMILES 得到同序类别, 识别不了的身份保持 None.
        env = os.environ.copy()
        ada = models.resolve() / "libraries/AdaLigand/Data_Preprocessing/Ori_Data"
        env["PYTHONPATH"] = os.pathsep.join([str(ada), *filter(None, [env.get("PYTHONPATH")])])
        python = env.get("LIGANDSEEK_PREPARE_PYTHON", sys.executable)
        types = subprocess.run([python, "-B", Path(__file__).with_name("chemical_types.py")], input=json.dumps(sample_data["identities"]),
                               text=True, stdout=subprocess.PIPE, env=env, check=True)
        for identity, tag in zip(sample_data["identities"], json.loads(types.stdout)):
            identity["type_tag"] = tag
    # (I,), 用户提供的身份轴. Find、Match、plan、Build 的同名字段保持这一顺序.
    smiles = [entry["smiles"] for entry in sample_data["identities"]]
    header = {"sample_id": sample_data["sample_id"], "identity_smiles": smiles}

    # 保存用户提供的样本来源与实际参数.
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = {"sample": sample, "output": output, "models": models, "candidates": candidates,
             "decisions": decisions, "reference": reference}
    options = {"command": command, "until": until, "k": k, "scope": scope, "find_seed": find_seed,
               "match_seed": match_seed, "build_seed": build_seed, "foreground_threshold": foreground_threshold}
    options.update({name: str(path.resolve()) if path is not None else None for name, path in paths.items()})
    write_json(output / "run.json", {"sample": sample_data, "command": options, "identity_inputs": original_inputs})

    # ———————— 自动端到端模式 ——————
    if command == "infer":
        blobs = execute("find", sample_data, [], models, output / "find", find_seed)
        write_json(output / "find.json", {**header, "seed": find_seed, "candidates": blobs})
        if until == "find":
            return
        # Match 保留完整候选顺序. 零候选时仍写出身份轴及模型的前景阈值, 以便评价该样本.
        blobs = execute("match", sample_data, blobs, models, output / "find", match_seed)
        for candidate in blobs:
            candidate["predicted_smiles"] = matched_smiles(candidate, smiles)
        write_json(output / "match.json", {**header, "seed": match_seed, "candidates": blobs})
        if until == "match":
            return
    else:
        blobs = previous["candidates"]
        # 身份规范化时已按原得分文件重排到当前轴. 这里仅接续属于本次样本的区域.
        for blob in blobs:
            blob.pop("reference", None)

    plan_identities = sample_data["identities"]
    # 窄范围允许显式传入某个沉积拷贝的资格. 无真值的正式批量推理使用 wide.
    if truth is not None:
        # 两张表用样本内的原编号连接: blob 17 -> occurrence 0 -> 这个拷贝的资格.
        occurrences = {r["occurrence_id"]: r for r in truth["occurrences"]}
        assignments = {r["blob_id"]: r["occurrence_id"] for r in truth["assignments"]}
        # 样本中缺少的身份类别可由 reference.identity_types 补足.
        reference_types = truth.get("identity_types", {})
        plan_identities = [
            {**identity, "type_tag": identity.get("type_tag") or reference_types.get(identity["smiles"])}
            for identity in plan_identities
        ]
        for candidate in blobs:
            occurrence = assignments.get(candidate["blob_id"])
            candidate["reference"] = occurrences[occurrence] if occurrence is not None else None

    # plan_docking 让人工决定优先, 并按样本预算给要提交的区域分配对接号.
    planned = plan_docking(blobs, plan_identities, k=k, scope=scope,
                           foreground_threshold=foreground_threshold, decisions=overrides, seed=build_seed)
    write_json(output / "plan.json", {**header, "scope": scope, "k": k, "seed": build_seed,
                                     "foreground_threshold": foreground_threshold, "candidates": planned})
    if command == "plan":
        return
    # 我们把中心、手选或预测身份及执行参数传给采样器; 参考坐标和真实归属留在 plan 中.
    request_fields = ("blob_id", "center_xyz_A", "chosen_smiles", "engine", "sampling_seed", "planned_docking_rank")
    # 长度 D 的请求列表, D 为实际安排的区域数. 每个中心按原实验图世界 XYZ 排列, 单位 Å.
    requests = [{field: item[field] for field in request_fields} for item in planned if item["plan_action"] == "dock"]
    results = execute("build", sample_data, requests, models, output / "build", build_seed) if requests else []
    write_json(output / "build.json", {**header, "results": results})
