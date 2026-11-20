"""用 SMI-TED Light 289M 编码配体. encode_records 在 output/smi_ted_289m 下按样本写 results.jsonl 与 candidate_<编号>.npz, 状态和向量按 candidate_id 对齐."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def encode_records(sample_id: str, records: list[dict], models: Path, output: Path) -> Path:
    """由单独的 SMI-TED 解释器编码本次样本, 返回样本目录的父目录 output/smi_ted_289m.

    - sample_id: 小写字符串, 标识记录所属的一个样本, 如 "demo".
    - records: 配体列表, 每项含 candidate_id (样本内唯一整数), object_key (化学身份来源), smiles (规范身份或 None), type_tag (已知类别或 None).
    - models: 发布模型根目录, libraries/AdaLigand 中含原编码器包装和向量写出函数.
    - output: 本次工作的目录. 子进程使用 LIGANDSEEK_SMITED_PYTHON 和 LIGANDSEEK_SMITED_MODEL 指定的解释器和模型.
    - smi_ted_289m/<sample_id>/results.jsonl: 每个 candidate_id 一条状态. status 为 encoded、no_smiles 或 model_failed, diagnostics 记录分词与有限性, output_file 为相对当前 JSONL 的 NPZ 文件名或 None.
    - smi_ted_289m/<sample_id>/candidate_<编号>.npz: embedding 是 (768,) float32 分子向量. 标量 candidate_id、pdb_id、object_key、prepared_smiles、model_smiles、model_name 描述身份和编码输入.
    """
    base = models.resolve() / "libraries/AdaLigand/Data_Preprocessing/Ori_Data_upgrade/ligand_language_model"
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(base), env.get("PYTHONPATH"))))
    python = env.get("LIGANDSEEK_SMITED_PYTHON", sys.executable)
    weights = env.get("LIGANDSEEK_SMITED_MODEL", str(models.resolve() / "smi_ted_289m"))
    command = [python, "-B", __file__, sample_id, weights, str(output.resolve())]
    output.mkdir(parents=True, exist_ok=True)
    log_path = output / "smi_ted.log"
    with log_path.open("w", encoding="utf-8") as log:
        print(json.dumps({"command": command}, ensure_ascii=False), file=log, flush=True)
        started = time.monotonic()
        completed = subprocess.run(command, input=json.dumps(records), text=True, env=env, stdout=log, stderr=subprocess.STDOUT)
        print(json.dumps({"elapsed_seconds": time.monotonic() - started, "returncode": completed.returncode}), file=log, flush=True)
    if completed.returncode:
        print(log_path.read_text(encoding="utf-8", errors="replace")[-8192:], file=sys.stderr)
        completed.check_returncode()
    return output / "smi_ted_289m"


# ================================================================================================
def encode_batch(sample_id: str, weight_dir: str, output: Path, records: list[dict]) -> None:
    """在 SMI-TED 解释器中调用原 run_smi_ted 并保存每种身份的状态.

    - sample_id: 一个样本的标识, 如 "demo".
    - weight_dir: SMI-TED Light 289M 官方代码、词表、检查点的本机目录.
    - output: 向量和状态写入 output/smi_ted_289m/<sample_id>/.
    - records: 当前输入身份或实例的列表, 字段与 encode_records 相同. 未知 type_tag 保留为 None, 科学计算只取 smiles.
    """
    from artifacts import ModelResult, PreparedLigand
    from encode_smi_ted import DEFAULT_BATCH_SIZE, MODEL_STAGE, run_smi_ted
    from model_results import ModelCandidate, ModelWork

    candidates, statuses = [], []
    for record in records:
        number, tag, smiles = record["candidate_id"], record["type_tag"], record["smiles"]
        if smiles is None:
            statuses.append(ModelResult(sample_id, number, tag, "no_smiles", None, None, None, None, {}))
            continue
        # 编码器读取本次字符串和编号. 候选身份的类别允许未知, 每个向量保留用户或结构模板给定的 object_key.
        prepared = PreparedLigand(sample_id, number, record["object_key"], "SMILES", tag, False, smiles, "smiles", {})
        candidates.append(ModelCandidate(sample_id, prepared))

    if candidates:
        source = Path(weight_dir) / "smi-ted/inference/smi_ted_light/load.py"
        spec = importlib.util.spec_from_file_location("ligandseek_official_smi_ted_light", source)
        official = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(official)
        model = official.load_smi_ted(folder=weight_dir, ckpt_filename="smi-ted-Light_40.pt", vocab_filename="bert_vocab_curated.txt")
        model.eval()
        # 原函数复用同一批的规范化、词表诊断、失败重试和 (B,768) 编码. B 是本次具有规范 SMILES 的配体数.
        work = ModelWork(candidates=candidates, results_by_pdb={sample_id: statuses},
                         expected_candidate_ids={sample_id: frozenset(r["candidate_id"] for r in records)},
                         active_pdb_ids=[sample_id], skipped_pdb_ids=[], source_failures=[])
        run_smi_ted(official, model, work, output, DEFAULT_BATCH_SIZE)

    directory = output / MODEL_STAGE / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    # 向量文件相对于本 JSONL 定位, 可把样本移到另一台机器.
    with (directory / "results.jsonl").open("w", encoding="utf-8") as stream:
        for status in sorted(statuses, key=lambda entry: entry.candidate_id):
            result = status.to_dict()
            if result["output_file"]:
                result["output_file"] = Path(result["output_file"]).name
            stream.write(json.dumps(result, ensure_ascii=False, allow_nan=False) + "\n")


if __name__ == "__main__":
    sample_id, weight_dir, output_dir = sys.argv[1:]
    encode_batch(sample_id, weight_dir, Path(output_dir), json.load(sys.stdin))
