"""为一个原始样本建立可迁移的 Infer 输入. prepare_sample 写 output/sample.json 及相对路径指向的受体、密度和语言向量, prepare.json 记录原文与来源."""

import json
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import mkdtemp


def prepare_sample(raw: Path, models: Path, output: Path) -> Path:
    """从 raw 指定的单样本 JSON 得到 output/sample.json. 一次准备全部成功后发布目录.

    - raw: UTF-8 JSON 路径. experimental_map 指原始实验图, smiles 是用户输入的字符串列表. receptor_cif 或按种类给出的 FASTA 决定受体来源.
    - models: 本次安装的模型资产根目录. libraries/AdaLigand 中保存所调用的密度与语言模型生产函数.
    - output: 尚不存在的新目录. 允许多个不同输出目录共享同一组只读输入.
    - sample.json: 单样本字典, sample_id 与 receptor_source 标识样本和受体条件. experimental_map、simulated_map、receptor_cif 指本目录中的两个 NPY 和一个 CIF. identities 按首次给出规范身份的顺序列出 smiles、embedding.
    - prepare.json: 记录原始图、文件路径、受体调用、分辨率来源和每个输入字符串对应的身份编号. 具体字段和失败示例见 INPUTS.md.
    - 返回值: 已发布的 sample.json 绝对路径, 可传给公开 Infer 入口.
    """
    raw, models, output = Path(raw).resolve(), Path(models).resolve(), Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    ada = models / "libraries/AdaLigand/Data_Preprocessing/Ori_Data"
    env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(ada), env.get("PYTHONPATH"))))
    python = env.get("LIGANDSEEK_PREPARE_PYTHON", sys.executable)
    # 失败时保留本次未完成目录及外部程序的真实输出. 并发调用拥有不同的工作目录, 只有完整样本会发布到 output.
    staging = Path(mkdtemp(dir=output.parent, prefix=f"{output.name}.unfinished-"))
    print(f"Prepare 工作目录 {staging}", flush=True)
    subprocess.run([python, "-B", __file__, str(raw), str(models), str(staging)], env=env, check=True)
    staging.rename(output)
    return output / "sample.json"


# ================================================================================================
def assemble_sample(raw_path: Path, models: Path, output: Path) -> None:
    """化学解释器的阶段入口. raw_path、models、output 均为绝对路径, output 是单次准备独占的目录.

    - raw_path: 用户样本 JSON. smiles 的下标是用户输入编号, 规范身份轴在本函数内按首次出现的顺序建立.
    - models: 已安装的模型资产根. 受体密度和 SMI-TED 使用各自的原生产入口.
    - output: 写 receptor.cif、exp/sim.npy 和几何 NPZ、smi_ted_289m 向量集、prepare.json 和 sample.json. 两个 JSON 的逐字段示例见 INPUTS.md.
    """
    from density import prepare_density
    from language import encode_records
    from receptor import prepare_receptor
    from ligandseek.identity import canonical_smiles
    from ligandseek.pipeline import read_json, write_json

    raw = read_json(raw_path)
    experiment = (raw_path.parent / raw["experimental_map"]).resolve()
    sample_id = str(raw.get("sample_id", experiment.name.split(".")[0])).lower()
    if not re.fullmatch(r"[a-z0-9][a-z0-9._-]*", sample_id):
        sys.exit("sample_id 需要以字母或数字开头, 并仅使用字母、数字、句点、下划线或短横线.")

    # 紧凑身份轴按首个有效写法出现的顺序建立. 无法解析的原文保留审计记录, 不占 Match 的模型身份槽.
    aliases, identity_index, records = [], {}, []
    original_smiles = [raw["smiles"]] if isinstance(raw["smiles"], str) else raw["smiles"]
    for original in original_smiles:
        smiles = canonical_smiles(original)
        if smiles and smiles not in identity_index:
            index = len(records)
            identity_index[smiles] = index
            records.append({"candidate_id": index, "object_key": smiles, "smiles": smiles, "type_tag": None})
        aliases.append({"input_smiles": original, "identity_index": identity_index.get(smiles)})

    sources = {"input": str(raw_path), "experimental_map": str(experiment), "models": str(models), "smiles_inputs": aliases,
               "smi_ted_model": os.environ.get("LIGANDSEEK_SMITED_MODEL", str(models / "smi_ted_289m"))}
    # 每个耗时步骤返回时即可在未完成目录查到它实际采用的参数. sample.json 最后写入.
    provenance = output / "prepare.json"
    write_json(provenance, sources)
    sources["receptor"] = receptor = prepare_receptor(raw, experiment, raw_path.parent, output)
    write_json(provenance, sources)
    sources["density"] = density = prepare_density(raw, experiment, output / "receptor.cif", raw_path.parent, output)
    write_json(provenance, sources)
    language = encode_records(sample_id, records, models, output) / sample_id
    statuses = [json.loads(line) for line in (language / "results.jsonl").read_text(encoding="utf-8").splitlines()]
    identities = []
    for status in statuses:
        # 保持向量与规范候选一一对齐. 官方编码失败或出现非有限分量时, 此槽保留名字, Match 返回空分数.
        embedding = None
        if status["status"] == "encoded" and status["diagnostics"]["all_finite"]:
            embedding = (language / status["output_file"]).relative_to(output).as_posix()
        identities.append({"smiles": status["prepared_smiles"], "embedding": embedding})
    sample = {"sample_id": sample_id, "experimental_map": "exp.npy", "simulated_map": "sim.npy", "receptor_cif": "receptor.cif",
              "receptor_source": receptor["source"], "map_convention": density["map_convention"], "identities": identities}
    write_json(output / "sample.json", sample)


if __name__ == "__main__":
    assemble_sample(*map(Path, sys.argv[1:]))
