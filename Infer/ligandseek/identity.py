"""统一配体身份字符串. canonical_smiles 返回三阶段共用的规范 SMILES, 此模块只在化学解释器内调用 RDKit."""

import json
import os
from pathlib import Path
import subprocess
import sys


def canonical_smiles(text: str) -> str | None:
    """从用户原文得到模型身份. RDKit 按固定原子顺序写出 SMILES, 省略四面体手性、双键构型和同位素标记.

    - text: 单种化学身份的用户输入, 如 C[C@H](O)C(=O)O.
    - 返回值: 无显式构型的字符串, 如 CC(O)C(=O)O. 不能解析或不含原子时返回 None.
    """
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(text)
    if molecule is None or not molecule.GetNumAtoms():
        return None
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    # 去同位素后显式氘可能成为可折叠的普通氢. 再读取一次使重复 Prepare、Infer 得到相同身份.
    if any(atom.GetAtomicNum() == 1 and atom.GetIsotope() for atom in molecule.GetAtoms()):
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(canonical), canonical=True, isomericSmiles=False)
    return canonical


# ================================================================================================
def normalize_identities(sample: dict, decisions: list[dict], previous: dict | None, truth: dict | None) -> list[dict]:
    """把本次读取的样本、手工选择、旧得分和可选归属统一到同一身份轴. 只修改内存, 原输入文件保留.

    - sample["identities"]: 原有 S 个身份字典, 每项的 smiles 是原文. 归并后按首次出现顺序保留 I 个身份, embedding、type_tag 和 ccd_id 分别取该身份首次给出的非空值.
    - decisions: 每项按 blob_id 指定手选身份. input_smiles 记录原文, smiles 变为模型身份, 解析失败时为 None.
    - previous: 可选的 find.json 或 match.json. 原身份轴 identity_smiles 与各区域的 identity_scores 一一对应. 旧轴中的同一规范身份取各有效分数的最大值, 映射到新轴时缺失填 None.
    - truth: 可选的窄范围归属, 仅规范 identity_types 的键和 occurrences 的 smiles, 每个真实拷贝继续独立保留.
    - 返回值: 长度 S 的列表, 保留原始 input_smiles 及其规范 smiles, 如 {"input_smiles":"C[C@H](F)Cl","smiles":"CC(F)Cl"}. 解析失败时 smiles 为 None.
    """
    # 同时规范样本身份、人工选择、历史得分轴与可选归属. 键只反映化学身份, 不合并真实拷贝.
    texts = [item["smiles"] for item in sample["identities"]]
    texts.extend(item["smiles"] for item in decisions if item.get("smiles") is not None)
    if previous is not None:
        texts.extend(previous["identity_smiles"])
    if truth:
        texts.extend(truth.get("identity_types", {}))
        texts.extend(item["smiles"] for item in truth["occurrences"] if item.get("smiles"))
    unique = list(dict.fromkeys(texts))
    keys = {}
    if unique:
        # 主推理环境可以没有 RDKit. 所有身份在一次子进程调用内用同一化学解释器规范化.
        python = os.environ.get("LIGANDSEEK_PREPARE_PYTHON", sys.executable)
        process = subprocess.run([python, "-B", Path(__file__)], input=json.dumps(unique), stdout=subprocess.PIPE, text=True, check=True)
        keys = dict(zip(unique, json.loads(process.stdout), strict=True))

    # 合并后身份轴按首次出现顺序. 优先取有向量的记录, 保留其余写法供 run.json 追溯.
    by_smiles = {}
    original_inputs = []
    for item in sample["identities"]:
        normalized = keys[item["smiles"]]
        original_inputs.append({"input_smiles": item["smiles"], "smiles": normalized})
        if normalized is None:
            continue
        if normalized not in by_smiles:
            by_smiles[normalized] = {**item, "smiles": normalized}
        else:
            for field in ("embedding", "type_tag", "ccd_id"):
                if not by_smiles[normalized].get(field) and item.get(field):
                    by_smiles[normalized][field] = item[field]
    sample["identities"] = list(by_smiles.values())
    for decision in decisions:
        if decision.get("smiles") is not None:
            decision["input_smiles"] = decision["smiles"]
            decision["smiles"] = keys[decision["smiles"]]
    if previous is not None:
        # 旧文件可以使用不同的候选顺序, 也可能给两个对映体分别打分. 先按原轴解读, 同一规范身份取已观察的最高分.
        source_axis = [keys[text] for text in previous["identity_smiles"]]
        for blob in previous["candidates"]:
            if blob["identity_scores"] is not None:
                scores = {}
                for key, score in zip(source_axis, blob["identity_scores"], strict=True):
                    if key and score is not None and (key not in scores or score > scores[key]):
                        scores[key] = score
                blob["identity_scores"] = [scores.get(key) for key in by_smiles]

    if truth:
        truth["identity_types"] = {keys[text]: tag for text, tag in truth.get("identity_types", {}).items() if keys[text]}
        for item in truth["occurrences"]:
            if item.get("smiles"):
                item["smiles"] = keys[item["smiles"]]
    return original_inputs


if __name__ == "__main__":
    # 标准输入与输出均为长度 I 的 JSON 列表. 每个原字符串在相同位置对应一个模型身份或 null.
    print(json.dumps([canonical_smiles(text) for text in json.load(sys.stdin)]))
