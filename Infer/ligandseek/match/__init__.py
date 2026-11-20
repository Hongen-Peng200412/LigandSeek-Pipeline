"""以所选受体、Find 的局部区域和用户指定的身份轴运行 strongest-1.

match 返回逐区域的前景概率及身份得分. match_model.json 记录原验证集阈值及本样本用不到的身份槽. match/chemistry/ 保存本次生成的图及对称映射.
"""

import json
from pathlib import Path


def match(sample: dict, candidates: list[dict], models: Path, find_directory: Path, seed: int) -> list[dict]:
    """从当前模型库加载 strongest-1, 对 B 个区域逐个计算前景概率及用户 I 个精确 SMILES 的得分.

    - sample: sample_id 标记一种受体条件. identities 长度 I, 每项的 smiles 为精确字符串, embedding 指向 float32 (768,) NPY 或含相同 embedding 字段的 NPZ.
    - candidates: 长度 B, centered_box_index 定位 Find 的局部框. 输出保留原 blob_id、顺序和 Find 字段.
    - models/match/BEST.pt: 训练权重. 使用其中 model_config 重建模型, model_state 恢复权重, foreground_threshold 仅供独立评估.
    - find_directory: centered/F2_centered.npz 保存 Find 的局部特征, inputs/parse/<sample_id小写>/receptor_tokens.npz 提供同一受体的原子特征. seed 是本次 Match 的整数随机种子.
    - 返回每个区域的 identity_scores 为长度 I 的原始 logit 列表, 例如 [0.4, None, 1.7], 无法编码的身份保留空槽. foreground_probability 为二分类 logit 的 sigmoid.
    - find_directory.parent/match_model.json: JSON 字典记录 checkpoint 路径、foreground_threshold、按输入位置列出的有效身份下标 valid_slots 和 unsupported_identities 中的 smiles、error.
    """
    import numpy as np
    import torch
    from matcher.data.smiles_graph import SmilesGraphStore
    from matcher.datasets.candidate_data import _centered_fields, _load_npz_fields
    from matcher.model.config import MatcherModelConfig
    from matcher.model.ligand_matcher import ClosedSetMatcher
    from matcher.training.checkpointing import load_checkpoint

    from .batch import candidate_arrays, collate_inputs
    from .chemistry import prepare_graphs

    checkpoint_path = models / "match" / "BEST.pt"
    checkpoint = load_checkpoint(checkpoint_path, map_location="cpu")
    identities = sample["identities"]
    chemistry = find_directory.parent / "match" / "chemistry"
    excluded = prepare_graphs([item["smiles"] for item in identities], chemistry)
    store = SmilesGraphStore(chemistry / "smiles_graphs_v1.npz")

    # valid_slots 按用户顺序枚举. 即使两种纳入模型的身份之间夹有无效身份, 写回后仍能区分每一槽的含义.
    slots, graphs, language = [], [], []
    for slot, item in enumerate(identities):
        if item["smiles"] not in store:
            continue
        if not item.get("embedding"):
            excluded.append({"smiles": item["smiles"], "error": "missing SMI-TED embedding"})
            continue
        embedding_path = Path(item["embedding"])
        if embedding_path.suffix == ".npz":
            with np.load(embedding_path, allow_pickle=False) as archive:
                vector = archive["embedding"]
        else:
            vector = np.load(embedding_path, allow_pickle=False)
        if vector.shape != (768,) or not np.isfinite(vector).all():
            excluded.append({"smiles": item["smiles"], "error": "expected a finite 768-dimensional SMI-TED vector"})
            continue
        slots.append(slot)
        graphs.append(store.get(item["smiles"]))
        language.append(vector)
    record = {"checkpoint": str(checkpoint_path), "foreground_threshold": float(checkpoint["foreground_threshold"]),
              "valid_slots": slots, "unsupported_identities": excluded}
    (find_directory.parent / "match_model.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8",
    )
    if not candidates:
        return []

    torch.manual_seed(seed)
    config = MatcherModelConfig.from_dict(checkpoint["model_config"])
    model = ClosedSetMatcher(config)
    model.load_state_dict(checkpoint["model_state"], strict=True)
    del checkpoint
    model = model.to("cuda").eval()
    # (I_valid, 768), 语言向量第一维与 graphs 及 valid_slots 一一对应.
    vectors = np.asarray(language, dtype=np.float32).reshape(-1, 768)
    centered = _load_npz_fields(find_directory / "centered" / "F2_centered.npz", _centered_fields(config.variant))
    positions = {int(box): pos for pos, box in enumerate(centered["centered_box_index"])}
    tokens = find_directory / "inputs" / "parse" / sample["sample_id"].lower() / "receptor_tokens.npz"
    with np.load(tokens, allow_pickle=False) as atom:
        coordinates = np.asarray(atom["coords"], dtype=np.float32)
        # [N,49] + [N,1] -> [N,50], 主链标志追加到原 49 维受体化学特征末尾.
        features = np.concatenate((atom["feat"], atom["is_backbone"][:, None]), axis=1).astype(np.float32)

    result = []
    with torch.inference_mode():
        for start in range(0, len(candidates), 32):
            selected = candidates[start:start + 32]
            regions = [candidate_arrays(centered, positions[c["centered_box_index"]], coordinates, features)
                       for c in selected]
            batch = {k: v.to("cuda") for k, v in collate_inputs(regions, graphs, vectors).items()}
            del regions
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                output = model(batch) if slots else model.candidate_encoder(batch)
            foreground = torch.sigmoid(output["foreground_logit"]).float().cpu().tolist()
            logits = output["identity_logits"].float().cpu().tolist() if slots else [[] for _ in selected]
            for region, probability, values in zip(selected, foreground, logits):
                scores = [None] * len(identities)
                for slot, value in zip(slots, values):
                    scores[slot] = value
                result.append({**region, "identity_scores": scores, "foreground_probability": probability})
            del output, batch
    return result
