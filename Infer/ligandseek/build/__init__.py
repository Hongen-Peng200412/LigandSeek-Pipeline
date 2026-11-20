"""在预测位点调用 DG-PocketXMol Centre 并按模型原规则对姿态排名.

build 按所选请求顺序发布 directory/<blob_id>/. result.json、candidates.json、poses.sdf 和 confidence.npz 沿用原采样器格式. ranking.json 保存逐姿态评分和模型首选编号.
"""

import json
from pathlib import Path
from tempfile import TemporaryDirectory


def build(sample: dict, requests: list[dict], models: Path, directory: Path) -> list[dict]:
    """给每个选定区域的精确 SMILES 生成世界坐标构象. 一次提交包含 50 个姿态, 每个姿态采样 100 步.

    - sample: sample_id 对应本次受体与实验图. 已有候选独立进入 Build 时, 用 experimental_map、receptor_cif、simulated_map 生成同网格输入.
    - requests: 按计划号排列的区域列表. blob_id 标识样本内区域, center_xyz_A 是 Find 的世界 XYZ 中心 (3,), 单位 Å, chosen_smiles 是自动或人工指定的身份, sampling_seed 是该请求独立的整数种子.
    - models/build/train.yml、step=24800.ckpt: 原 Centre 配置及训练参数. 采样使用训练模型的 RA 全原子受体、中心口袋、32 近邻与局部密度条件.
    - directory/<blob_id>/result.json: 采样器原记录, 包含种子、成功数、采样步数、峰值显存及输出文件名. scientific_input.json 保存实际中心、身份与种子. candidates.json 按原 sample_index 保留每个尝试, 包括失败.
    - poses.sdf: 成功生成的重原子世界 XYZ, Å, 由零基 sdf_index 索引. confidence.npz 的 sample_index 连接原采样编号, confidence_pos_traj 形状 (成功姿态数, 重原子数, 步数).
    - ranking.json: candidate_metrics 沿原采样编号, 每项含 cfd_traj、碰撞、立体、自评分. 按 self_ranking 降序及 sample_index 升序决定 top1_sample_index, 无可评分姿态时该编号为 null.
    - 返回按提交次序排列的字典列表. blob_id 是原区域编号, docking_rank 从 1 起, status 为 success、partial 或 failed. pose_file 是相对 directory 的路径或 null, failure_reason 保存第一个失败原因或 null.
    """
    from easydict import EasyDict
    from rdkit import Chem

    from docking.assets import read_receptor
    from docking.dataset import ConditionedDockingDataset
    from docking.evaluation import build_receptor_molecule, rank_candidate_metrics, score_saved_candidates
    from docking.sampling import build_sampling_noiser, load_sampling_runtime, sample_occurrence
    from docking.smiles import read_smiles_graph
    from ligandseek.match.chemistry import prepare_graphs

    if not requests:
        return []
    sample_id = sample["sample_id"].lower()
    inputs = directory.parent / "find" / "inputs"
    if not inputs.is_dir():
        from ligandseek.find.inputs import prepare_inputs

        inputs = directory / "inputs"
        prepare_inputs(sample, inputs, models / "ccd")

    # 自动链复用 Match 本次生成的图. 独立 Build 或手工指定集合外身份时, 为所有实际请求重新构图.
    directory.mkdir(parents=True, exist_ok=True)
    chemistry = directory.parent / "match" / "chemistry"
    supplied = {identity["smiles"] for identity in sample["identities"]}
    if not chemistry.is_dir() or any(r["chosen_smiles"] not in supplied for r in requests):
        chemistry = directory / "chemistry"
        prepare_graphs([r["chosen_smiles"] for r in requests], chemistry)
    config = EasyDict(
        model_name="local_cov-C", train_config=str(models / "build" / "train.yml"),
        checkpoint=str(models / "build" / "step=24800.ckpt"), receptor_branch="RA",
        device="cuda", num_candidates=50, batch_size=50, num_steps=100, split="inference",
        dataset=EasyDict(root=str(inputs), smiles_root=str(chemistry), knn=32, pocket_mode="center"),
    )
    training, model_config, model, featurizer, transforms, sampling = load_sampling_runtime(config)
    noiser = build_sampling_noiser(config, training, sampling, featurizer)
    records = [{
        "pdb_id": sample_id, "candidate_id": request["blob_id"], "prepared_smiles": request["chosen_smiles"],
        "given_center_xyz_A": request["center_xyz_A"], "sampling_seed": request["sampling_seed"],
        "views": [], "center_offset_xyz_A": [0.0, 0.0, 0.0],
    } for request in requests]
    dataset = ConditionedDockingDataset(config.dataset, records, transforms, "RA", density_config=model_config.density)
    # 评碰撞用当前所选受体的完整标准重原子, 即使某原子不在用于生成的局部口袋内也保留.
    receptor = build_receptor_molecule(read_receptor(inputs / "parse" / sample_id / "receptor_tokens.npz"))
    results = []
    # 原采样器的输出按 split/protocol/PDB/blob 嵌套. 在同一文件系统内完成写盘与自评分, 再一次移动完整的 blob 目录.
    with TemporaryDirectory(prefix=".sampling-", dir=directory) as staging:
        config.output_root = staging
        for index, (request, record) in enumerate(zip(requests, records)):
            result = sample_occurrence(dataset, index, model, noiser, featurizer, config, "local-c1")
            native = Path(staging) / config.split / "local-c1" / sample_id / str(request["blob_id"])
            (native / "scientific_input.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
            candidates = json.loads((native / result["candidate_file"]).read_text(encoding="utf-8"))
            template, poses = None, []
            if result["pose_file"]:
                # 无构象的模板仅确定本次预测身份声明的立体标记, 坐标由原输出 SDF 提供.
                template = read_smiles_graph(str(chemistry), record["prepared_smiles"])["mol"]
                poses = list(Chem.SDMolSupplier(str(native / result["pose_file"]), sanitize=True, removeHs=True))
            metrics = score_saved_candidates(candidates, poses, receptor, template, rmsd_reference=None)
            ranked = rank_candidate_metrics(metrics)
            top = ranked[0]["sample_index"] if ranked else None
            ranking = {"pdb_id": sample_id, "source_blob_index": request["blob_id"], "stage": "local-c1",
                       "candidate_metrics": metrics, "top1_sample_index": top}
            (native / "ranking.json").write_text(json.dumps(ranking, ensure_ascii=False, allow_nan=False, indent=2), encoding="utf-8")
            target = directory / str(request["blob_id"])
            native.rename(target)
            pose_file = (Path(target.name) / result["pose_file"]).as_posix() if result["pose_file"] else None
            results.append({
                "blob_id": request["blob_id"], "docking_rank": index + 1, "status": result["status"],
                "pose_file": pose_file, "top1_sample_index": top,
                "failure_reason": next((c["error"] for c in candidates if c["error"]), None),
            })
    return results
