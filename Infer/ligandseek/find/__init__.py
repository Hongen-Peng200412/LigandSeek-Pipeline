"""在完整实验图上定位配体候选, 按本次样本保存原始概率、连通区域和局部特征.

find 依次调用 Pocket Plus 的整图滑窗、连通区域和 centered 前向. directory 决定三个原生 NPZ 及输入网格的位置.
"""

import json
import random
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def calibration_directory(models: Path, sample: dict) -> Path:
    """返回本次受体的两份冻结校准文件所在目录. 新包使用 find/calibration/ca2 或 deposited, 旧 CA2 模型包使用 find/."""
    root = models / "find"
    source = sample.get("receptor_source", "ca2")
    if source == "ca2" and not (root / "calibration" / "ca2").is_dir():
        return root
    return root / "calibration" / source


# ================================================================================================
def find(sample: dict, models: Path, directory: Path, seed: int) -> list[dict]:
    """跑 Find：对单个样本产出 candidate blobs, 以及相关的特征和文件.

    输入:
        - sample: 含 sample_id、experimental_map、receptor_cif、simulated_map 的样本. 用同一受体和网格执行整图与局部前向, 如 sample_id="9ter".
        - models: 资产根. find/config.yaml、model.ckpt 是 Find_1 的训练配置和权重. runtime.yaml 指定滑窗批量、精度和设备. F2_semantic.json 及 F2_gaussian.json 分别提供冻结概率阈值和 Gaussian 区域选择参数.
        - directory: 本次样本的 Find 目录, 如 results/9ter/find.
        - seed: 整数. 本阶段 Python、NumPy 和 PyTorch 随机种子.

    输出文件:
        - inputs/: 输入图和受体的训练格式, 详见 prepare_inputs.
        - probability/probability_map.npz: probability_map 为 float32 (Z,Y,X) 完整配体概率. origin_xyz、voxel_size_xyz 为 float32 (3,) 世界 XYZ 网格角点和间距, 单位 Å.
        - blobs/F2_blobs.npz: 原概率过阈值的全部 26 连通区域. voxel_offsets (B+1,) 同步切分 voxel_index_global_zyx (V,3) 与 source_probability (V,), B 为区域数, V 为这些区域的完整体素总数. source_probability_mean (B,) 为原区域的平均概率.
        - centered/F2_centered.npz: 上述 B 个区域中达到局部前向最小体积者的原生特征, 按 centered_box_index 指向局部框、source_blob_index 指向完整区域. 增写逐框 score (C,) 与 selected (C,), C 为局部候选数. 文档列出特征和稀疏切片关系.

    返回 C 个区域字典, 同一列表传给 Match:
        - blob_id、centered_box_index: 样本内原区域和局部框的整数编号, 如 17 和 4.
        - center_xyz_A: 长度 3 的浮点列表, 完整区域内体素中心的世界 XYZ 算术均值, 单位 Å.
        - selected、probability_mean: 布尔选择标志及原区域平均概率, 如 False 和 0.74.
        - identity_scores、foreground_probability: 当前为 None, Match 填入身份得分与前景概率.
    """
    import numpy as np
    import torch
    from omegaconf import OmegaConf
    from src.datasets.stage1_dataset import Stage1Dataset
    from src.datasets.stage1_requests import ResolvedStage1Crop
    from src.inference.artifacts import publish_stage1_artifact
    from src.inference.blobs import extract_probability_blobs
    from src.inference.centered import infer_centered_boxes
    from src.inference.checkpoint import load_stage1_wrapper
    from src.inference.full_map import infer_full_map
    from src.inference.scoring import score_centered_candidates

    from .inputs import prepare_inputs

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    sample_id = sample["sample_id"].lower()
    assets = models / "find"
    calibration = calibration_directory(models, sample)
    inputs = directory / "inputs"
    prepare_inputs(sample, inputs, models / "ccd")

    runtime = OmegaConf.load(assets / "runtime.yaml")
    wrapper, training = load_stage1_wrapper(assets / "model.ckpt", assets / "config.yaml", "cpu", "current_workspace")
    wrapper.to(torch.device(runtime.device))
    # 真实滑窗及 centered 起点由两个原生前向函数生成. 内存请求只标识本样本及其首个滑窗, require_targets=False.
    first_window = ResolvedStage1Crop(pdb_id=sample_id, box_start_zyx=(0, 0, 0), require_targets=False, role="sliding")
    dataset_options = OmegaConf.to_container(training.dataset, resolve=True)
    dataset_options.pop("_target_")
    dataset_options.update(all_data_path=str(inputs), split_file=[first_window], mode="full_map", box_pool_root=None, enable_random_rotation=False)
    dataset = Stage1Dataset(**dataset_options)

    windows = runtime.window
    full = infer_full_map(dataset, dataset.collate_fn, wrapper, sample_id, str(runtime.device),
                          windows.stride_zyx, windows.gaussian_sigma, windows.batch_size, windows.workers,
                          windows.prefetch_batches, windows.precision, windows.pending_fusion_batches)
    # (Z,Y,X) 完整概率供本次连通区域提取和后续体素评估; 两个 (3,) 世界 XYZ 向量给出完整图几何.
    publish_stage1_artifact(directory / "probability" / "probability_map.npz", {
        "probability_map": full.probability_map, "origin_xyz": full.origin_xyz, "voxel_size_xyz": full.voxel_size_xyz,
    }, None)
    semantic = json.loads((calibration / "F2_semantic.json").read_text(encoding="utf-8"))
    blobs = extract_probability_blobs(full.probability_map, float(semantic["threshold_value"]))
    publish_stage1_artifact(directory / "blobs" / "F2_blobs.npz", blobs, None)

    # 原 centered 前向的最小体素数为 8. 最终 selected 还要求区域体积和 Gaussian 分数分别过冻结阈值.
    selection = json.loads((calibration / "F2_gaussian.json").read_text(encoding="utf-8"))
    local = runtime.centered
    with ThreadPoolExecutor(max_workers=1) as packer:
        packed, _ = infer_centered_boxes(dataset=dataset, collator=dataset.collate_fn, wrapper=wrapper, pdb_id=sample_id,
                                         producer="Find_1", blobs=blobs, full_probability=full.probability_map,
                                         origin_xyz=full.origin_xyz, voxel_size_xyz=full.voxel_size_xyz, device=str(runtime.device),
                                         precision=str(local.precision), centered_batch_size=local.batch_size,
                                         centered_workers=local.workers, prefetch_batches=local.prefetch_batches,
                                         pending_cpu_batches=local.pending_cpu_batches, forward_min_voxels=8, packer=packer)
        centered = packed.result()
    centered["score"] = score_centered_candidates(centered, score_mode=selection["score_mode"], score_parameters=selection["score_parameters"]).astype(np.float32)
    # (C,), True 表示进入 Find 最终选择集合. 未选中的局部候选仍参与 Match.
    count = centered["source_voxel_count"]
    centered["selected"] = (centered["score"] >= np.float32(selection["score_threshold"])) & (count >= selection["prefiltered_min_voxel"]) & (count >= selection["min_voxels"])
    publish_stage1_artifact(directory / "centered" / "F2_centered.npz", centered, None)

    candidates = []
    for index, blob in enumerate(centered["source_blob_index"]):
        # (V_blob,3), 此原区域的完整 ZYX 体素索引. 以所有体素算中心, 包括未进入局部框的体素.
        start, stop = blobs["voxel_offsets"][blob:blob + 2]
        voxels = blobs["voxel_index_global_zyx"][start:stop]
        center_xyz = full.origin_xyz + (voxels.mean(axis=0)[::-1] + 0.5) * full.voxel_size_xyz
        candidates.append({
            "blob_id": int(blob), "centered_box_index": int(centered["centered_box_index"][index]),
            "center_xyz_A": center_xyz.tolist(), "selected": bool(centered["selected"][index]),
            "probability_mean": float(blobs["source_probability_mean"][blob]), "identity_scores": None, "foreground_probability": None,
        })
    return candidates
