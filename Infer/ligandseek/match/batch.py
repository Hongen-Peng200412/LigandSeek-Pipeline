"""为 strongest-1 构造一个受体条件下的区域批次.

candidate_arrays 产生单个区域的局部体素与四类 token. collate_inputs 用原补齐函数组装区域并把同一样本的身份图复制到每个区域, 返回模型前向需要的张量字典.
"""

import numpy as np
import torch

from matcher.data.density import build_density_channels
from matcher.datasets.batching import _pad_float, _pad_integer
from matcher.datasets.candidate_data import (
    _build_full_a_tokens, _build_p_tokens, _build_pp_tokens, _build_v_tokens,
    _manual_features, _scatter_receptor,
)


def candidate_arrays(centered, position, receptor_coords, receptor_features):
    """构造一个预测区域的 strongest-1 输入.

    - centered: 当前 Find 的 F2_centered.npz 原数组. position 索引 centered_box_index 所标记的一个区域.
    - receptor_coords: (N, 3), N 个受体原子在实验图中的世界 XYZ 坐标, Å. receptor_features 为 (N, 50), 前 49 列为原子特征, 末列为主链标志.
    - 返回 voxel_input 为 float16 (106, 48, 48, 48), 空间三轴 ZYX. a_tokens、p_tokens、v_tokens、pp_tokens 的末维为 230、116、68、4, 首维各为该区域原子、P 点、来源 blob 体素、概率点的数量.
    - v_index_zyx 和 pp_index_zyx: int64 (对应点数, 3), 在 48³ 网格中的 ZYX 下标.
    - candidate_manual_features: (12,), 依次为 A、P、V、PP 四组的数量、概率均值和最大值. candidate_sampling_features 为 (4,), 依次记录框裁切、中心偏移 L2 和 Linf 长度、保留的完整 blob 体素比例.
    """
    crop = centered["crop_start_local_zyx"][position].astype(np.int64)
    # (50, 48, 48, 48) 和 (48, 48, 48), 所选受体的原子特征和原子占据掩膜, 沿 ZYX 对齐本预测区域.
    grid, occupied = _scatter_receptor(
        receptor_coords, receptor_features, centered["box_origin_world"][position].astype(np.float32),
        centered["voxel_size_world"][position].astype(np.float32), crop,
    )
    # (56, 48, 48, 48), 固定密度通道由本样本实验图、受体模拟图和受体占据掩膜共同构造.
    density = build_density_channels(
        centered["experimental_density_48"][position], centered["simulated_density_48"][position], occupied,
    )
    atom = _build_full_a_tokens(centered, position, crop)
    pocket = _build_p_tokens(centered, position, crop)
    voxels, voxel_indices = _build_v_tokens(centered, position, crop)
    # 原 strongest-1 配置保留局部概率排名前 2048 个空间点.
    probability, probability_indices = _build_pp_tokens(centered["source_probability_48"][position], 2048)
    offset = centered["crop_center_offset_zyx"][position]
    return {
        "voxel_input": np.concatenate((density, grid)).astype(np.float16),
        "a_tokens": atom, "p_tokens": pocket, "v_tokens": voxels, "pp_tokens": probability,
        "v_index_zyx": voxel_indices, "pp_index_zyx": probability_indices,
        "candidate_manual_features": _manual_features(atom, pocket, voxels, probability),
        "candidate_sampling_features": np.asarray([
            np.any(centered["crop_clipped_axis_mask"][position]), np.linalg.norm(offset),
            np.max(np.abs(offset)), len(voxels) / int(centered["source_voxel_count"][position]),
        ], dtype=np.float32),
    }


# ================================================================================================
def collate_inputs(regions, graphs, embeddings):
    """按原候选顺序拼接一个样本的 B 个区域. graphs 和 float32 (I,768) embeddings 均沿本样本有效身份顺序排列.

    - voxel_input: (B,106,48,48,48). candidate_manual_features、candidate_sampling_features 依次为 (B,12)、(B,4).
    - 各模态 tokens: (B, 该模态批内最多点数, 特征数). 同名 mask 的 True 表示真实点. v_index_zyx 和 pp_index_zyx 的末维为三轴体素下标.
    - language_embeddings、identity_mask: (B,I,768)、(B,I). 同一样本内 B 个区域使用相同的 I 个有效身份.
    - graph_to_sample、graph_to_identity: (B*I,), 每个图在区域轴与身份轴的下标. graph_node_graph_index 是每个重原子在连接图轴的归属.
    - graph_element、graph_charge、graph_atom_in_ring: 所有图的原子序数、形式电荷、3 至 6 元环标记. graph_edge_index 为 (2, 所有有向键数), 端点索引拼接后的原子轴, graph_edge_features 为同序的 9 维键特征.
    """
    batch = {key: torch.from_numpy(np.stack([region[key] for region in regions]))
             for key in ("voxel_input", "candidate_manual_features", "candidate_sampling_features")}
    for name in ("a", "p", "v", "pp"):
        batch[f"{name}_tokens"], batch[f"{name}_mask"] = _pad_float([r[f"{name}_tokens"] for r in regions])
    for name in ("v_index_zyx", "pp_index_zyx"):
        batch[name], _ = _pad_integer([r[name] for r in regions])

    # B 个区域共享此样本的 I 个候选身份, 向量只沿区域轴广播. I=0 时只运行模型的区域前景分支.
    count, identities = len(regions), len(graphs)
    batch["identity_mask"] = torch.ones((count, identities), dtype=torch.bool)
    batch["language_embeddings"] = torch.from_numpy(embeddings).expand(count, -1, -1)
    if graphs:
        # 原 collate_candidates 按区域再按身份连接各图. 节点与有向边顺序保持单图读取器的顺序.
        copies = graphs * count
        nodes = np.asarray([len(graph.element) for graph in copies])
        starts = np.cumsum(nodes) - nodes
        batch["graph_element"] = torch.from_numpy(np.concatenate([g.element for g in copies]))
        batch["graph_charge"] = torch.from_numpy(np.concatenate([g.charge for g in copies]))
        batch["graph_atom_in_ring"] = torch.from_numpy(np.concatenate([g.atom_in_ring for g in copies]))
        batch["graph_edge_features"] = torch.from_numpy(np.concatenate([g.edge_features for g in copies]))
        batch["graph_edge_index"] = torch.from_numpy(np.concatenate(
            [g.edge_index + start for g, start in zip(copies, starts)], axis=1,
        ))
        batch["graph_node_graph_index"] = torch.from_numpy(np.repeat(np.arange(len(copies), dtype=np.int64), nodes))
        batch["graph_to_sample"] = torch.arange(count).repeat_interleave(identities)
        batch["graph_to_identity"] = torch.arange(identities).repeat(count)
    return batch
