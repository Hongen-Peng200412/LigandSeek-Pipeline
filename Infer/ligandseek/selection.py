"""将一个样本的 Find 区域和 Match 身份转成 Build 计划.

plan_docking 返回全部候选区域, 包括跳过的区域. 用原 blob_id 连接 Find 和 Match.
返回列表按考察顺序排列, planned_docking_rank 只给计划对接的区域编号.
"""


def matched_smiles(candidate: dict, smiles_axis: list[str]) -> str | None:
    """从单个区域的 Match 得分中取第一名. I 是本样本不同候选身份的数量.

    - candidate["identity_scores"]: (I,), 与 smiles_axis 同序的原始得分, 如 [1.2, 0.7]. Find 单独运行时为 None.
    - smiles_axis: (I,), 用户给定的精确 SMILES, 如 ["CCO", "CCN"]. 同分取列表中靠前的身份.
    - 返回值: 最高分的 SMILES, 如 "CCO". 未匹配或身份轴为空时返回 None.
    """
    scores = candidate["identity_scores"]
    if not scores:
        return None
    # 无法从某个 SMILES 构图时, Match 对应槽为 None. 数值同分时 max 保留靠前的身份.
    available = [i for i, score in enumerate(scores) if score is not None]
    return smiles_axis[max(available, key=scores.__getitem__)] if available else None


# ================================================================================================
def plan_docking(
    candidates: list[dict], identities: list[dict], *, k: int | None, scope: str,
    foreground_threshold: float | None, decisions: list[dict], seed: int,
) -> list[dict]:
    """按 Find 的考察顺序, 为一个样本的每个区域选择身份并分配对接预算.

    输入参数:
        - candidates: 长度 B 的区域列表, B 是 Find 区域数. 每项的 blob_id 是原区域编号, 如 17.
            - center_xyz_A: (3,), 完整区域的体素中心均值, 实验图世界 XYZ 坐标, 单位 Å.
            - selected: bool, Find 的 Gaussian 选择标记. True 的区域优先考察.
            - probability_mean: 浮点数, 原区域内的平均概率, 如 0.74.
            - identity_scores: (I,), I 是 identities 中的身份数, 同位置属于同一身份. Find 单独运行时为 None.
            - foreground_probability: 浮点数或 None, Match 判断为配体前景的概率, 如 0.31.
            - reference: 可选的归属实例字典. occurrence_id 指向具体拷贝, smiles 标识该拷贝, 各资格为 True、False 或未知.
        - identities: 长度 I 的身份列表, 每项包含精确 smiles 和可选的 type_tag, 如 "small_molecule".
        - k: 单样本最多安排的请求数. None 表示不限, 0 仍然保留全部区域.
        - scope: wide 按 small molecule 筛选, narrow 加入同一实例的已知资格, all 接受各类别.
        - foreground_threshold: 自动选择的可选前景概率下限. None 关闭, 等于阈值或概率未知时放行.
        - decisions: 按 blob_id 指定 action 的列表, 如 [{"blob_id":17,"action":"dock","smiles":"CCN","seed":123}], 可能为空表示全自动.
        - seed: uint32 范围内的 Build 基础种子, 如 23.

    返回值:
        - 计划列表: 长度 B, 每项保留原区域字段, 并增加下列字段.
            - probability_rank: 从 1 起的原概率名次, 同分取较小 blob_id.
            - inspection_rank: 从 1 起的考察名次, selected 分组优先, 组内沿用原概率顺序.
            - predicted_smiles: 模型最高分身份, 如 "CCO". 未得到身份时为 None.
            - chosen_smiles: 选给 Build 的身份, 如手选的 "CCN". 手动跳过或没有身份时为 None.
            - manual_input_smiles: 手工指定的 SMILES 原文, 如 "NCC". 没有人为指定身份时为 None.
            - selection_source: automatic 或 manual. 仅覆盖种子的 auto 仍为 automatic.
            - plan_action: dock 表示计划对接, skip 表示跳过.
            - reason: 决定动作的原因, 如 budget_exhausted 表示预算已用完.
            - engine: 计划使用的引擎名 dg, 跳过时为 None.
            - sampling_seed: 本请求实际安排的采样种子, 跳过时为 None.
            - planned_docking_rank: 从 1 起的计划提交号, 跳过时为 None. Build 另记实际提交号.
    """
    smiles_axis = [identity["smiles"] for identity in identities]  # (I,), 顺序与每个区域的身份得分一致.
    identity_types = {identity["smiles"]: identity.get("type_tag") for identity in identities}
    decisions_by_blob = {decision["blob_id"]: decision for decision in decisions}

    # 两个长度 B 的排列. 稳定排序只把 Find 选中的一组移到前面, 保留各组原概率名次.
    by_probability = sorted(candidates, key=lambda blob: (-blob["probability_mean"], blob["blob_id"]))
    probability_ranks = {blob["blob_id"]: rank for rank, blob in enumerate(by_probability, 1)}
    inspection_order = sorted(by_probability, key=lambda blob: not blob["selected"])

    plan = []
    planned_count = 0
    for inspection_rank, blob in enumerate(inspection_order, 1):
        blob_id = blob["blob_id"]
        # 显式给出的身份优先. predicted 始终保存自动模型答案, 供后续独立计费.
        predicted = matched_smiles(blob, smiles_axis)
        decision = decisions_by_blob.get(blob_id, {})
        action = decision.get("action", "auto")
        chosen = decision.get("smiles", predicted)

        # 只有所选身份正好对应本区域归属的真实拷贝, 才能在窄范围使用该拷贝的沉积重原子完整性等条件.
        occurrence = blob.get("reference") or {}
        if scope != "narrow" or occurrence.get("occurrence_id") is None or occurrence.get("smiles") != chosen:
            occurrence = {}
        identity_type = identity_types.get(chosen) or occurrence.get("type_tag")
        foreground = blob["foreground_probability"]

        if action == "skip":
            chosen = None
            reason = "manual_skip"
        elif chosen is None:
            reason = "no_identity"
        elif action == "dock":
            reason = "manual_dock"
        elif foreground_threshold is not None and foreground is not None and foreground < foreground_threshold:
            reason = "below_foreground_threshold"
        elif scope != "all" and identity_type not in (None, "small_molecule"):  # 基于 identity_type 的过滤
            reason = "known_non_small_molecule"
        else:
            reason = "eligible"  # 基于 occurrence 的过滤
            for rule in ("eligible", "single_residue_ccd", "noncovalent", "complete_heavy_atoms", "template_compatible", "language_compatible"):
                if occurrence.get(rule) is False:
                    reason = "reference_failed_" + rule
                    break

        # 一次请求使用一个预算单位. 跳过者保留考察名次, 对接号与采样种子为空.
        docking_rank = sampling_seed = None
        if reason in ("eligible", "manual_dock"):
            if k is None or planned_count < k:
                planned_count += 1
                docking_rank = planned_count
                sampling_seed = decision.get("seed", (seed + blob_id) % 2**32)
            else:
                reason = "budget_exhausted"

        plan.append({
            **blob,
            "probability_rank": probability_ranks[blob_id],
            "inspection_rank": inspection_rank,
            "planned_docking_rank": docking_rank,
            "predicted_smiles": predicted,
            "chosen_smiles": chosen,
            "selection_source": "manual" if "smiles" in decision or action != "auto" else "automatic",
            "manual_input_smiles": decision.get("input_smiles", decision.get("smiles")),
            "plan_action": "dock" if docking_rank is not None else "skip",   # docking 或跳过
            "reason": reason,

            "engine": "dg" if docking_rank is not None else None,
            "sampling_seed": sampling_seed,
        })
    return plan
