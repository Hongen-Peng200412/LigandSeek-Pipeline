# 从预测区域和指定 SMILES 生成姿态

[build](./__init__.py) 执行公开推理命令中已排定的对接请求。先看 [Infer 的命令示例](../../README.md)，再读[请求与姿态字段](../../SELECTION.md)。例如请求中 `blob_id=17`、`chosen_smiles="CCO"`、`center_xyz_A=[12.5,-3.0,8.25]` 时，引擎以该世界坐标为中心选择受体口袋，使用乙醇的规范 SMILES、固定种子和局部密度生成构象。

代码按实际执行顺序排列。函数开始复用本次 Find 准备的地图与受体，根据用户所选身份生成分子图。随后加载 `models/build/train.yml` 和 `models/build/step=24800.ckpt`，用原训练模型、`ConditionedDockingDataset` 与 `sample_occurrence` 对每个请求采样 50 个姿态，每个姿态 100 步。最后用 `score_saved_candidates` 计算原置信度、碰撞与立体检查的综合评分，按 `rank_candidate_metrics` 的顺序选择首位。

下列路径均相对一次推理的输出目录。树中的 `build/17/` 是一个实际目录，分两段列出其中的科学文件和采样记录。

```text
## <科学产物>
build.json                    # 每个实际提交请求的区域编号、提交次序、状态和首位姿态
build/17/
  scientific_input.json        # 当前请求实际使用的身份、位点中心及种子
  candidates.json              # 按原采样编号保存包括失败在内的全部姿态记录
  poses.sdf                    # 成功构建的重原子构象，实验世界 XYZ，Å
  confidence.npz               # 成功姿态对应的模型原置信度
  ranking.json                 # 逐姿态综合评分与首位姿态编号
## <其他文件>
build/17/result.json           # 原采样器的耗时、显存与完成计数（运行统计）
build/chemistry/               # 独立 Build 或手选集合外身份时生成的图和对称排列（可再生工具输入）
build/inputs/                  # 独立运行 Build 且需重新准备地图与受体时生成（可再生工具输入）
```

### 请求身份及提交结果

`scientific_input.json` 的 `candidate_id` 为 Find 区域编号，`pdb_id` 为本次样本标识，`prepared_smiles`、`given_center_xyz_A` 和 `sampling_seed` 记录实际用于采样的化学身份、三维中心和该请求的随机种子。上层写出的 `build.json` 中，`docking_rank=1` 指本样本第一个实际提交。该请求失败时仍保留其名次。全请求失败记录例如 `{"blob_id":17,"docking_rank":1,"status":"failed","pose_file":null,"top1_sample_index":null,"failure_reason":"ValueError: unsupported_smiles_element_or_bond"}`，这里仅演示字段格式。

### candidates.json、poses.sdf 和 confidence.npz

`candidates.json` 是按原 `sample_index` 排列的 50 个记录组成的列表。`sdf_index` 为成功姿态在 SDF 中的零基位置。例如前两个姿态失败时，第 3 个姿态的 `sample_index=2`、`sdf_index=0`。SDF 中每个分子的原子顺序与所选 SMILES 图对应。置信度 NPZ 的首轴依 SDF 排列，`sample_index` 数组保留对应的原采样编号。`confidence_pos_traj` 的三轴是成功姿态数、配体重原子数、采样步数。没有成功姿态时 `result.json` 的 `pose_file` 和 `confidence_file` 均为空。

### ranking.json

`candidate_metrics` 与 `candidates.json` 顺序相同，新增碰撞检查 `no_clashes`、立体检查 `stereo`、置信度与这两项布尔指示之和 `self_ranking`。无可用评分时 `self_ranking=null`。取前 K 个姿态时先按综合评分降序排列，同分取较小原采样编号。`top1_sample_index` 是首位的原编号，例如 `2`。推理阶段的 `rmsd_A` 为 `null`。

### result.json

`status` 是整次采样的 `success`、`partial` 或 `failed`。`success_count` 和 `failed_count` 的和为 `num_candidates`。生成统计还包括 `num_steps`、`batch_size`、模型原点 `model_origin_world_xyz_A`、推理秒数 `inference_seconds` 和 GPU 峰值 `peak_memory_allocated_bytes`。完整的原采样记录与每个姿态的错误信息均保留在对应区域目录。
