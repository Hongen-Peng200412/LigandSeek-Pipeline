# 从 Find 区域预测前景和配体身份

`ligandseek.match.match` 接收同一轮 Find 的输出、所选受体和样本的候选身份列表。公开命令和完整样本示例见[推理用法](../../README.md)，阶段文件的字段见[推理契约](../../SELECTION.md)。

按以下顺序读实现。

1. [match](./__init__.py) 按用户的规范 SMILES 构造身份轴，从 `models/match/BEST.pt` 读取模型配置、权重和原验证集前景阈值。每批最多 32 个同一样本的区域，前向沿用原 bfloat16 设置。
2. [prepare_graphs](./chemistry.py) 用 `LIGANDSEEK_PREPARE_PYTHON` 指定的 RDKit 解释器调用 Matcher 的 `build_arrays`，得到原子、键、环标记和原子对称排列。省略该环境变量时沿用当前解释器。
3. [candidate_arrays 和 collate_inputs](./batch.py) 从完整 Find 局部特征和本次受体组装 106 通道网格，再连接区域 token 与身份图。比如 `centered_box_index=4` 表示在局部框编号数组中查找值 4 的位置。图内的原子下标只在批次连接时加偏移。

下列路径均相对一次推理的输出目录。

```text
## <科学产物>
match.json                    # 完整用户身份轴及每个 Find 区域的得分与前景概率
match_model.json              # 验证集前景阈值、进入前向的身份位置及排除原因
## <其他文件>
match/chemistry/
  smiles_graphs_v1.npz         # 按当前输入重建的分子图（可再生工具输入）
  smiles_symmetries_v1.npz     # 分子图对应的原子等价排列（可再生工具输入）
```

### match.json

上层推理入口写一个样本的 `identity_smiles` 和 `candidates`。每个区域的 `foreground_probability` 是前景 logit 的 sigmoid，`identity_scores` 保留用户给出的全部身份槽。举例来说，`[0.2, null, 1.4]` 表示中间的身份没有可用得分，第一和第三项可参与身份选择。空身份表仍能输出区域前景概率。

### match_model.json

`checkpoint` 为实际参数文件路径。`foreground_threshold` 为原验证集选出的固定分类阈值，例如 `0.359375`。`valid_slots` 是有效身份在规范身份轴中的零基位置，`[0, 2]` 表示使用第一个和第三个身份。`unsupported_identities` 逐个给出无法编码的规范 `smiles` 和 `error`，例如构造示例 `{"smiles":"CCO","error":"missing SMI-TED embedding"}`。

### smiles_graphs_v1.npz 与 smiles_symmetries_v1.npz

两份数组均按成功解析的规范 SMILES 字符串排序。`atom_offsets` 和 `bond_offsets` 从 0 到全表原子数和无向键数，分别切分逐原子和逐键特征。`bond_index` 形状为 `(2, 总键数)`，每列是一个键在该分子内的两个零基原子编号。`matches_offsets` 从 0 到展平排列元素总数，切出的区间按照该身份的 `matches_shape` 恢复为二维整数数组。例如 `[1,0]` 表示只有恒等排列，没有需要记录的可交换原子列。
