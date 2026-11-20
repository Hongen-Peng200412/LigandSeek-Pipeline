# 从实验图定位配体

对一个样本运行 `python -m ligandseek infer sample.json --models models --output result --until find` 后，可从 `result/find.json` 找到候选及其中心。每个中心是完整连通区域的体素中心在实验世界坐标中的算术均值，单位为 Å。可从 [公开输入契约](../../SELECTION.md) 核对样本文件和候选字段。

阅读代码从 [find](./__init__.py) 开始。它先调用 [prepare_inputs](./inputs.py) 生成样本资产，再依次调用 Pocket Plus 的 `infer_full_map`、`extract_probability_blobs`、`infer_centered_boxes` 和 `score_centered_candidates`。这四个函数分别计算整图概率、提取过阈值的连通区域、计算区域附近特征和给区域评分。阶段运行时须能导入本次发布随附的 Pocket Plus 和 AdaLigand 预处理库。

```text
## <科学产物>
result/
  find.json                       # 当前样本的区域编号、选择状态、中心和分数
  find/
    probability/probability_map.npz  # 完整实验网格上的配体概率
    blobs/F2_blobs.npz               # 冻结概率阈值产生的全部 26 连通区域

## <其他文件>
result/find/
  inputs/density/<sample>/{exp,sim}.npy      # 单通道实验密度及受体模拟密度（计算中间文件）
  inputs/density/<sample>/{exp,sim}.npz      # 对应密度的网格几何（计算中间文件）
  inputs/parse/<sample>/receptor_tokens.npz  # 当前受体的原子表和键（计算中间文件）
  centered/F2_centered.npz                   # Match 使用的局部特征及预测（计算中间文件）
```

两部分中的 `result/find/` 是同一个目录。每次调用处理一个样本，`<sample>` 是其小写 `sample_id`，如 `9ter`。

## 样本和固定资产

实验图与同一受体的模拟图均可使用 MRC/map（可 gzip 压缩），也可使用模型网格上的 NPY 与同名 NPZ 几何。受体 CIF 使用首模型，按原适配器选定异构位置并提取聚合物重原子。

`sample.json` 的可选字段 `map_convention` 指定两张图的原点解释。省略时，模拟图为 MRC 选 `mrc`，模拟图为 NPY 选 `adaligand`。已经物化的 NPY/NPZ 逐字节传递。

| 取值 | 适用文件和世界坐标 |
| --- | --- |
| `mrc` | 适用于正交实空间密度图。采用 ChimeraX 的采样点约定。非零 `ORIGIN=(x,y,z)` 给出第零采样点的世界坐标（Å）。三项全零时，把列、行、层的 `NSTART` 按 `MAPC/MAPR/MAPS` 转为 XYZ 体素下标，再分别乘三轴间距。若提供者以首体素外角点记原点，导出时应将该角点加上半个原体素间距，得到首采样点。 |
| `adaligand` | 历史 AdaLigand 实验 MRC 调用训练数据的原始导入函数。它将头中的 XYZ 原点加上已转轴的 NSTART，再逐轴乘体素间距。历史模拟 MRC 使用原函数的模拟图分支。既有实验 MRC 配 `sim.npy`、`sim.npz` 的样本自动沿用此约定。 |

两种约定均复用训练时的密度转轴和 Fourier 重采样。`mrc` 会在奇数长度的轴前端补一个原体素。补齐后的第零采样点在重采样前后具有相同世界坐标。为供 Find 用 `origin+(i+0.5)*voxel_size` 定位，写 NPZ 时从该点扣除一半**输出网格实际间距**，得到下角点。

MRC2014 说明和 MRC2000/CCP4 使用文献对原点有外角点、首体素等不同表述。本包的 `mrc` 约定对齐 ChimeraX 的 `GridData.ijk_to_xyz((0,0,0))`。例如三轴间距 `[2,3,4]` Å、`MAPC/MAPR/MAPS=[3,1,2]`、`NSTART=[2,5,-3]` 且 ORIGIN 全零时，首采样点为 `[10,-9,8]` Å。XYZ 下标 `[3,2,1]` 的密度峰在重采样前后均位于 `[16,-3,12]` Å。[CCP-EM MRC2014](https://www.ccpem.ac.uk/mrc-format/mrc2014/) · [Cheng 等，第 4.6 节](https://pmc.ncbi.nlm.nih.gov/articles/PMC4642651/) · [Wriggers，补充材料第 6 节](https://journals.iucr.org/d/issues/2012/04/00/ba5170/ba5170sup1.pdf) · [ChimeraX GridData](https://chimerax.readthedocs.io/en/latest/modules/map/map.html)

模型根下的 `find/config.yaml` 和 `find/model.ckpt` 提供网络结构与权重。`find/runtime.yaml` 提供窗口步长、批量、线程数和设备。Find 依据 `receptor_source` 从 `find/calibration/ca2/` 或 `find/calibration/deposited/` 读取校准。`F2_semantic.json` 的 `threshold_value` 用于提取连通区域。`F2_gaussian.json` 的 `score_parameters`、`score_threshold`、`prefiltered_min_voxel` 和 `min_voxels` 决定区域评分和体积筛选。早期 CA2 模型包的两份文件直接位于 `find/`。`ccd/` 提供受体残基的键和原子名。

模型和化学预处理工具安装在不同环境时，设置 `LIGANDSEEK_PREPARE_PYTHON` 为含 gemmi、RDKit、mrcfile、SciPy 和 NumPy 的 Python 解释器。已有完整环境时 Find 直接使用当前解释器。

### 完整密度和受体

- `exp.npy`、`sim.npy` 是 float32 `(1,Z,Y,X)` 数组。两个 NPZ 的 `origin` 和 `voxel_size` 均为 float32 `(3,)` 世界 XYZ 向量。以体素 ZYX `[2,4,6]` 为例，其中心为 `origin + [6.5,4.5,2.5] * voxel_size`。元数据 `schema_version=2`，实验图还写入整数 `canonical_shape_zyx`，如 `[192,192,192]`。
- `receptor_tokens.npz` 的原子轴长度为 N。`coords` 是 float32 `(N,3)` 世界 XYZ Å，`feat` 是 float32 `(N,49)` 原子特征。`element` 是原子序数，`res_type` 是残基类别，两者为 uint8 `(N,)`。`atom_name` 为 `S4 (N,)` 的原子名，如 `CA`。`is_backbone` 是 bool `(N,)` 的主链标志。`res_index`、`chain_index` 是 int32 `(N,)`，相同的值标识同一残基或链。`bond_index` 的每一列是两个原子的零基编号，形状 int32 `(2,E)`；`bond_type` 为 uint8 `(E,)`，如 0 表示单键。

### 整图概率和连通区域

- `probability_map.npz` 保存 float32 `probability_map (Z,Y,X)`，每个值为对应体素属于配体的预测概率，如 0.72。`origin_xyz`、`voxel_size_xyz` 与实验几何对应。
- `F2_blobs.npz` 以内含端点的概率阈值 `source_threshold_value`（float32，长度 1）提取 B 个区域，按平均概率降序、最小线性体素编号升序排列。`blob_index` 和 `voxel_count` 是长度 B 的 int32 数组。`source_probability_mean` 为 float32 `(B,)`。`voxel_offsets` 是 int64 `(B+1,)`，从 0 到所有区域体素总数，切分全局 ZYX 下标 `voxel_index_global_zyx`（int32，三列）与同序 `source_probability`（float32）。例如 offsets 为 `[0,10,25]` 时第二个区域占 `[10:25]`。`fits_centered_box` 是 bool `(B,)`，表示完整区域能否放入一个 80³ 局部框。`centered_box_start_zyx` 是 int32 `(B,3)` 的全局框起点，无法完整容纳时为 `[-1,-1,-1]`。

### Match 使用的局部特征

`F2_centered.npz` 包含 C 个局部框。C 是达到原前向门槛 8 个体素的区域数，`selected=False` 的局部框也保留。跨框的变长数组都用长度 C+1 的 int64 offsets 切片，首值 0、末值是该组数组的总长度。

| 局部框中的实体 | 索引、字段和数值的解释 |
| --- | --- |
| 区域与局部框 | `centered_box_index`、`source_blob_index`（int32，长度 C）分别是局部框编号和对应的原区域编号。`box_start_zyx`、`box_shape_zyx` 为 int32 (C,3)，框边长为 80。`box_origin_world`、`voxel_size_world` 为 float32 (C,3)，世界 XYZ，单位 Å。`source_probability_mean`、`source_threshold_value`、`source_blob_fits_centered_box`、`source_voxel_count` 均与当前文件的 C 个框同序。例如 `source_blob_index[4] == 17` 时，该框的原概率均值为 `source_probability_mean[4]`，其完整区域在 `F2_blobs.npz` 内的编号为 17。`score`（float32，长度 C）为 Gaussian 分数，`selected`（bool，长度 C）结合分数与体积门槛确定入选框。 |
| 框内区域体素 | `voxel_offsets` 切分 `voxel_index_local_zyx`、`source_probability`、`centered_probability` 和 `voxel_final`。下标是从框起点算起的 int16 三元组，如 `[12,31,5]`；两种 float32 概率分别来自整图和局部前向。`voxel_final` 为二维 float16 特征，第一维按上述体素排列。 |
| 受体占据体素 | `voxel_aux_offsets` 切分 int16 三列局部下标 `voxel_aux_index_local_zyx` 与 float32 `voxel_aux_probability`，每个值对应该占据体素的受体辅助预测。 |
| 受体原子 A | `A_offsets` 切分 `A_global_index`、`A_coord_local_xyz`、`A_coord_centered_world`、`A_probability`、`A_feat_L0` 至 `A_feat_L3`。全局编号指向受体表，如 52；局部 XYZ 坐标单位为体素，相对框中心的 XYZ 位移单位为 Å，均为 float32 三列。`A_probability` 是逐原子概率；L0 为 float32 的 49 维原子特征加主链标志，L1–L3 为 float16 模型特征。 |
| 预测锚点 P | `P_offsets` 切分 `P_coord_local_xyz`、`P_probability`、`P_feat_L2`、`P_feat_L3`。每个锚点有三维局部体素坐标、一个概率及两层模型特征，坐标和概率均为 float32，模型特征为 float16。 |
| 48³ 密度裁块 | `experimental_density_48`、`simulated_density_48`、`source_probability_48` 都为 float32 `(C,48,48,48)`。每框的 `v_centroid_local_zyx` 为框内区域体素下标质心，`crop_start_local_zyx` 为 48³ 裁块相对 80³ 框的起点，`crop_center_offset_zyx` 为质心减去裁块中心。`crop_clipped_axis_mask` 为 bool `(C,3)`，如 `[true,false,false]` 表示沿 Z 轴把裁块移回了框内。 |

体积过滤后 C 可以为零。此时局部特征的候选轴为空，offsets 为 `[0]`，`find.json` 的 `candidates` 为 `[]`。有候选时，面向 Build 的 `center_xyz_A` 始终由 `F2_blobs.npz` 的完整区域计算。
