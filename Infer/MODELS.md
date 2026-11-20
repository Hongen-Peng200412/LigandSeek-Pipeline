# 模型、运行库与单样本输入

`--models` 指向一份可供所有样本共用的本地目录。推理从这份目录加载原训练权重、冻结校准和当前科学库。样本本身的地图、受体与身份向量通过 [sample.json](SELECTION.md) 提供。

| 模型目录下的文件 | 读取者与用途 |
| --- | --- |
| `find/config.yaml`、`find/model.ckpt` | Pocket Plus 按 Find_1 的训练配置加载原 `TOP_epoch_04_score_0.6654.ckpt`，计算整图概率和局部特征。 |
| `find/runtime.yaml`、`find/calibration/{ca2,deposited}/F2_semantic.json`、同目录 `F2_gaussian.json` | 滑窗、局部批量和精度，生成完整区域的语义概率阈值，选择 Find 定位结果的 Gaussian 评分及阈值。 |
| `match/BEST.pt` | Matcher strongest-1 权重，按用户给出的精确 SMILES 顺序预测前景概率和候选身份得分。 |
| `build/train.yml`、`build/step=24800.ckpt` | DG-PocketXMol Centre 使用全原子受体和局部密度，从预测中心采样。每个实际提交采样 50 个姿态，每姿态 100 步。 |
| `libraries/{Pocket_Plus,Matcher,PocketXMol,AdaLigand}/` | 保留各原库的包目录及科学配置，分别供 Find、Match、DG 采样和地图／化学准备调用。 |
| `ccd/`、`assets.json` | 受体准备使用的通用 CCD 化学组分字典；研究用资产包另记录上述文件的原始来源。 |

`sample.receptor_source` 决定使用哪一对阈值文件，如 `deposited` 读取沉积受体冻结的校准。已有 CA2 资产包将这两个 JSON 放在 `find/` 时仍可读取。

启动 `python -m ligandseek` 的环境承担 Find 和 Match，需要 Pocket Plus 和 Matcher 的运行依赖。受体与地图预处理可用环境变量 `LIGANDSEEK_PREPARE_PYTHON` 指向具备 AdaLigand、Gemmi、RDKit 等依赖的解释器。`LIGANDSEEK_BUILD_PYTHON` 可指向安装 PocketXMol 运行依赖的解释器。两项缺省时都使用启动推理的 Python。

发行版安装时设置 `LIGANDSEEK_MODELS`，每个样本可共用这份资产。用户显式提供的 `--models` 优先。每个调用公开包的 Python 环境安装同版 Infer。

`ccd/<组分号>.pkl` 保存本地 CCD 原子名、键和参考构象，例如 `ATP.pkl`。当本样本出现尚未缓存的组分时，化学环境从同目录的固定 `components.cif.gz` 读取组分，调用 pdbeccdutils 生成同格式模板并写入缓存。词典需收录该组分，缓存目录需允许写入新模板。

一个新样本可以使用 CryoAtom2 给出的受体及其按实验图网格、分辨率生成的模拟密度。为每个候选 SMILES 提供 SMI-TED 表征。模型的图特征、受体体素、局部候选和各次采样输入均从本次样本构造。已知另一受体来源时，可显式指定受体和它在相同实验网格上的模拟密度。
