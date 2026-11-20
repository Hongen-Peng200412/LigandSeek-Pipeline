# 从实验密度图预测配体位点、身份和姿态

完整 Linux 安装提供 `ligandseek run raw.json --output result`。它先在 `result/prepared/` 写准备产物，再将三个模型的预测写到 `result/predictions/`。只准备本次样本可用 `ligandseek prepare raw.json --output prepared`；原始输入与准备产物见 [Prepare](../Prepare/README.md)。

已有实验图、同网格模拟图、受体和身份向量时，可以单独使用本目录。在本目录安装 `python -m pip install .`，或直接在本目录运行模块。准备样本 JSON 和[模型资产](MODELS.md)后，以下命令将三步结果写入新的 output 目录。

```sh
python -m ligandseek infer sample.json --models weights --output output
```

Find 在实验图上生成完整区域和局部特征，Match 对用户列出的每种可编码身份打分，Build 按计划使用预测位点和所选 SMILES 生成姿态。要先审阅候选，可以使用 `--until match`。仅制定计划时，Python 3.10 或更新版本可调用本机 RDKit 或由 `LIGANDSEEK_PREPARE_PYTHON` 指定的化学解释器。

## 将区域 17 从乙醇改选为乙胺

把下面三个构造示例放到当前目录。`plan` 只需要这三个 JSON，可以直接照此运行。

**sample.json** 描述实验输入和可选配体身份。列表顺序规定 Match 得分的顺序。

```json
{
  "sample_id": "demo-ca2",
  "experimental_map": "data/experiment.mrc",
  "receptor_cif": "data/ca2.cif",
  "simulated_map": "data/ca2_density.mrc",
  "identities": [
    {"smiles": "CCO", "embedding": "data/ethanol.npz", "type_tag": "small_molecule"},
    {"smiles": "CCN", "embedding": "data/ethylamine.npz", "type_tag": "small_molecule"}
  ]
}
```

**match.json** 是已经得到的逐区域预测。区域 17 在 `CCO` 和 `CCN` 上分别得 1.2 和 0.7，自动预测为 `CCO`。

```json
{
  "sample_id": "demo-ca2",
  "identity_smiles": ["CCO", "CCN"],
  "seed": 0,
  "candidates": [{
    "blob_id": 17, "centered_box_index": 4,
    "center_xyz_A": [12.5, -3.0, 8.25], "selected": true, "probability_mean": 0.74,
    "identity_scores": [1.2, 0.7], "foreground_probability": 0.31, "predicted_smiles": "CCO"
  }]
}
```

**decisions.json** 通过 Find 的原区域编号指定动作。这里选择另一候选身份 `CCN`，同时固定这次对接的种子。

```json
[{"blob_id": 17, "action": "dock", "smiles": "CCN", "seed": 123}]
```

```sh
python -m ligandseek plan sample.json --candidates match.json --decisions decisions.json --k 1 --output reviewed
```

打开 `reviewed/plan.json`，找到 `blob_id=17`。它的 `predicted_smiles` 是 `CCO`，`chosen_smiles` 是 `CCN`，`planned_docking_rank` 是 1，`sampling_seed` 是 123。`reviewed/run.json` 保存输入路径和实际参数。若要比较自动选择，省略 `--decisions` 并选一个新输出目录。

同一个入口也可在 Python 中调用。

```python
from pathlib import Path
from ligandseek.pipeline import run_sample

run_sample(Path("sample.json"), Path("reviewed-python"), command="plan",
           candidates=Path("match.json"), decisions=Path("decisions.json"), k=1)
```

## 选择入口和参数

| 命令 | 读取什么 | 停在哪里 |
| --- | --- | --- |
| prepare | 一个原始样本的 raw.json | 写出样本准备目录及可单独递交推理的 sample.json。 |
| run | 同样的原始输入 | 先写 prepared/，随后写 predictions/ 中的指定阶段结果。 |
| plan | 本样本的 Find 或 Match 候选 | 写出全部区域的计划。 |
| build | 同样的已保存候选及模型目录 | 先写计划，再提交计划中的请求。每次实际提交保留种子、状态及可用的构象。 |
| infer | 实验图、受体、同网格模拟图、候选身份和冻结模型 | 依次进行 Find、Match、选择和 Build，可用 --until 提前结束。 |

Python 参数使用表中的同名下划线写法，例如 `--foreground-threshold` 对应 `foreground_threshold`。路径传 `pathlib.Path`。Python 默认为 `command="infer"`，成功写出指定阶段后返回 `None`。

| 参数 | 适用命令 | 取值、默认值和作用 |
| --- | --- | --- |
| sample | 全部 | 必填，单样本输入文件，如 sample.json。文件内相对路径以该文件所在目录为基准。 |
| --output | 全部 | 必填，输出目录，如 reviewed。命令行的相对路径从当前目录解析。 |
| --candidates | plan、build | 必填，同一样本的 Find 或 Match JSON，如 match.json。已有得分按该文件的 identity_smiles 解释，再映射到当前样本的规范身份顺序。 |
| --models | 全部命令 | 权重、配置与科学库所在目录。安装环境通过 LIGANDSEEK_MODELS 指定默认位置，可用此选项覆盖。plan 可省略模型目录；传入时会识别未标注类别的候选。 |
| --until | infer、run | 最后运行的阶段，取 find、match 或默认的 build。 |
| --scope | 选择请求时 | 默认 wide 按已知小分子类别筛选。narrow 再看同一参考实例已知的 Stage3 资格。all 接受各类别。缺少某项判据时放行。 |
| --k | 选择请求时 | 非负整数，每个样本最多安排几次对接。省略表示不限，0 表示只考察和排序。 |
| --foreground-threshold | 选择请求时 | 可选，[0,1] 内的前景概率下限，如 0.4。默认关闭；相等或概率未知时通过。 |
| --decisions | 选择请求时 | 可选，按区域指定 skip、dock 或 auto。dock 接受候选列表里的任意身份，也可直接提供新的 SMILES，先规范化再对接，仍占同一预算。 |
| --reference | narrow 选择时 | 可选，本样本参考类别、实例资格及区域归属，如 reference.json。每次审阅显式选用当次参考。 |
| --find-seed | infer、run | 默认 0，控制 Find 的随机步骤，如 7。 |
| --match-seed | infer、run | 默认 0，控制 Match 的随机步骤，如 11。 |
| --build-seed | 选择请求时 | 默认 0，生成各区域的采样种子，如 23。决定文件中的单区域种子优先。 |

先保存 Match，再根据逐区域预测手选身份。也可在 Find 后直接手选身份，使用 `--until find` 和 `review/find.json`。

```sh
python -m ligandseek infer sample.json --models weights --until match --output review
python -m ligandseek build sample.json --models weights --candidates review/match.json --decisions decisions.json --output poses
```

输入坐标、身份轴、各次序及空值见[字段契约](SELECTION.md)。完整阶段产物分别见 [Find 网格与区域](ligandseek/find/README.md)、[Match 身份得分](ligandseek/match/README.md)和 [Build 姿态](ligandseek/build/README.md)。环境与权重见[模型来源](MODELS.md)。

命令行将单样本参数交给 `run_sample`。它依次通过 `runtime.execute` 在对应环境中调用 Find、Match 和 Build，并在 Build 前根据区域、身份和用户选择确定对接次序。[命令行](ligandseek/__main__.py) · [单样本流程](ligandseek/pipeline.py) · [执行环境](ligandseek/runtime.py) · [对接选择](ligandseek/selection.py)
