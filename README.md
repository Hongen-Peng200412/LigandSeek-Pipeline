# LigandSeek

LigandSeek 从冷冻电镜实验密度图定位配体结合区域，比较用户提供的候选 SMILES，再用所选身份在预测位点生成并排列三维姿态。

## 用自己的样本运行

在 Linux x86-64、NVIDIA GPU 服务器上克隆本仓库，按[安装说明](install/README.md)配置解释器、权重和外部工具。安装器可以复用已有软件。使用默认安装位置时，运行以下命令加载环境。

```sh
source ~/.local/ligandseek/environment.sh
```

在与数据同一目录创建 `raw.json`。下面的构造示例提供用户已知的受体和两种候选身份，地图与受体坐标位于同一世界坐标系。

```json
{
  "sample_id": "demo",
  "experimental_map": "experiment.mrc",
  "receptor_cif": "receptor.cif",
  "resolution": 3.0,
  "smiles": ["CCO", "CCN"]
}
```

```sh
ligandseek run raw.json --output result
```

已有蛋白质或核酸序列时，可以用 `protein_fasta`、`rna_fasta`、`dna_fasta` 替代 `receptor_cif`。这些字段也可全部省略，由 CryoAtom2 从地图预测受体。候选 SMILES 在 Prepare 中统一去除显式立体标记，原始写法保存在准备记录中。

| 查看的结果 | 文件与用途 |
| --- | --- |
| 输入和准备产物 | `result/prepared/sample.json` 关联本次受体、同网格实验图与模拟图、逐身份语言向量。 |
| 候选位点和身份 | `result/predictions/match.json` 包含区域世界 XYZ 中心、前景概率和候选配体身份得分。 |
| 对接选择和姿态 | `result/predictions/plan.json` 保存自动或手工决定，`build.json` 索引实际提交。每个成功请求的姿态保存在 `build/<blob_id>/poses.sdf`。 |

样本字段和三种受体来源见 [Prepare 使用说明](Prepare/README.md)及[输入与产物契约](Prepare/INPUTS.md)。分阶段运行、人工审阅和对接预算见 [Infer 使用说明](Infer/README.md)及[选择字段契约](Infer/SELECTION.md)。软件、权重和通用词典在安装时准备一次，后续样本直接使用本机安装。

自有代码采用 [Apache-2.0](LICENSE)，自有模型权重采用 [CC BY 4.0](WEIGHTS_LICENSE.md)。随附及单独安装的第三方组件保留各自的许可。
