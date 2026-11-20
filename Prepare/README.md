# 从原始冷冻电镜图准备 LigandSeek 输入

Prepare 为一个样本生成受体、共用实验采样网格的两张密度图和候选配体语言向量。完成后可把整个输出目录移动到另一台机器，由 Infer 读取其中的 `sample.json`。

## 准备一个样本

例如，把下面的构造示例写成 `raw.json`。相对输入路径以此文件所在目录为基准。

```json
{
  "sample_id": "demo",
  "experimental_map": "maps/demo.mrc.gz",
  "protein_fasta": "sequences/protein.fa",
  "smiles": ["OCC", "CCO", "CCN"],
  "resolution": 3.2
}
```

在已完成模型与工具安装的环境中运行以下命令。第一条产生准备产物，第二条读取这些产物完成 Find、Match 和 Build。一次运行全部阶段也可以直接使用公开的 `ligandseek run` 命令。

```sh
ligandseek prepare raw.json --output prepared
ligandseek infer prepared/sample.json --output predictions
```

如已设置 `LIGANDSEEK_MODELS`，Prepare 可省略 `--models`。`--output` 指定新目录，只在本次准备成功后发布。Python 入口是 `prepare_sample(Path("raw.json"), Path("models"), Path("prepared"))`，可从 `ligandseek_prepare` 导入。

| 原始输入字段 | 用户填写的内容 | 构造示例 |
| --- | --- | --- |
| `experimental_map` | 必填。原始实空间 MRC/map，也支持 gzip 压缩。 | `"maps/demo.mrc.gz"` |
| `smiles` | 必填。可填单个 SMILES 字符串，或依次列出各候选的字符串。空列表对应空身份轴。 | `"CCO"` 或 `["OCC", "CCO", "CCN"]` |
| `sample_id` | 样本标识，缺省取实验图文件名首个句点之前的部分，转为小写。首字符是字母或数字，其余可含字母、数字、句点、下划线和短横线。 | `"demo"` |
| `protein_fasta`、`rna_fasta`、`dna_fasta` | 各自对应一份 FASTA，按实际存在的序列填写。全部省略且没有受体时，CryoAtom2 走官方无序列模式。 | `"sequences/protein.fa"` |
| `receptor_cif` | 已有且与实验图处于同一坐标系的未压缩受体 mmCIF。给出时直接使用，即使同时提供了 FASTA。 | `"structures/receptor.cif"` |
| `receptor_source` | 直接提供 CIF 时标明受体条件。缺省为已知受体 `deposited`；先前由 CryoAtom2 生成的 CIF 填 `ca2`。Infer 据此选择对应的 Find 校准。 | `"ca2"` |
| `resolution` | 模拟受体密度时采用的正分辨率，单位 Å。 | `3.2` |
| `emdb_metadata`、`rcsb_metadata` | 同一实验的 EMDB 或 RCSB 条目元数据 JSON。用户给出分辨率时优先使用手填值。 | `"maps/emdb.json"` |
| `map_convention` | 缺省 `mrc`，按 MRC 首采样点解释世界坐标。使用 AdaLigand 历史网格约定的样本可指定 `adaligand`。 | `"mrc"` |

分辨率依次取手填值、同图元数据中 EMDB 最终重建或 RCSB 报告的分辨率、MRC 标签中带 Å 单位的 `resolution`，最后缺省为 3.0 Å。体素间距单独保存在输出几何中。若原图只给出 `resolution=3.2 A` 标签，`prepare.json` 会记录数值 `3.2` 与来源 `mrc_label`。

## 阅读代码与产物

先看 `ligandseek_prepare/pipeline.py` 的 `prepare_sample` 与 `assemble_sample`。前者调用化学环境，后者按受体、密度、配体向量的顺序生成样本。然后可以按需要阅读下表中的一个模块。

| 模块 | 入口和科学计算 |
| --- | --- |
| `receptor.py` | `prepare_receptor` 选择已有受体，或按实际存在的序列调用 CryoAtom2 的 `build`。使用官方最终过滤的 CIF。 |
| `density.py` | `prepare_density` 用原 Fourier 网格变换及 Chimera Classic 的 `molmap … onGrid` 生成受体模拟密度。采用首个模型、原替代构象选择规则和 `ATOM` 重原子。 |
| `language.py` | `encode_records` 将每个规范身份交给 SMI-TED 的独立解释器，`encode_batch` 使用原 `run_smi_ted` 记录编码状态、诊断和 768 维向量。 |

输出文件、数组形状、身份编号和坐标示例见[单样本产物契约](INPUTS.md)。

运行开始时会打印与目标相邻的 `<output>.unfinished-*` 工作目录。若某个外部程序失败，该目录保留已完成的阶段文件、工具日志及尚未结束的工具工作目录。全部完成后将其改名为请求的输出目录。

## 复用本机工具

安装器管理各解释器和模型资产。已有工具可通过 `LIGANDSEEK_CA2_PYTHON`、`LIGANDSEEK_PREPARE_PYTHON`、`LIGANDSEEK_SMITED_PYTHON` 分别提供 CryoAtom2、化学准备和 SMI-TED 解释器。例如可将第一项指到 CryoAtom2 2.1.1 环境的 `bin/python`。Chemistry 环境应可导入已安装的 `ligandseek`。

`LIGANDSEEK_CHIMERA` 指向 Chimera Classic 1.19 的可执行文件，缺省从 PATH 查找 `chimera`。`LIGANDSEEK_SMITED_MODEL` 指向同一目录下包含 `smi-ted-Light_40.pt`、`bert_vocab_curated.txt` 和 `smi-ted/inference/smi_ted_light/load.py` 的官方模型副本，缺省为模型根下的 `smi_ted_289m/`。三个 Python 路径分别缺省为启动相应步骤的解释器。
