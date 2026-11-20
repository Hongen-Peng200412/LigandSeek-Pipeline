# 第三方运行组件与许可

项目自有发布代码采用 Apache-2.0，自有的 Find、Match、Build 训练权重采用 CC BY 4.0。对第三方实现、工具和数据，应按其来源使用相应许可。特别是完整运行 Find 所用的 PTv3 实现和 CryoAtom2 路径所用的 Chimera 各有自己的使用条件。模型来源、版本和运行时取得方式如下。

| 组件 | 本版本实际使用的内容 | 来源和许可 |
| --- | --- | --- |
| CryoAtom2 2.1.1 | 固定 Git 提交 `856e250df7b784b854b892f1b619d32d51188cef` 的源码，官方 v2.1 的三份权重。按原流程取得 ESM-2 650M 与 RNA-FM 所需权重。 | [Yang Lab 源码与原许可](https://github.com/YangLab-SDU/CryoAtom/tree/856e250df7b784b854b892f1b619d32d51188cef)、[ESM](https://github.com/facebookresearch/esm)、[RNA-FM](https://github.com/ml4bio/RNA-FM)。这些权重从各自官方服务获取。 |
| SMI-TED Light 289M | IBM 修订版 `414c3ea0a8603ef49d1c5bb3db336e09877c01ce` 提供的 encoder、词表和 `smi-ted-Light_40.pt`。 | [IBM 模型卡与 Apache-2.0 许可](https://huggingface.co/ibm-research/materials.smi-ted/tree/414c3ea0a8603ef49d1c5bb3db336e09877c01ce) |
| UCSF Chimera Classic 1.19 | Linux x86-64 OSMesa build 42556。安装时从 UCSF 取得，或使用用户已获许可的安装。 | [UCSF 下载及非商业许可](https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64_osmesa/chimera-1.19-linux_x86_64_osmesa.bin) |
| Point Transformer V3 | Find 所用的改写网络与点序列化实现。对应源码随科学库发布。 | [上游模型和 CC BY-NC 4.0 许可](https://github.com/Pointcept/PointTransformerV3/blob/main/LICENSE) |
| PocketXMol、PoseBusters、GVP | DG-PocketXMol、距离碰撞与分子身份评分，以及依赖的几何网络实现。来源与随包许可保存在科学库归档。 | [PocketXMol](https://github.com/pengxingang/PocketXMol)、[PoseBusters](https://github.com/maabuu/posebusters)、[GVP](https://github.com/drorlab/gvp-pytorch) |
| Emap2lig 化学模板相关代码 | AdaLigand 化学准备中沿用的配体化学模板结构。修改和引用见科学库内随附说明。 | [Emap2lig 与 GPL-3.0](https://github.com/KiharaLab/Emap2lig) |
| FlashAttention 2.8.3 | `native.tar.gz` 内原生扩展，匹配 PyTorch 2.4.1+cu121 与 CPython 3.10。附原 BSD 3-Clause 许可。 | [上游源码](https://github.com/Dao-AILab/flash-attention/tree/v2.8.3)、[随包许可](notices/flash-attention.BSD-3-Clause.txt) |
| pytorch-fast-transformers 0.4.0 | `native.tar.gz` 内 Idiap 的 CPU 构建，与 SMI-TED 使用的算子一致。附原 MIT 许可。 | [上游源码](https://github.com/idiap/fast-transformers/tree/v0.4.0)、[随包许可](notices/fast-transformers.MIT.txt) |
| wwPDB Chemical Component Dictionary | `ccd.tar.gz` 中的通用模板及 `components.cif.gz` 完整词典。 | [wwPDB 词典](https://www.wwpdb.org/data/ccd)、[wwPDB 数据使用条款](https://www.wwpdb.org/about/usage-policies) |

Conda 和 pip 依赖从各自发布渠道下载。每个二进制包保留自身携带的许可与署名。科学库归档中各源文件对应的作者、来源及修改范围见[科学库文件来源](notices/runtime-sources.md)。
