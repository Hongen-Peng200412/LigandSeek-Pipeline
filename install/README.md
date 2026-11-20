# 在 Linux GPU 服务器上安装

在 Linux x86-64 服务器克隆同一个发布标签，首次安装后用一条命令处理一个样本。Find 使用 CUDA 12.1 的 PyTorch，Build 使用 CUDA 12.4 的 PyTorch，带序列与无序列的 CryoAtom2 均使用其原有 CUDA 11.8 环境。需要能运行 CUDA 12.4 应用程序的 NVIDIA 驱动，包内的原生扩展以 CPython 3.10 和 H200 为验证目标。

```bash
git clone --branch v1.0.0 https://github.com/Hongen-Peng200412/LigandSeek-Pipeline.git
cd LigandSeek-Pipeline
bash install.sh --accept-chimera-noncommercial-license
source "$HOME/.local/ligandseek/environment.sh"
ligandseek run raw.json --output result
```

上述首次安装命令适用于已阅读并接受 UCSF 非商业许可的用户。已有许可安装时使用 `bash install.sh --chimera /opt/UCSF-Chimera64-1.19/bin/chimera`。安装器直接向 UCSF 取得无图形界面的 Chimera Classic 1.19。用户可在一台能够访问 GitHub、conda-forge、PyTorch、PyG、Hugging Face、Yang Lab、Meta 和 UCSF 下载服务的机器上完成安装。系统命令需要 Bash、curl、tar、gzip、bzip2、unzip、flock（util-linux）和 GNU coreutils。首次安装会下载多个 GPU 运行库和模型，安装区和下载缓存都需要充足空间。[UCSF 许可与下载](https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?file=linux_x86_64_osmesa/chimera-1.19-linux_x86_64_osmesa.bin)

## 使用已有模型或解释器

Hugging Face 访问受限时，可在安装命令前设置标准环境变量 `HF_ENDPOINT`，例如 `HF_ENDPOINT=https://hf-mirror.com bash install.sh --accept-chimera-noncommercial-license`。安装器仍固定 IBM 修订号并核验权重、词表、编码脚本三个 SHA-256。

`--prefix` 指定本次安装和缓存的位置，默认 `$HOME/.local/ligandseek`。在另外的终端使用这一套安装时，载入该目录的 `environment.sh`。可以混用本地依赖和自动安装，例如下面的安装复用 Find 的权重与现有 Chimera，其余模型和工具从公开来源补齐。

```bash
bash install.sh --prefix "$HOME/opt/ligandseek" --find /data/models/find --chimera /opt/chimera/bin/chimera
source "$HOME/opt/ligandseek/environment.sh"
```

| 可复用的内容 | 安装选项 | 给定路径中的实际文件 |
| --- | --- | --- |
| Find | `--find`、`--find-python` | Find 权重、配置、两种受体的冻结校准目录，以及可选的 Find/Match 解释器，例如 `/opt/find/bin/python`。 |
| Match | `--match` | 含 `BEST.pt` 的目录。Match 与 Find 共用 PyTorch 环境。 |
| Build | `--build`、`--build-python` | 含 `train.yml` 与 `step=24800.ckpt` 的目录，以及可选的 DG-PocketXMol 解释器。 |
| CryoAtom2 | `--cryoatom-python` | 已安装 2.1.1 的解释器或启动脚本，能够读取三份结构权重，以及 ESM-2 和 RNA-FM 的语言权重。 |
| UCSF Chimera | `--chimera` | Classic 1.19 build 42556 可执行文件，例如 `/opt/chimera/bin/chimera`。 |
| SMI-TED | `--smi-ted-model`、`--smi-ted-python` | IBM 模型目录中的 `smi-ted-Light_40.pt`、`bert_vocab_curated.txt` 和 `smi-ted/inference/smi_ted_light/load.py`。独立解释器须能导入 RDKit 2026.03.3、PyTorch 2.4.1+cu121、fast-transformers 0.4.0 和 transformers 4.57.6。 |
| 实验图与化学结构处理 | `--chemistry-python` | 安装 Gemmi、mrcfile、RDKit 和 pdbeccdutils 的化学解释器，例如 `/opt/chemistry/bin/python`。 |
| 已展开的完整模型包 | `--models` | 公用父目录，下含 `find/`、`match/`、`build/`、`libraries/` 和 `ccd/`。单独给出的阶段目录优先。 |
| 已下载的发布附件 | `--assets` | 同一发布的附件文件夹，含 `SHA256SUMS` 及本机尚缺的对应归档，例如 `find.tar.gz`。每个实际读取的附件会校验 SHA-256。 |

目录参数可以是绝对路径或相对当前目录的路径。工具解释器优先采用显式指定值，其次检查本前缀中的适配环境和 PATH 中的兼容解释器，缺失时按锁定构建安装。只有从安装器获取的模型目录版式才能直接作为 `--find`、`--match` 或 `--build` 参数。将本机其他工程的输出整理好以后，可在同一次安装中指定其中任意几个目录。各环境执行的版本检测见 [environments.sh](environments.sh)。

当本机已具备全部工具和资产时，可以离线安装。传入含完整词典的 `--models`、五个 `--*-python`、`--chimera`、`--smi-ted-model`。安装器用指定解释器安装公开包到自己的 `packages/` 目录。若仅提前取得了发布附件，创建 Python 环境以及下载外部模型仍需访问各官方源。

## 安装目录

两个目录树均对应 `--prefix` 指向的同一位置。表中目录名与发布附件的顶层目录一致。

```text
## <科学产物>
$HOME/.local/ligandseek/
├── models/
│   ├── find/           # Find 权重与两种受体的分类、区域选择校准
│   ├── match/          # 候选身份预测权重
│   ├── build/          # 姿态生成权重及配置
│   ├── libraries/      # 四个领域库本次推理需要的代码与配置
│   └── ccd/            # 通用化学组分模板和固定版本 components.cif.gz
└── tools/              # 安装器获取的 CryoAtom2 与 SMI-TED 权重、Chimera 程序和 torch 网络缓存

## <其他文件>
$HOME/.local/ligandseek/
├── bin/                # ligandseek 命令及本机工具入口（可再生辅助导出）
├── environment.sh      # 本机解释器和资产路径（可再生辅助导出）
├── packages/           # 已安装的 Prepare、Infer（可再生辅助导出）
├── envs/               # 自动创建的各解释器环境（可再生辅助导出）
├── native/             # 与 CUDA 和 Python ABI 配套的两个扩展 wheel（可再生工具输入）
└── downloads/          # Conda、pip 及固定发布附件缓存（可再生工具输入）
```

### `environment.sh`

加载它会把 `bin/` 加入当前 Bash 的 PATH，为后续子进程提供本次安装的包路径、各解释器、Chimera 和权重路径。例如 `LIGANDSEEK_MODELS=$HOME/.local/ligandseek/models` 指向科学资产根；`LIGANDSEEK_CA2_PYTHON` 指向 CryoAtom2 专用解释器或带独立模型缓存的启动脚本。`$prefix/bin/ligandseek` 也会自行加载该配置，可以从其他工作目录直接调用。

### `models/ccd/components.cif.gz`

该发布固定的 wwPDB 化学组分词典。每个 CIF 数据块定义一个组分的原子名、元素和化学键，CCD 三字母标准组分例如 `ALA`。程序先复用已有模板，缺项从这份词典本地转换。提供的 CCD 缓存会复制到可写安装区，方便处理新样本。[wwPDB CCD](https://www.wwpdb.org/data/ccd)

### `locks/*.txt` 与第三方来源

`install/locks/{find,chemistry,build,cryoatom}.txt` 是 Linux x86-64 的 Conda 显式文件，每条记录给出一个精确包构建的公共下载地址和 MD5。安装器随后补入对应的 PyTorch/PyG wheel。SMI-TED 在 Find 解释器之上建立单独的 Python 环境，连接同一份 RDKit 化学库。CryoAtom2、IBM 和 Chimera 的原始下载步骤见 [external.sh](external.sh)，第三方条款见 [第三方组件说明](THIRD_PARTY.md)。
