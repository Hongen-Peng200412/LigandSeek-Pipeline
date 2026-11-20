# 安装 Find/Match、化学输入、Build 和 SMI-TED 所需的解释器. 由仓库根 install.sh 读取.
# 锁文件逐包指定 Linux x86-64 构建, 使 RDKit 可在没有系统编译器的环境安装.
conda_prefix() {
    local role=$1 target="$prefix/envs/$1" manager operation=create
    local lock="$root/install/locks/$1.txt"
    if [[ -f "$target/conda-meta/history" ]]; then
        if cmp -s "$lock" "$target/.ligandseek-conda.txt"; then return; fi
        operation=install
    fi
    if command -v micromamba >/dev/null; then manager=$(command -v micromamba)
    elif [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then manager=$CONDA_EXE
    elif command -v conda >/dev/null; then manager=conda
    else
        manager="$prefix/bin/micromamba"
        if [[ ! -x "$manager" ]]; then
            fetch 'https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.10.0-0.tar.bz2' "$scratch/micromamba.tar.bz2"
            printf '%s  %s\n' 9a97976244929bbae28fb1fae504e74a3e7e6d1932870e3a03321f01b7fbd07b "$scratch/micromamba.tar.bz2" | sha256sum -c -
            tar -xjf "$scratch/micromamba.tar.bz2" -C "$scratch" bin/micromamba
            mv "$scratch/bin/micromamba" "$manager"
        fi
    fi
    MAMBA_ROOT_PREFIX="$prefix/mamba" "$manager" "$operation" -y -p "$target" --file "$lock"
    cp "$lock" "$target/.ligandseek-conda.txt"
}

# 显式路径优先, 其次检查本安装目录或 PATH 中的 Python. 工具包名和版本来自上线时使用的运行环境.
choose_python() {
    local supplied=$1 role=$2 probe=$3 candidate
    if [[ -n "$supplied" ]]; then
        "$supplied" -c "$probe" >&2 || return 1
        realpath -s "$(command -v "$supplied")"
        return
    fi
    for candidate in "$prefix/envs/$role/bin/python" "$(command -v python3 || true)"; do
        if [[ -x "$candidate" ]] && "$candidate" -c "$probe" >/dev/null 2>&1; then
            printf '%s\n' "$candidate"
            return
        fi
    done
}
find_probe='import torch, torchvision, flash_attn, torch_cluster, torch_scatter, lightning, timm, addict, einops, mrcfile; assert torch.__version__=="2.4.1+cu121" and flash_attn.__version__=="2.8.3"'
chem_probe='import gemmi, scipy, numpy, mrcfile, pdbeccdutils, requests; from rdkit import rdBase; assert rdBase.rdkitVersion=="2026.03.3" and gemmi.__version__=="0.7.5"'
build_probe='import torch, torch_cluster, torch_scatter, torch_sparse, torch_geometric, lmdb, PeptideBuilder, easydict, pytorch_lightning, pandas; from rdkit import rdBase; assert torch.__version__=="2.6.0+cu124" and rdBase.rdkitVersion=="2023.09.3"'
smi_probe='import torch, transformers, fast_transformers, pandas, regex; from fast_transformers.attention import AttentionLayer; from rdkit import rdBase; assert torch.__version__=="2.4.1+cu121" and transformers.__version__=="4.57.6" and rdBase.rdkitVersion=="2026.03.3"'
find_python=$(choose_python "$find_python" find "$find_probe")
chem_python=$(choose_python "$chem_python" chemistry "$chem_probe")
build_python=$(choose_python "$build_python" build "$build_probe")
smi_python=$(choose_python "$smi_python" smi "$smi_probe")
if [[ -z "$find_python" ]]; then
    conda_prefix find
    find_python="$prefix/envs/find/bin/python"
    directory_asset native "$prefix"
    "$find_python" -m pip install --index-url https://download.pytorch.org/whl/cu121 torch==2.4.1 torchvision==0.19.1
    "$find_python" -m pip install --only-binary=:all: --no-index --no-deps -f 'https://data.pyg.org/whl/torch-2.4.0+cu121.html' torch-scatter==2.1.2+pt24cu121 torch-cluster==1.6.3+pt24cu121
    "$find_python" -m pip install lightning==2.2.5 pytorch-lightning==2.2.5 torchmetrics==1.3.2 timm==0.9.16 addict==2.4.0 einops==0.8.2
    "$find_python" -m pip install --no-deps "$prefix"/native/flash_attn-2.8.3-*.whl
    "$find_python" -c "$find_probe"
fi
if [[ -z "$chem_python" ]]; then
    conda_prefix chemistry
    chem_python="$prefix/envs/chemistry/bin/python"
    "$chem_python" -m pip install mrcfile==1.5.4
    "$chem_python" -c "$chem_probe"
fi
if [[ -z "$build_python" ]]; then
    conda_prefix build
    build_python="$prefix/envs/build/bin/python"
    "$build_python" -m pip install --index-url https://download.pytorch.org/whl/cu124 torch==2.6.0
    "$build_python" -m pip install --only-binary=:all: --no-index --no-deps -f 'https://data.pyg.org/whl/torch-2.6.0+cu124.html' torch-scatter==2.1.2+pt26cu124 torch-cluster==1.6.3+pt26cu124 torch-sparse==0.6.18+pt26cu124
    "$build_python" -m pip install torch-geometric==2.6.1 lmdb==1.7.5 PeptideBuilder==1.1.0 einops==0.8.1 pytorch-lightning==2.6.0 torchmetrics==1.9.0
    "$build_python" -c "$build_probe"
fi
if [[ -z "$smi_python" ]]; then
    directory_asset native "$prefix"
    smi_python="$prefix/envs/smi/bin/python"
    "$find_python" -m venv --system-site-packages "$prefix/envs/smi"
    "$smi_python" -m pip install transformers==4.57.6 huggingface-hub==0.36.0 tokenizers==0.22.2
    "$smi_python" -m pip install --no-deps "$prefix"/native/pytorch_fast_transformers-0.4.0-*.whl
    # 两个解释器均为 CPython 3.10. 只连接化学包, SMI-TED 的 numpy、pandas 与 PyTorch 来自 Find 环境.
    chem_site=$("$chem_python" -c 'from pathlib import Path; import rdkit; print(Path(rdkit.__file__).parent.parent)')
    smi_site=$("$smi_python" -c 'import sysconfig; print(sysconfig.get_path("purelib"))')
    for pkg in "$chem_site/rdkit" "$chem_site"/rdkit-*.dist-info; do
        [[ -e "$smi_site/${pkg##*/}" ]] || ln -s "$pkg" "$smi_site/"
    done
    "$smi_python" -c "$smi_probe"
fi
