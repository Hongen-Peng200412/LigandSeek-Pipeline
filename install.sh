#!/usr/bin/env bash
# 从同一份源码和固定附件安装公开的 Prepare + Infer, 用法见 install/README.md.
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
prefix=${LIGANDSEEK_PREFIX:-"$HOME/.local/ligandseek"}
assets= find_dir= match_dir= build_dir= models_root=${LIGANDSEEK_MODELS:-}
find_python=${LIGANDSEEK_FIND_PYTHON:-}
chem_python=${LIGANDSEEK_PREPARE_PYTHON:-}
build_python=${LIGANDSEEK_BUILD_PYTHON:-}
ca2_python=${LIGANDSEEK_CA2_PYTHON:-}
smi_python=${LIGANDSEEK_SMITED_PYTHON:-}
smi_model=${LIGANDSEEK_SMITED_MODEL:-}
chimera=${LIGANDSEEK_CHIMERA:-}
accept_license=0
while (($#)); do
    case "$1" in
        --prefix) prefix=$2 ;;
        --assets) assets=$2 ;;
        --models) models_root=$2 ;;
        --find) find_dir=$2 ;;
        --match) match_dir=$2 ;;
        --build) build_dir=$2 ;;
        --find-python) find_python=$2 ;;
        --chemistry-python) chem_python=$2 ;;
        --build-python) build_python=$2 ;;
        --cryoatom-python) ca2_python=$2 ;;
        --smi-ted-python) smi_python=$2 ;;
        --smi-ted-model) smi_model=$2 ;;
        --chimera) chimera=$2 ;;
        --accept-chimera-noncommercial-license) accept_license=1; shift; continue ;;
        -h|--help) cat "$root/install/README.md"; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
    shift 2
done
mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)
exec 9>"$prefix/.install.lock"
flock -x 9
release=v1.0.0
base_url="https://github.com/Hongen-Peng200412/LigandSeek-Pipeline/releases/download/$release"
cache="$prefix/downloads/$release"
models="$prefix/models"
mkdir -p "$cache" "$models" "$prefix/envs" "$prefix/bin"
export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-"$prefix/downloads/conda"}
export PIP_CACHE_DIR=${PIP_CACHE_DIR:-"$prefix/downloads/pip"}
scratch=$(mktemp -d "$prefix/.install.XXXXXXXX")
trap 'rm -rf -- "$scratch"' EXIT

# 下载、核验与解包只发生在安装阶段. 每个归档内的首级目录与安装目录名相同.
fetch() {
    local url=$1 destination=$2
    curl --fail --show-error --silent --location --connect-timeout 20 --speed-time 120 --speed-limit 512 \
        --retry 3 --output "$scratch/download" "$url"
    mv -f -- "$scratch/download" "$destination"
}
asset() {
    local name=$1 checksum
    if [[ -z "${inventory:-}" ]]; then
        if [[ -n "$assets" ]]; then inventory="$(cd "$assets" && pwd)/SHA256SUMS"
        else inventory="$scratch/SHA256SUMS"; fetch "$base_url/SHA256SUMS" "$inventory"; fi
    fi
    checksum=$(awk -v f="$name" '$2==f || $2=="*"f {print $1}' "$inventory")
    [[ $checksum =~ ^[0-9a-fA-F]{64}$ ]] || { echo "Missing release checksum: $name" >&2; return 1; }
    if [[ -n "$assets" ]]; then archive="$(dirname "$inventory")/$name"
    else
        archive="$cache/$name"
        if [[ ! -f "$archive" ]] || ! printf '%s  %s\n' "$checksum" "$archive" | sha256sum --status -c -; then
            fetch "$base_url/$name" "$archive"
        fi
    fi
    printf '%s  %s\n' "$checksum" "$archive" | sha256sum --status -c -
}
directory_asset() {
    local role=$1 dest=$2
    [[ -d "$dest/$role" ]] && return 0
    asset "$role.tar.gz"
    tar -xzf "$archive" -C "$scratch"
    mv -T -- "$scratch/$role" "$dest/$role"
}
source "$root/install/environments.sh"
source "$root/install/external.sh"
for role in find match build libraries ccd; do
    variable="${role}_dir"
    source_dir=${!variable:-}
    if [[ -z "$source_dir" && -n "$models_root" ]]; then source_dir="$models_root/$role"; fi
    if [[ -n "$source_dir" ]]; then
        source_dir=$(cd "$source_dir" && pwd -P)
        if [[ "$role" == ccd && "$source_dir" != "$models/ccd" ]]; then
            if [[ ! -d "$models/ccd" ]]; then
                cp -rL -- "$source_dir" "$scratch/ccd"
                mv -T -- "$scratch/ccd" "$models/ccd"
            fi
        elif [[ -d "$models/$role" ]]; then [[ $(cd "$models/$role" && pwd -P) == "$source_dir" ]]
        else ln -s "$source_dir" "$models/$role"; fi
    else directory_asset "$role" "$models"; fi
done
if [[ ! -s "$models/ccd/components.cif.gz" ]]; then
    asset components.cif.gz
    cp "$archive" "$scratch/components.cif.gz"
    mv "$scratch/components.cif.gz" "$models/ccd/"
fi

# 公开的两个纯 Python 包集中放在一个目录, 各解释器可以独立复用, 无须更改用户已有环境.
cp -R -- "$root/Infer" "$root/Prepare" "$scratch/"
"$find_python" -m pip install --no-index --no-deps --no-build-isolation --upgrade \
    --target "$prefix/packages" "$scratch/Infer" "$scratch/Prepare"
configuration="$prefix/environment.sh"
{
    printf 'export LIGANDSEEK_MODELS=%q\n' "$models"
    printf 'export LIGANDSEEK_FIND_PYTHON=%q\n' "$find_python"
    printf 'export LIGANDSEEK_PREPARE_PYTHON=%q\n' "$chem_python"
    printf 'export LIGANDSEEK_BUILD_PYTHON=%q\n' "$build_python"
    printf 'export LIGANDSEEK_CA2_PYTHON=%q\n' "$ca2_python"
    printf 'export LIGANDSEEK_SMITED_PYTHON=%q\n' "$smi_python"
    printf 'export LIGANDSEEK_SMITED_MODEL=%q\n' "$smi_model"
    printf 'export LIGANDSEEK_CHIMERA=%q\n' "$chimera"
    printf 'export PYTHONPATH=%q\n' "$prefix/packages"
    printf 'export PYTHONNOUSERSITE=1\n'
    printf 'export PATH=%q:$PATH\n' "$prefix/bin"
} > "$configuration"
printf '#!/usr/bin/env bash\nsource %q\nexec %q -m ligandseek "$@"\n' \
    "$configuration" "$find_python" > "$prefix/bin/ligandseek"
chmod +x "$prefix/bin/ligandseek"
printf 'Installed. To use this installation: source %q\n' "$configuration"
