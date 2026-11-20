# 官方工具和预训练参数. 本文件由 install.sh 读取, 每次公开运行使用安装时记录的本机路径.
verified_download() {
    local url=$1 destination=$2 sha=$3
    if [[ -s "$destination" ]] && printf '%s  %s\n' "$sha" "$destination" | sha256sum --status -c -; then return; fi
    mkdir -p "$(dirname "$destination")"
    fetch "$url" "$scratch/candidate"
    printf '%s  %s\n' "$sha" "$scratch/candidate" | sha256sum -c -
    mv "$scratch/candidate" "$destination"
}

if [[ -z "$chimera" ]]; then
    for candidate in "$prefix/tools/chimera/bin/chimera" "$(command -v chimera || true)"; do
        if [[ -x "$candidate" ]] && "$candidate" --version 2>&1 | grep '1.19 (build 42556)' >/dev/null; then chimera=$candidate; break; fi
    done
fi
if [[ -z "$chimera" ]]; then
    if (( ! accept_license )); then
        echo 'Provide --chimera with a licensed Chimera Classic 1.19 executable or accept the UCSF non-commercial license with --accept-chimera-noncommercial-license.' >&2
        exit 2
    fi
    ucsf='https://www.cgl.ucsf.edu'
    binary="$cache/chimera-1.19-linux_x86_64_osmesa.bin"
    md5=aace74cbbbb3dfc5acc7d88b2851ae44
    if ! printf '%s  %s\n' "$md5" "$binary" | md5sum --status -c - 2>/dev/null; then
        curl --fail --show-error --silent --retry 3 --cookie-jar "$scratch/cookies" \
            --data-urlencode 'file=linux_x86_64_osmesa/chimera-1.19-linux_x86_64_osmesa.bin' --data-urlencode 'choice=Accept' \
            -o "$scratch/notice" "$ucsf/chimera/cgi-bin/secure/chimera-get.py"
        link=$(sed -n 's/.*content="[0-9][0-9]*;url=\([^"]*\)".*/\1/p' "$scratch/notice" | head -n 1 | sed 's/&amp;/\&/g')
        case "$link" in https://*) ;; /*) link="$ucsf$link" ;; *) echo 'The UCSF license form supplied no download URL.' >&2; exit 1 ;; esac
        curl --fail --show-error --silent --location --retry 3 --cookie "$scratch/cookies" -o "$scratch/chimera.bin" "$link"
        [[ $(wc -c < "$scratch/chimera.bin") -eq 149448194 ]]
        printf '%s  %s\n' "$md5" "$scratch/chimera.bin" | md5sum -c -
        mv "$scratch/chimera.bin" "$binary"
    fi
    mkdir -p "$prefix/tools"
    chmod u+x "$binary"
    printf '%s\n0\n\n' "$prefix/tools/chimera" | "$binary"
    chimera="$prefix/tools/chimera/bin/chimera"
fi
"$chimera" --version 2>&1 | grep '1.19 (build 42556)'
chimera=$(realpath -s "$(command -v "$chimera")")

provided_smi=${smi_model:+yes}
smi_model=${smi_model:-"$prefix/tools/smi-ted"}
hf="${HF_ENDPOINT:-https://huggingface.co}/ibm-research/materials.smi-ted/resolve/414c3ea0a8603ef49d1c5bb3db336e09877c01ce"
while read -r sha file; do
    if [[ -z "$provided_smi" ]]; then verified_download "$hf/$file" "$smi_model/$file" "$sha"
    else (cd "$smi_model" && printf '%s  %s\n' "$sha" "$file" | sha256sum -c -); fi
done < "$root/install/locks/smi-ted.sha256"
smi_model=$(cd "$smi_model" && pwd)

ca2_probe='import CryoAtom2, torch; from pathlib import Path; assert CryoAtom2.__version__=="2.1.1"; assert all((Path(CryoAtom2.__file__).parent/"checkpoint"/n).is_file() for n in ("CryoNet.pth","CryoNet_no_seq.pth","RUNet.pth")); assert all((Path(torch.hub.get_dir())/"checkpoints"/n).is_file() for n in ("RNA-FM_pretrained.pth","esm2_t33_650M_UR50D.pt","esm2_t33_650M_UR50D-contact-regression.pt"))'
if [[ -z "$ca2_python" && -x "$prefix/bin/cryoatom-python" ]]; then ca2_python="$prefix/bin/cryoatom-python"; fi
ca2_python=$(choose_python "$ca2_python" cryoatom "$ca2_probe")
if [[ -z "$ca2_python" ]]; then
    conda_prefix cryoatom
    ca2_python="$prefix/envs/cryoatom/bin/python"
    "$ca2_python" -m pip install einops==0.6.0 scipy==1.7.3 biopython==1.81 matplotlib==3.5.3 mrcfile==1.4.3 \
        pandas==1.3.5 fair-esm==2.0.0 pyhmmer==0.7.1 joblib==1.0.1 ptflops==0.6.6 scikit-learn==0.24.0 \
        threadpoolctl==2.2.0 rna-fm==0.2.2 openpyxl==3.1.5
    revision=856e250df7b784b854b892f1b619d32d51188cef
    ca2_dir="$prefix/tools/CryoAtom-$revision"
    mkdir -p "$prefix/tools"
    if [[ ! -d "$ca2_dir" ]]; then
        fetch "https://codeload.github.com/YangLab-SDU/CryoAtom/tar.gz/$revision" "$scratch/source.tar.gz"
        tar xzf "$scratch/source.tar.gz" -C "$scratch"
        mv -T "$scratch/CryoAtom-$revision" "$ca2_dir"
    fi
    yanglab='https://yanglab.qd.sdu.edu.cn/CryoAtom/download'
    checkpoint="$ca2_dir/CryoAtom2/checkpoint"
    mkdir -p "$checkpoint"
    if ! (cd "$checkpoint" && sha256sum --status -c "$root/install/locks/cryoatom-weights.sha256"); then
        fetch "$yanglab/checkpoints_v2.1.zip" "$scratch/checkpoints.zip"
        unzip -joq "$scratch/checkpoints.zip" '*CryoNet.pth' '*CryoNet_no_seq.pth' '*RUNet.pth' -d "$checkpoint"
        (cd "$checkpoint" && sha256sum -c "$root/install/locks/cryoatom-weights.sha256")
    fi
    chmod +x "$ca2_dir/CryoAtom2/utils/getp"
    "$ca2_python" -m pip install --no-deps --no-build-isolation --editable "$ca2_dir"
    hub="$prefix/tools/torch/hub/checkpoints"
    verified_download "$yanglab/RNA-FM_pretrained.pth" "$hub/RNA-FM_pretrained.pth" 5b5d7d87b37c291ef42c140ef9edf7aea29f255fa2a4fd435f776c52e93d5e99
    verified_download 'https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt' "$hub/esm2_t33_650M_UR50D.pt" ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c
    verified_download 'https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt' "$hub/esm2_t33_650M_UR50D-contact-regression.pt" 8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575
    printf '#!/usr/bin/env bash\nexport TORCH_HOME=%q\nexec %q "$@"\n' "$prefix/tools/torch" "$ca2_python" > "$prefix/bin/cryoatom-python"
    chmod +x "$prefix/bin/cryoatom-python"
    ca2_python="$prefix/bin/cryoatom-python"
    "$ca2_python" -c "$ca2_probe"
fi
