"""确定本次受体. prepare_receptor 写出 output/receptor.cif, 返回给 prepare.json 的受体来源、序列种类和实际命令."""

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def prepare_receptor(raw: dict, experiment: Path, base: Path, output: Path) -> dict:
    """沿用用户提供的受体或调用 CryoAtom2 生成最终受体.

    - raw: 单样本输入. receptor_cif 优先, receptor_source 可注明已生成的 ca2 受体. protein_fasta、rna_fasta、dna_fasta 分别指定三种序列文件.
    - experiment: 原始实验图 MRC/map 或 gzip 路径, 处于 receptor_cif 的世界 XYZ 坐标系.
    - base: raw.json 所在目录, 用于解析相对输入路径.
    - output: 写入 receptor.cif. 调用 CryoAtom2 时还保存 cryoatom2.log, 包含完整 stdout 和 stderr.
    - 返回值: source 为 ca2 或 deposited, mode 为 supplied、with_sequences 或 without_sequences, sequences 为种类到 FASTA 路径的映射, command 为本次实际参数列表或 None, version 为实际 CryoAtom2 版本或 None, input_cif 是提供的受体路径或 None.
    """
    if raw.get("receptor_cif"):
        cif = (base / raw["receptor_cif"]).resolve()
        shutil.copyfile(cif, output / "receptor.cif")
        return {"source": raw.get("receptor_source", "deposited"), "mode": "supplied", "input_cif": str(cif),
                "sequences": {}, "command": None, "version": None}

    import mrcfile
    import numpy as np

    sequences = {kind: str((base / raw[f"{kind}_fasta"]).resolve()) for kind in ("protein", "rna", "dna") if raw.get(f"{kind}_fasta")}
    python = os.environ.get("LIGANDSEEK_CA2_PYTHON", sys.executable)
    version = subprocess.check_output([python, "-c", "from importlib.metadata import version; print(version('CryoAtom2'))"], text=True).strip()
    work = output / "cryoatom2_work"
    work.mkdir()
    native = experiment
    if experiment.suffix.lower() == ".gz":
        native = work / "experiment.mrc"
        with gzip.open(experiment, "rb") as compressed, native.open("wb") as unpacked:
            shutil.copyfileobj(compressed, unpacked)
    if raw.get("map_convention", "mrc") == "mrc":
        # CA2 2.1.1 会将 ORIGIN 与按轴重排的 NSTART 相加再乘以体素间距. 标准 MRC 的非零 ORIGIN 本身已是首采样点的 XYZ Å.
        with mrcfile.open(native, header_only=True) as handle:
            origin = np.array([handle.header.origin.x, handle.header.origin.y, handle.header.origin.z])
            spacing = np.array([handle.voxel_size.x, handle.voxel_size.y, handle.voxel_size.z])
        if np.any(origin):
            if native == experiment:
                native = work / "experiment.mrc"
                shutil.copyfile(experiment, native)
            # 只在交给 CA2 的副本中改写 header 单位. 原密度、轴序和所有采样值保持不变.
            with mrcfile.mmap(native, mode="r+") as handle:
                handle.header.origin = tuple(origin / spacing)
                handle.header.nxstart = handle.header.nystart = handle.header.nzstart = 0

    receptor = work / "cryoatom2"
    command = [python, "-u", "-m", "CryoAtom2", "build", "--map-path", str(native), "--output-dir", str(receptor), "--device", "cuda:0"]
    for kind, fasta in sequences.items():
        command.extend([f"--{kind}-sequence-path", fasta])
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join((str(Path(shutil.which(python) or python).parent), env.get("PATH", "")))
    started = time.monotonic()
    with (output / "cryoatom2.log").open("w", encoding="utf-8") as log:
        print(json.dumps({"command": command, "version": version}, ensure_ascii=False), file=log, flush=True)
        completed = subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, env=env)
        elapsed = time.monotonic() - started
        print(json.dumps({"elapsed_seconds": elapsed, "returncode": completed.returncode}), file=log, flush=True)
    if completed.returncode:
        print((output / "cryoatom2.log").read_text(encoding="utf-8", errors="replace")[-8192:], file=sys.stderr)
        completed.check_returncode()
    # 最终过滤后的 CIF 名称由官方入口取输出目录的末段得到.
    shutil.copyfile(receptor / "cryoatom2.cif", output / "receptor.cif")
    shutil.rmtree(work)
    return {"source": "ca2", "mode": "with_sequences" if sequences else "without_sequences", "input_cif": None,
            "sequences": sequences, "command": command, "version": version, "elapsed_seconds": elapsed}
