"""在实验采样网格上生成受体密度. prepare_density 写 exp.npy、sim.npy 和各自的世界 XYZ 几何 NPZ, 返回分辨率和外部程序记录."""

import json
import os
from pathlib import Path
import re
import shutil


def prepare_density(raw: dict, experiment: Path, receptor: Path, base: Path, output: Path) -> dict:
    """使用 AdaLigand 的 Fourier 重采样与 Chimera Classic 的 molmap onGrid.

    - raw: 原始输入参数. resolution 单位 Å, emdb_metadata、rcsb_metadata 指向同一实验的元数据 JSON. map_convention 可取 mrc 或 adaligand.
    - experiment: 所给 MRC/map(.gz) 的路径. receptor 是世界 XYZ 坐标相同的 mmCIF.
    - base: raw.json 所在目录, 用于解析元数据相对路径.
    - output: exp.npy 和 sim.npy 均为 float32 (1,Z,Y,X), 第零维是单通道. 两个同名 NPZ 的 origin、voxel_size 为 float32 (3,) 世界 XYZ 下角点和实际间距, 单位 Å. canonical_shape_zyx 是三维网格形状.
    - 返回值: resolution 是生成模拟图采用的 Å 值, resolution_source 指明实际依据, map_convention 是原点解释规则, chimera 是实测版本, atom_count 是模拟使用的 ATOM 重原子数.
    """
    import mrcfile
    import numpy as np
    from adaligand_preprocessing.external_tools.chimera import ChimeraRunner
    from adaligand_preprocessing.external_tools.model_cif import write_normalized_model_cif
    from adaligand_preprocessing.geometry.mrc import MapGrid, load_map, make_canonical_grid, write_canonical_mrc
    from ligandseek.find.inputs import mrc_first_sample

    resolution, source = raw.get("resolution"), "user"
    if resolution is None and (raw.get("emdb_metadata") or raw.get("rcsb_metadata")):
        from adaligand_preprocessing.stages.stage_a import extract_resolution_info

        metadata = {kind: json.loads((base / raw[kind]).read_text(encoding="utf-8")) if raw.get(kind) else {}
                    for kind in ("rcsb_metadata", "emdb_metadata")}
        info = extract_resolution_info(metadata["rcsb_metadata"], metadata["emdb_metadata"])
        resolution, source = info["selected"], info["selected_source"]
    if resolution is None:
        # 只读取明确标为 resolution 且带 Å 单位的 header 标签. 体素间距属于另一物理量.
        with mrcfile.open(experiment, header_only=True) as handle:
            labels = [bytes(s).decode("utf-8", errors="replace").strip() for s in handle.header.label[:int(handle.header.nlabl)]]
        for label in labels:
            match = re.search(r"\bresolution\s*[:=]?\s*(\d+(?:\.\d+)?)\s*(?:angstroms?|Å|A)\b", label, flags=re.I)
            if match and float(match[1]) > 0:
                resolution, source = float(match[1]), "mrc_label"
                break
    if resolution is None:
        resolution, source = 3.0, "default"

    convention = raw.get("map_convention", "mrc")
    native = load_map(experiment, multiply_global_origin=convention == "adaligand")
    grid = make_canonical_grid(native, target_voxel_size=1.0)
    spacing = grid.voxel_size
    origin = grid.origin
    chimera_origin = origin
    if convention == "mrc":
        # (3,), 奇数轴在较小索引处补的原始体素数, 按 XYZ 排列. Fourier 重采样保留补齐后第一个采样点的相位.
        pad_xyz = np.asarray(native.grid.shape[::-1]) % 2
        first = mrc_first_sample(experiment, native.voxel_size) - pad_xyz * native.voxel_size
        origin = np.asarray(first - 0.5 * spacing, dtype=np.float32)
        # 标准 MRC 向 Chimera 传首采样点; 历史 adaligand 条件沿用生产适配器写出的原点.
        chimera_origin = np.asarray(first, dtype=np.float32)
    del native
    geometry = {"origin": origin, "voxel_size": spacing, "canonical_shape_zyx": np.asarray(grid.grid.shape, dtype=np.int32),
                "schema_version": np.uint16(2)}
    np.save(output / "exp.npy", grid.grid[None], allow_pickle=False)
    np.savez(output / "exp.npz", **geometry)

    runner = ChimeraRunner([os.environ.get("LIGANDSEEK_CHIMERA", "chimera")])
    work = output / "chimera_work"
    work.mkdir()
    stats = write_normalized_model_cif(receptor, work / "receptor.cif", atom_only=True)
    write_canonical_mrc(work / "grid.mrc", MapGrid(grid.grid, spacing, chimera_origin))
    version = runner.probe(work)
    runner.molmap_on_grid(work / "receptor.cif", work / "grid.mrc", work / "density.mrc", resolution=float(resolution), scratch_dir=work)
    simulated = load_map(work / "density.mrc", multiply_global_origin=False)
    # 两个场必须在同一组物理采样点上. 比较的是 Chimera 实际生成的网格, 然后给两份数组保存相同的角点和间距.
    np.testing.assert_array_equal(simulated.grid.shape, grid.grid.shape)
    np.testing.assert_allclose(simulated.voxel_size, spacing, rtol=0, atol=1e-5)
    np.testing.assert_allclose(simulated.origin, chimera_origin, rtol=0, atol=1e-5)
    np.save(output / "sim.npy", simulated.grid[None], allow_pickle=False)
    np.savez(output / "sim.npz", **geometry, resolution=np.float32(resolution), chimera_version=np.asarray(version))
    shutil.rmtree(work)
    return {"resolution": float(resolution), "resolution_source": source, "map_convention": convention, "chimera": version, "atom_count": stats["n_atoms"]}
