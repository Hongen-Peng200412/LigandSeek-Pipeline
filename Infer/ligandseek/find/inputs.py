"""从本次样本建立三个模型共用的受体表和密度图.

prepare_inputs 将实验图和同受体模拟图写入 inputs_root/density/<样本>/, 受体重原子表写入 inputs_root/parse/<样本>/receptor_tokens.npz.
"""

import os
import subprocess
import sys
from pathlib import Path
from shutil import copyfile


def mrc_first_sample(path: Path, voxel_xyz):
    """读正交实空间 MRC 首个密度采样点的世界 XYZ 坐标, 返回长度 3 的 Å 向量.

    path 可指向 gzip 文件. voxel_xyz 是 XYZ 三轴间距 (Å/体素). 本包采用 ChimeraX 的采样点约定.
    非零 ORIGIN 给出 XYZ Å. 当三个分量均为零时, 用列、行、层的整数 NSTART 与 MAPC/MAPR/MAPS 计算起点.
    """
    import mrcfile
    import numpy as np

    with mrcfile.open(path, header_only=True) as handle:
        h = handle.header
        origin = np.array([h.origin.x, h.origin.y, h.origin.z], dtype=np.float64)
        if np.any(origin):
            return origin
        axis = np.array([h.mapc, h.mapr, h.maps]) - 1
        start_xyz = np.empty(3, dtype=np.int64)
        start_xyz[axis] = [h.nxstart, h.nystart, h.nzstart]
        return start_xyz * voxel_xyz


def prepare_density(experimental: Path, simulated: Path, destination: Path, convention: str | None = None):
    """把一对用户密度图写为同一受体条件的模型资产. destination 是已创建的单样本密度目录.

    experimental 与 simulated 可为 MRC/map(.gz) 或带同名几何 NPZ 的 NPY. convention 可选 mrc 或 adaligand.
    缺省时模拟 NPY 表示沿用训练导入规则, 模拟 MRC 表示两张原图都采用 MRC 采样点约定.
    输出 exp/sim.npy 为 float32 (1,Z,Y,X). 同名 NPZ 中 origin 是网格下角点, voxel_size 是实际间距, 均为 XYZ (3,) Å.
    """
    import numpy as np
    from adaligand_preprocessing.geometry.mrc import load_map, make_canonical_grid

    convention = convention or ("adaligand" if simulated.suffix.lower() == ".npy" else "mrc")
    for name, source in (("exp", experimental), ("sim", simulated)):
        if source.suffix.lower() == ".npy":
            # 已在模型网格上的 NPY/NPZ 直接复制, 包括浮点几何和原有元数据.
            copyfile(source, destination / f"{name}.npy")
            copyfile(source.with_suffix(".npz"), destination / f"{name}.npz")
            continue
        legacy_origin = convention == "adaligand" and name == "exp"
        native = load_map(source, multiply_global_origin=legacy_origin)
        resample = True
        if name == "sim":
            with np.load(destination / "exp.npz", allow_pickle=False) as geometry:
                exp_shape = tuple(geometry["canonical_shape_zyx"])
                exp_spacing, exp_origin = geometry["voxel_size"], geometry["origin"]
            resample = not (native.grid.shape == exp_shape and np.allclose(native.voxel_size, exp_spacing, rtol=0, atol=1e-5))
        grid = make_canonical_grid(native, target_voxel_size=1.0) if resample else native
        origin, spacing = grid.origin, grid.voxel_size
        if convention == "mrc":
            # make_cubic 在奇数长度的轴前方补一个原体素. Fourier 变换保留补齐后第零采样点的相位.
            pad_xyz = np.asarray(native.grid.shape[::-1]) % 2 if resample else np.zeros(3, dtype=int)
            first_xyz = mrc_first_sample(source, native.voxel_size) - pad_xyz * native.voxel_size
            # Find 从下角点出发, 用 origin+(i+0.5)*spacing 定位采样点, 此处使用重采样后的实际间距.
            origin = (first_xyz - 0.5 * spacing).astype(np.float32)
        if name == "sim":
            # 同网格图的单精度舍入按原受体适配器的 1e-5 Å 容差归一. 实际位置偏移保留在输出几何中.
            spacing = exp_spacing if np.allclose(spacing, exp_spacing, rtol=0, atol=1e-5) else spacing
            origin = exp_origin if np.allclose(origin, exp_origin, rtol=0, atol=1e-5) else origin
        np.save(destination / f"{name}.npy", grid.grid[None], allow_pickle=False)
        np.savez(destination / f"{name}.npz", origin=origin, voxel_size=spacing,
                 schema_version=np.uint16(2), canonical_shape_zyx=np.asarray(grid.grid.shape, dtype=np.int32))


# ================================================================================================
def prepare_inputs(sample: dict, inputs_root: Path, ccd_directory: Path) -> None:
    """在训练使用的实验网格上保存实验密度、对应模拟密度和受体原子表.

    输入:
        - sample: 单样本字典. sample_id 标识其受体条件, 如 "9ter". experimental_map 与 simulated_map 是实验图与对应模拟图, receptor_cif 提供指定受体. 可选 map_convention 为 mrc 或 adaligand, 见 prepare_density.
        - inputs_root: 本次运行的资产目录, 如 find/inputs. 使用小写 sample_id 作为子目录名.
        - ccd_directory: 固定化学组分词典的 RDKit 模板缓存, 如 models/ccd, 用于受体残基的键和原子名.

    环境:
        - LIGANDSEEK_PREPARE_PYTHON: 可选的化学预处理 Python 解释器. 子进程接收两张图、原点约定、受体、样本标识、模板目录与输出目录.

    输出:
        - density/<样本>/{exp,sim}.npy: float32 (1,Z,Y,X), 单通道实验/模拟密度. 两图共用实验世界坐标.
        - density/<样本>/{exp,sim}.npz: schema_version=2, origin 和 voxel_size 均为 float32 (3,), 对应世界 XYZ 网格下角点和实际间距, 单位 Å. exp 的 canonical_shape_zyx 为整图 ZYX 形状.
        - parse/<样本>/receptor_tokens.npz: N 个重原子的 coords (N,3) 为世界 XYZ Å, feat (N,49) 为模型输入, element、res_type、is_backbone、atom_name、res_index、chain_index (N,) 均按 coords 顺序. bond_index (2,E) 每列给出键两端的原子下标, bond_type (E,) 同序.
    """
    interpreter = os.environ.get("LIGANDSEEK_PREPARE_PYTHON")
    if interpreter and os.path.abspath(interpreter) != os.path.abspath(sys.executable):
        subprocess.run([interpreter, "-B", __file__, sample["sample_id"], sample["experimental_map"], sample["receptor_cif"],
                        sample["simulated_map"], inputs_root, ccd_directory, sample.get("map_convention") or ""], check=True)
        return

    import gemmi
    import numpy as np
    from adaligand_preprocessing.stages.stage_c.pipeline import category_rows, optional_category_rows, selected_atom_rows, split_candidate_atoms
    from adaligand_preprocessing.stages.stage_c.receptor import build_receptor_arrays

    sample_id = sample["sample_id"].lower()
    density = inputs_root / "density" / sample_id
    receptor = inputs_root / "parse" / sample_id
    density.mkdir(parents=True, exist_ok=True)
    receptor.mkdir(parents=True, exist_ok=True)

    # 原解析器选择首模型和规范异构位置, 受体表只保留重原子.
    block = gemmi.cif.read(str(sample["receptor_cif"])).sole_block()
    atoms = selected_atom_rows(category_rows(block, "_atom_site."))
    entity_types = {str(row["id"]): str(row["type"]).lower() for row in optional_category_rows(block, "_entity.")}
    if entity_types:
        _, receptor_atoms = split_candidate_atoms(atoms, entity_types)
    else:
        # CryoAtom2 最终 CIF 通常省略 _entity. 原受体适配器按 ATOM 记录识别聚合物.
        receptor_atoms = [atom for atom in atoms if str(atom["group_PDB"]).upper() == "ATOM"]
    from ligandseek.ccd import materialize_components

    materialize_components({atom["label_comp_id"] for atom in receptor_atoms}, ccd_directory)
    tokens = build_receptor_arrays(receptor_atoms, optional_category_rows(block, "_struct_conn."), ccd_directory, allow_ccd_fetch=False)
    np.savez_compressed(receptor / "receptor_tokens.npz", **tokens)

    prepare_density(Path(sample["experimental_map"]), Path(sample["simulated_map"]), density, sample.get("map_convention"))


if __name__ == "__main__":
    sample_id, experiment, receptor, simulation, inputs_root, ccd_directory, *conventions = sys.argv[1:]
    sample = {"sample_id": sample_id, "experimental_map": experiment, "receptor_cif": receptor, "simulated_map": simulation}
    if conventions:
        sample["map_convention"] = conventions[0]
    prepare_inputs(sample, Path(inputs_root), Path(ccd_directory))
