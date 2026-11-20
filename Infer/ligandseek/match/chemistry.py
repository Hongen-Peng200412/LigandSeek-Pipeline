"""将本样本指定的精确 SMILES 转为 Matcher 与 DG-PocketXMol 共同读取的图和对称排列.

prepare_graphs 在 directory 内生成 smiles_graphs_v1.npz 及 smiles_symmetries_v1.npz. 两者的 smiles 轴均为成功解析的唯一身份, 排序与查找沿用原图读取器.
"""

import json
import os
from pathlib import Path
import subprocess
import sys


def prepare_graphs(smiles: list[str], directory: Path) -> list[dict]:
    """从当前给出的精确字符串构建重原子图, 使用 RDKit 的原子顺序.

    - smiles: 长度 I 的身份列表, 如 ["CCO", "[Mg+2]"]. 同一身份的多个请求共用该图.
    - directory/smiles_graphs_v1.npz: smiles 为排序后的身份字符串, atom_offsets 切分逐原子的 element、charge、atom_in_ring, bond_offsets 切分 bond_index 第二维及逐键的 bond_type、bond_in_ring.
    - directory/smiles_symmetries_v1.npz: atom_count 是逐身份重原子数, matches_offsets 切分展平的 matches_iso, 每图按 matches_shape 还原为 (排列数, 可变原子数).
    - 返回被化学生成器排除的身份列表, 每项的 smiles 是输入字符串, error 是解析或规范化原因. 例如非法字符串会留下排除原因.
    """
    # Find/Match 的 CUDA 环境可独立于 RDKit. 将精确身份作为 JSON 数组交给安装 RDKit 的解释器, 仅返回排除原因.
    interpreter = os.environ.get("LIGANDSEEK_PREPARE_PYTHON", sys.executable)
    process = subprocess.run(
        [interpreter, "-B", __file__, str(directory)],
        input=json.dumps(smiles), stdout=subprocess.PIPE, encoding="utf-8", check=True,
    )
    return json.loads(process.stdout)


if __name__ == "__main__":
    import numpy as np
    from rdkit import Chem
    from ops.matcher.build_smiles_assets import BOND_TYPES, build_arrays

    # 原构图器可记录无效字符串, 对词表外的键则会终止整批. 将这类身份从图输入移出, 如四重键 [Mo]$[Mo].
    supported = []
    excluded_bonds = []
    for smiles in sorted(set(json.load(sys.stdin))):
        mol = Chem.MolFromSmiles(smiles)
        missing = {str(b.GetBondType()) for b in mol.GetBonds() if b.GetBondType() not in BOND_TYPES} if mol else set()
        if missing:
            excluded_bonds.append({"smiles": smiles, "error": "unsupported bond: " + ", ".join(sorted(missing))})
        else:
            supported.append(smiles)
    # 每张图的原子端点和对称排列沿用 Matcher 的原始顺序. 排除的身份仍在外部 Match 结果保留原槽位.
    graph, symmetry, excluded = build_arrays(supported)
    directory = Path(sys.argv[1])
    directory.mkdir(parents=True, exist_ok=True)
    np.savez(directory / "smiles_graphs_v1.npz", **graph)
    np.savez(directory / "smiles_symmetries_v1.npz", **symmetry)
    print(json.dumps(excluded + excluded_bonds))
