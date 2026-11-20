"""从本地固定版 wwPDB 组分词典补充 CCD 缓存. materialize_components 为受体和完整真值使用相同的原子与键模板."""

import os
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory


def materialize_components(components: set[str], cache: Path) -> None:
    """把缺少的本次组分转换为历史预处理器读取的 RDKit 分子, 使用同版解析器与 sanitize=False.

    - components: 本次受体或真实配体中出现的组分编号集合, 如 {"ATP", "ALA"}.
    - cache/components.cif.gz: 安装时固定的 wwPDB 词典, 每个 data 块对应一个组分.
    - cache/<组分>.pkl: RDKit Mol, 原子名、键及参考构象来自该 CCD 块, PDB_NAME 为大写组分号. 已有模板直接复用.
    """
    import gemmi
    from pdbeccdutils.core import ccd_reader
    from rdkit import Chem

    missing = sorted(code.upper() for code in components if not (cache / (code.upper() + ".pkl")).is_file())
    if not missing:
        return
    # 只有遇到缺失模板时读取压缩词典. 每次以独立目录暂存解析结果, 相同组分并发写入也只发布完整文件.
    dictionary = gemmi.cif.read(str(cache / "components.cif.gz"))
    with TemporaryDirectory(dir=cache, prefix=".ccd-") as scratch:
        for code in missing:
            block = dictionary.find_block(code)
            if block is None:
                continue
            source = Path(scratch) / "component.cif"
            block.write_file(str(source))
            molecule = ccd_reader.read_pdb_cif_file(str(source), sanitize=False).component.mol
            if molecule.GetNumAtoms():
                molecule.SetProp("PDB_NAME", code)
                packed = Path(scratch) / "component.pkl"
                with packed.open("wb") as stream:
                    # RDKit 默认 pickle 会丢失 CCD 原子名. 仅对这份输出启用 AllProps, 读取结果仍是普通 RDKit Mol.
                    writer = pickle.Pickler(stream)
                    writer.dispatch_table = {type(molecule): lambda mol: (Chem.Mol, (mol.ToBinary(Chem.PropertyPickleOptions.AllProps),))}
                    writer.dump(molecule)
                os.replace(packed, cache / (code + ".pkl"))
