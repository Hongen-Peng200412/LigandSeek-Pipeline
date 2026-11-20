"""从候选自身的 CCD 身份和分子图判断可识别的配体类别. 命令读入 SMILES 身份列表, 返回同序类别.

外部已提供的 type_tag 优先. 无法确定的类别为 null, 由自动对接范围按未知条件放行.
"""

import json
import sys


def identify(identities: list[dict]) -> list[str | None]:
    """返回与输入 identities 同序的 I 个类别, 如 [Mg+2] 为 ion, CCD NAG 为 sugar.

    - identities: 长度 I 的候选列表. smiles 是精确身份, ccd_id 是可选通用组分名, type_tag 是可选用户类别.
    - 返回值: 可确定的 sugar、ion、peptide_like、nucleotide_like 或用户声明类别. None 表示缺少分类证据.
    """
    import gemmi
    from rdkit import Chem
    from adaligand_preprocessing.stages.stage_c.constants import METAL_ELEMENTS

    tags = []
    # Gemmi 的通用残基表按 CCD 名称查询化学类别, 独立于本次样本的沉积结构.
    kind_to_tag = {"PYR": "sugar", "KET": "sugar", "RNA": "nucleotide_like", "DNA": "nucleotide_like",
                   "AA": "peptide_like", "AAD": "peptide_like", "MAA": "peptide_like", "PAA": "peptide_like"}
    for identity in identities:
        tag = identity.get("type_tag")
        if tag is None and identity.get("ccd_id"):
            tag = kind_to_tag.get(gemmi.find_tabulated_residue(identity["ccd_id"]).kind.name)
        if tag is None:
            mol = Chem.MolFromSmiles(identity["smiles"])
            if mol is not None and mol.GetNumAtoms() == 1 and mol.GetAtomWithIdx(0).GetSymbol().upper() in METAL_ELEMENTS:
                tag = "ion"
            elif mol is not None:
                # 仅识别由含氧糖环和羟基侧链组成的结构. 含磷酸基、芳环或非糖骨架的结合物保留未知.
                rings = mol.GetRingInfo().AtomRings()
                if rings and all(a.GetSymbol() in {"C", "O"} for a in mol.GetAtoms()):
                    sugar_rings = []
                    for ring in rings:
                        atoms = [mol.GetAtomWithIdx(i) for i in ring]
                        oxygens = sum(a.GetSymbol() == "O" for a in atoms)
                        hydroxyl_sites = sum(any(n.GetIdx() not in ring and n.GetSymbol() == "O" for n in a.GetNeighbors()) for a in atoms)
                        sugar_rings.append(len(ring) in (5, 6) and oxygens == 1 and hydroxyl_sites >= 3
                                          and all(not a.GetIsAromatic() for a in atoms))
                    ring_atoms = {i for ring in rings for i in ring}
                    carbons_outside = [a for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetIdx() not in ring_atoms]
                    if all(sugar_rings) and all(any(n.GetSymbol() == "O" for n in a.GetNeighbors()) for a in carbons_outside):
                        tag = "sugar"
        tags.append(tag)
    return tags


if __name__ == "__main__":
    print(json.dumps(identify(json.load(sys.stdin))))
