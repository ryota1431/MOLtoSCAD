import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def mol_to_scad(mol_file, scad_file):
    # molファイルを読み込む
    mol = Chem.MolFromMolFile(mol_file)
    
    # 既存の3D座標がない場合は生成
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
    
    # 3D座標の最適化
    AllChem.MMFFOptimizeMolecule(mol)
    
    # SCADファイルを作成
    with open(scad_file, 'w') as f:
        f.write("// Generated SCAD file from MOL\n\n")
        f.write("$fn = 32;\n\n")
        
        # 原子を球として表現
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            element = atom.GetSymbol()
            radius = get_atom_radius(element)
            
            f.write(f"translate([{pos.x}, {pos.y}, {pos.z}]) ")
            f.write(f"color(\"{get_atom_color(element)}\") ")
            f.write(f"sphere(r = {radius});\n")
        
        # 結合を円柱として表現
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            pos1 = mol.GetConformer().GetAtomPosition(atom1.GetIdx())
            pos2 = mol.GetConformer().GetAtomPosition(atom2.GetIdx())
            
            f.write(f"color(\"grey\") ")
            f.write(f"hull() {{\n")
            f.write(f"  translate([{pos1.x}, {pos1.y}, {pos1.z}]) sphere(r = 0.1);\n")
            f.write(f"  translate([{pos2.x}, {pos2.y}, {pos2.z}]) sphere(r = 0.1);\n")
            f.write(f"}}\n")

def get_atom_color(element):
    # 原子の色を定義（一部の元素のみ）
    colors = {
        'C': 'grey',
        'O': 'red',
        'N': 'blue',
        'H': 'white',
        'S': 'yellow'
    }
    return colors.get(element, 'purple')  # デフォルトは紫

def get_atom_radius(element):
    # 原子の半径を定義（一部の元素のみ）
    radii = {
        'H': 0.3,
        'C': 0.5,
        'N': 0.45,
        'O': 0.4,
        'F': 0.35,
        'Cl': 0.6,
        'Br': 0.7,
        'I': 0.8,
        'S': 0.55,
        'P': 0.55,
    }
    return radii.get(element, 0.5)  # デフォルトは0.5

# 使用例
mol_to_scad('flarene80.mol', 'output.scad')
