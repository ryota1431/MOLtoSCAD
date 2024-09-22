import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def select_mol_file():
    root = tk.Tk()
    root.withdraw()  # メインウィンドウを表示しない
    file_path = filedialog.askopenfilename(
        title="MOLファイルを選択してください",
        filetypes=[("MOL files", "*.mol"), ("All files", "*.*")]
    )
    return file_path

def add_hydrogens(mol):
    """
    分子に明示的な水素原子を追加し、3D座標を生成します。
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def mol_to_scad(mol_file, scad_file):
    # MOLファイルを読み込む
    with open(mol_file, 'r') as f:
        mol_block = f.read()
    mol = Chem.MolFromMolBlock(mol_block)
    
    if mol is None:
        print(f"エラー: MOLファイル '{mol_file}' を読み込めません。ファイルが存在し、正しい形式であることを確認してください。")
        return

    # 水素原子を追加し、3D座標を生成または最適化
    mol = add_hydrogens(mol)
    
    # SCADファイルを作成
    with open(scad_file, 'w') as f:
        f.write("// Generated SCAD file from MOL\n\n")
        
        # 変数の定義
        f.write("// Adjustable Parameters\n")
        f.write("atom_scale = 1;  // Adjusting the size of atoms\n")
        f.write("bond_radius = 0.2;  // Coupling Radius\n")
        f.write("\n")
        
        # 原子の色と半径の定義
        f.write("// Atomic color\n")
        for element, color in get_atom_color.colors.items():
            f.write(f"{element}_color = \"{color}\";\n")
        f.write("\n")
        
        f.write("// Atomic radius\n")
        for element, radius in get_atom_radius.radii.items():
            f.write(f"{element}_radius = {radius};\n")
        f.write("\n")
        
        # 原子を球として表現
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            element = atom.GetSymbol()
            
            f.write(f"translate([{pos.x}, {pos.y}, {pos.z}]) ")
            f.write(f"color({element}_color) ")
            f.write(f"sphere(r = {element}_radius * atom_scale);\n")
        
        # 結合を円柱として表現
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            pos1 = mol.GetConformer().GetAtomPosition(atom1.GetIdx())
            pos2 = mol.GetConformer().GetAtomPosition(atom2.GetIdx())
            
            f.write(f"color(\"grey\") ")
            f.write(f"hull() {{\n")
            f.write(f"  translate([{pos1.x}, {pos1.y}, {pos1.z}]) sphere(r = bond_radius);\n")
            f.write(f"  translate([{pos2.x}, {pos2.y}, {pos2.z}]) sphere(r = bond_radius);\n")
            f.write(f"}}\n")

def get_atom_color(element):
    return get_atom_color.colors.get(element, 'pink')  # デフォルトを桃色に変更

get_atom_color.colors = {
    'H': 'white',
    'C': 'black',
    'N': 'navy',
    'O': 'red',
    'F': 'green',
    'Cl': 'green',
    'Br': 'brown',
    'I': 'darkviolet',
    'He': 'cyan',
    'Ne': 'cyan',
    'Ar': 'cyan',
    'Xe': 'cyan',
    'Kr': 'cyan',
    'P': 'orange',
    'S': 'yellow',
    'B': 'pink',
    'Li': 'purple',
    'Na': 'purple',
    'K': 'purple',
    'Rb': 'purple',
    'Cs': 'purple',
    'Be': 'darkgreen',
    'Mg': 'darkgreen',
    'Ca': 'darkgreen',
    'Sr': 'darkgreen',
    'Ba': 'darkgreen',
    'Ra': 'darkgreen',
    'Ti': 'gray',
    'Fe': 'orange',
}

def get_atom_radius(element):
    return get_atom_radius.radii.get(element, 0.5)

get_atom_radius.radii = {
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

def open_scad_file(scad_file):
    # OpenSCADのパスを指定します。環境に合わせて変更してください。
    openscad_path = r"C:\Program Files (x86)\OpenSCAD\openscad.exe"
    
    # SCADファイルを開く
    try:
        subprocess.Popen([openscad_path, scad_file])
        print(f"OpenSCADで {scad_file} を開きました。")
    except FileNotFoundError:
        print("OpenSCADが見つかりません。パスを確認してください。")
    except Exception as e:
        print(f"エラーが発生しました: {e}")

def open_freecad_with_scad(scad_file):
    # FreeCADのパスを指定します。環境に合わせて変更してください。
    freecad_path = r"C:\Program Files\FreeCAD 0.21\bin\FreeCAD.exe"
    
    # FreeCADでSCADファイルを開く
    try:
        subprocess.Popen([freecad_path, scad_file])
        print(f"FreeCADで {scad_file} を開きました。")
    except FileNotFoundError:
        print("FreeCADが見つかりません。パスを確認してください。")
    except Exception as e:
        print(f"エラーが発生しました: {e}")

if __name__ == "__main__":
    # MOLファイルの選択
    mol_file = select_mol_file()
    
    if mol_file:
        # 出力ファイル名の生成（入力ファイルと同じ名前で拡張子を.scadに変更）
        base_name = os.path.splitext(os.path.basename(mol_file))[0]
        scad_file = f"{base_name}.scad"
        
        # MOLファイルからSCADファイルを生成
        mol_to_scad(mol_file, scad_file)
        print(f"SCADファイルが生成されました: {scad_file}")
        
        # ユーザーに開くアプリケーションを選択させる
        root = tk.Tk()
        root.withdraw()
        choice = messagebox.askyesnocancel("アプリケーション選択", "生成されたSCADファイルを開きますか？\n\nYes: OpenSCADで開く\nNo: FreeCADで開く\nCancel: 開かない")
        
        if choice is True:
            open_scad_file(scad_file)
        elif choice is False:
            open_freecad_with_scad(scad_file)
        else:
            print("ファイルは生成されましたが、自動的には開きません。")
    else:
        print("ファイルが選択されませんでした。")