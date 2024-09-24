import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

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
    複数の方法を試みて失敗した場合は例外を発生させます。
    """
    mol = Chem.AddHs(mol)
    
    # 方法1: ETKDG
    result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    if result == 0:
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    # 方法2: ETDG
    result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True, useBasicKnowledge=False)
    if result == 0:
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    # 方法3: スケーリングされたETKDG
    result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True, boxSizeMult=2.0)
    if result == 0:
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    # 方法4: 距離ジオメトリ
    result = AllChem.EmbedMolecule(mol, randomSeed=42, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    if result == 0:
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    raise ValueError("すべての3D構造生成方法が失敗しました。")

def mol_to_scad(mol_file, scad_file, max_atoms_per_file=1000):
    # MOLファイルを読み込む
    with open(mol_file, 'r') as f:
        mol_block = f.read()
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
    
    if mol is None:
        print(f"エラー: MOLファイル '{mol_file}' を読み込めません。ファイルが存在し、正しい形式であることを確認してください。")
        return False
    
    try:
        # 分子の前処理
        mol = Chem.RemoveHs(mol)  # 既存の水素を除去
        mol = Chem.AddHs(mol)  # 水素を再追加
        
        # 3D構造の生成
        mol = add_hydrogens(mol)
    except Exception as e:
        print(f"エラー: 3D構造の生成に失敗しました: {e}")
        return False
    
    # 3D座標が生成されたことを確認
    conf = mol.GetConformer()
    if not conf.Is3D():
        print("エラー: 3D座標が生成されませんでした。")
        print("SCADファイルの生成をスキップします。")
        return False
    
    atom_count = mol.GetNumAtoms()
    file_count = (atom_count - 1) // max_atoms_per_file + 1
    
    for file_index in range(file_count):
        start_atom = file_index * max_atoms_per_file
        end_atom = min((file_index + 1) * max_atoms_per_file, atom_count)
        
        if file_count > 1:
            current_scad_file = f"{os.path.splitext(scad_file)[0]}_{file_index + 1}.scad"
        else:
            current_scad_file = scad_file
        
        try:
            f = open(current_scad_file, 'w')
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
            for atom_idx in range(start_atom, end_atom):
                atom = mol.GetAtomWithIdx(atom_idx)
                pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
                element = atom.GetSymbol()
                
                f.write(f"translate([{pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f}]) ")
                f.write(f"color({element}_color) ")
                f.write(f"sphere(r = {element}_radius * atom_scale);\n")
            
            # 結合を円柱として表現
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                if start_atom <= atom1_idx < end_atom or start_atom <= atom2_idx < end_atom:
                    atom1 = mol.GetAtomWithIdx(atom1_idx)
                    atom2 = mol.GetAtomWithIdx(atom2_idx)
                    pos1 = mol.GetConformer().GetAtomPosition(atom1.GetIdx())
                    pos2 = mol.GetConformer().GetAtomPosition(atom2.GetIdx())
                    
                    f.write(f"color(\"grey\") ")
                    f.write(f"hull() {{\n")
                    f.write(f"  translate([{pos1.x:.2f}, {pos1.y:.2f}, {pos1.z:.2f}]) sphere(r = bond_radius);\n")
                    f.write(f"  translate([{pos2.x:.2f}, {pos2.y:.2f}, {pos2.z:.2f}]) sphere(r = bond_radius);\n")
                    f.write(f"}}\n")
            
            f.flush()
        finally:
            f.close()
    
    sys.stdout.flush()

    if file_count > 1:
        print(f"{file_count}個のSCADファイルが生成されました。")
    else:
        print(f"SCADファイルが生成されました: {scad_file}")

    return True

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