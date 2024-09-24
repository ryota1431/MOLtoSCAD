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
    original_radius = get_atom_radius.original_radii.get(element, 0.77)  # デフォルト値を炭素の半径に設定
    return original_radius * get_atom_radius.normalization_factor

# オリジナルの半径データ 参考文献:https://qiita.com/chemweb000/items/83efb778cfdfe6d63ab4
# 参考文献の参考文献 CrystalMakerのHP https://www.hulinks.co.jp/support/c-maker/qa_005.html
get_atom_radius.original_radii = {
    'H': 0.37, 'He': 0.32,
    'Li': 1.34, 'Be': 0.90, 'B': 0.82, 'C': 0.77, 'N': 0.75, 'O': 0.73, 'F': 0.71, 'Ne': 0.69,
    'Na': 1.54, 'Mg': 1.30, 'Al': 1.18, 'Si': 1.11, 'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'Ar': 0.97,
    'K': 1.96, 'Ca': 1.74, 'Sc': 1.44, 'Ti': 1.36, 'V': 1.25, 'Cr': 1.27, 'Mn': 1.39, 'Fe': 1.25,
    'Co': 1.26, 'Ni': 1.21, 'Cu': 1.38, 'Zn': 1.31, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.19, 'Se': 1.16,
    'Br': 1.14, 'Kr': 1.10,
    'Rb': 2.11, 'Sr': 1.92, 'Y': 1.62, 'Zr': 1.48, 'Nb': 1.37, 'Mo': 1.45, 'Tc': 1.56, 'Ru': 1.26,
    'Rh': 1.35, 'Pd': 1.31, 'Ag': 1.53, 'Cd': 1.48, 'In': 1.44, 'Sn': 1.41, 'Sb': 1.38, 'Te': 1.35,
    'I': 1.33, 'Xe': 1.30,
    'Cs': 2.11, 'Ba': 1.92, 'Ce': 1.62, 'Pr': 1.48, 'Nd': 1.37, 'Pm': 1.45, 'Sm': 1.56, 'Eu': 1.26,
    'Gd': 1.35, 'Tb': 1.31, 'Dy': 1.53, 'Ho': 1.48, 'Er': 1.44, 'Tm': 1.41, 'Yb': 1.38, 'Lu': 1.35,
    'Hf': 1.33, 'Ta': 1.30, 'W': 1.30, 'Re': 1.30, 'Os': 1.30, 'Ir': 1.30, 'Pt': 1.30, 'Au': 1.30,
    'Hg': 1.30, 'Tl': 1.30, 'Pb': 1.30, 'Bi': 1.30, 'Po': 1.30, 'At': 1.30, 'Rn': 1.30,
    'Fr': 2.11, 'Ra': 1.92, 'Ac': 1.62, 'Pa': 1.48, 'U': 1.37, 'Np': 1.45, 'Pu': 1.56, 'Am': 1.26,
    'Cm': 1.35, 'Bk': 1.31, 'Cf': 1.53, 'Es': 1.48, 'Fm': 1.44, 'Md': 1.41, 'No': 1.38, 'Lr': 1.35
}

# 炭素の半径
carbon_radius = get_atom_radius.original_radii['C']

# 規格化係数を計算（炭素の半径を0.5にするための係数）
get_atom_radius.normalization_factor = 0.5 / carbon_radius

# 規格化された半径を計算して新しい辞書に格納
get_atom_radius.radii = {element: get_atom_radius(element) for element in get_atom_radius.original_radii}

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
        if mol_to_scad(mol_file, scad_file):
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
            print("SCADファイルの生成に失敗しました。")
    else:
        print("ファイルが選択されませんでした。")