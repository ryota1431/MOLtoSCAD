import os
import tkinter as tk
from tkinter import filedialog, messagebox
from moltoscad import mol_to_scad, open_scad_file, open_freecad_with_scad


def select_mol_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="MOLファイルを選択してください",
        filetypes=[("MOL files", "*.mol"), ("All files", "*.*")],
    )
    return file_path


if __name__ == "__main__":
    mol_file = select_mol_file()

    if mol_file:
        base_name = os.path.splitext(os.path.basename(mol_file))[0]
        scad_file = f"{base_name}.scad"

        try:
            mol_to_scad(mol_file, scad_file)
        except Exception as e:
            print(f"SCADファイルの生成に失敗しました: {e}")
        else:
            print(f"SCADファイルが生成されました: {scad_file}")
            root = tk.Tk()
            root.withdraw()
            choice = messagebox.askyesnocancel(
                "アプリケーション選択",
                "生成されたSCADファイルを開きますか？\n\nYes: OpenSCADで開く\nNo: FreeCADで開く\nCancel: 開かない",
            )
            if choice is True:
                try:
                    open_scad_file(scad_file)
                except FileNotFoundError as e:
                    print(e)
            elif choice is False:
                try:
                    open_freecad_with_scad(scad_file)
                except FileNotFoundError as e:
                    print(e)
            else:
                print("ファイルは生成されましたが、自動的には開きません。")
    else:
        print("ファイルが選択されませんでした。")
