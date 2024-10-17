# MOL to SCAD コンバーター

このプロジェクトは、MOLファイル形式の分子構造データをOpenSCAD用のSCADファイルに変換するPythonスクリプトです。

## 機能

- MOLファイルから分子構造を読み込む
- 3D座標を生成または最適化
- 原子を球として、結合を円柱としてSCADファイルに出力
- CPK配色に基づいて元素ごとに異なる色と大きさで表現
- 生成されたSCADファイルをOpenSCADまたはFreeCADで自動的に開くオプション

## 必要条件

- Python 3.6以上
- RDKit
- OpenSCAD（オプション）
- FreeCAD（オプション）

## Pythonのインストールと仮想環境の構築

### Pythonのインストール

1. [Python公式ウェブサイト](https://www.python.org/downloads/)から、最新版のPythonをダウンロードします。
2. ダウンロードしたインストーラーを実行し、指示に従ってインストールを完了します。
3. インストール時に「Add Python to PATH」オプションにチェックを入れてください。

### 仮想環境の構築

1. コマンドプロンプト（Windows）またはターミナル（macOS/Linux）を開きます。

2. プロジェクトのディレクトリに移動します：
   ```
   cd path/to/your/project
   ```

3. 仮想環境を作成します：
   ```
   python -m venv venv
   ```

4. 仮想環境を有効化します：
   - Windows:
     ```
     venv\Scripts\activate
     ```
   - macOS/Linux:
     ```
     source venv/bin/activate
     ```

## インストール

1. このリポジトリをクローンまたはダウンロードします。

2. 仮想環境を有効化した状態で、必要なパッケージをインストールします：
   ```
   pip install rdkit
   ```

## OpenSCADとFreeCADのパス設定

スクリプトがOpenSCADとFreeCADを正しく起動できるように、それぞれのアプリケーションのパスを設定する必要があります。

1. `main.py` ファイルをテキストエディタで開きます。

2. `open_scad_file` 関数内の `openscad_path` 変数を、お使いのシステムのOpenSCADの実行ファイルのパスに変更します：

   ```python
   openscad_path = r"C:\Program Files\OpenSCAD\openscad.exe"  # Windowsの場合
   # openscad_path = "/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD"  # macOSの場合
   # openscad_path = "/usr/bin/openscad"  # Linuxの場合
   ```

3. 同様に、`open_freecad_with_scad` 関数内の `freecad_path` 変数を、FreeCADの実行ファイルのパスに変更します：

   ```python
   freecad_path = r"C:\Program Files\FreeCAD 0.21\bin\FreeCAD.exe"  # Windowsの場合
   # freecad_path = "/Applications/FreeCAD.app/Contents/MacOS/FreeCAD"  # macOSの場合
   # freecad_path = "/usr/bin/freecad"  # Linuxの場合
   ```

4. 変更を保存します。

## 使用方法

1. コマンドプロンプトまたはターミナルを開きます。

2. プロジェクトのディレクトリに移動します：
   ```
   cd path/to/your/project
   ```

3. 仮想環境を有効化します（まだ有効化していない場合）：
   - Windows:
     ```
     venv\Scripts\activate
     ```
   - macOS/Linux:
     ```
     source venv/bin/activate
     ```

4. 以下のコマンドでスクリプトを実行します：
   ```
   python main.py
   ```

5. ファイル選択ダイアログが表示されるので、変換したいMOLファイルを選択します。

6. スクリプトが自動的にMOLファイルを処理し、SCADファイルを生成します。

7. 変換が完了すると、以下のオプションが表示されます：
   - OpenSCADでSCADファイルを開く
   - FreeCADでSCADファイルを開く
   - ファイルを開かずに終了する

8. 選択したオプションに応じて、対応するアプリケーションが起動します（OpenSCADまたはFreeCAD）。

注意：
- OpenSCADとFreeCADのパスが正しく設定されていることを確認してください。
- 大きな分子の場合、処理に時間がかかる場合があります。
- 生成されたSCADファイルは、プロジェクトのディレクトリに保存されます。

## カスタマイズ

- `get_atom_color` 関数を編集して、元素の色を変更できます。
- `get_atom_radius` 関数を編集して、元素の大きさを調整できます。
- OpenSCADとFreeCADのパスを環境に合わせて変更してください。

## 注意事項

- 大きな分子の場合、生成されるSCADファイルが非常に大きくなる可能性があります。
- 複雑な分子構造の場合、3D座標の生成に時間がかかることがあります。
- OpenSCADまたはFreeCADが正しくインストールされていることを確認してください。

## ライセンス

このプロジェクトは [MITライセンス](LICENCE) の下で公開されています。

## 貢献

バグ報告や機能リクエストは、GitHubのIssueを通じてお願いします。プルリクエストも歓迎します。

## プログラムの詳細な仕組み

### 1. プログラムの概要

このプログラムは、分子構造を表すMOLファイルを3Dプリント用のSCADファイルに変換するツールです。主な機能は以下の通りです：

1. MOLファイルの選択
2. 3D構造の生成と最適化
3. SCADファイルの生成
4. OpenSCADまたはFreeCADでの自動表示

### 2. 主要な関数と処理の詳細

#### 2.1 MOLファイルの選択

```python
def select_mol_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="MOLファイルを選択してください",
        filetypes=[("MOL files", "*.mol"), ("All files", "*.*")]
    )
    return file_path
```

この関数は、Tkinterを使用してGUIでMOLファイルを選択するダイアログを表示します。

#### 2.2 3D構造の生成と最適化

```python
def add_hydrogens(mol):
    mol = Chem.AddHs(mol)
    
    # 4つの異なる方法で3D構造の生成を試みます
    methods = [
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True, "useBasicKnowledge": False}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useRandomCoords": True, "boxSizeMult": 2.0}),
        (AllChem.EmbedMolecule, {"randomSeed": 42, "useExpTorsionAnglePrefs": True, "useBasicKnowledge": True})
    ]

    for method, params in methods:
        result = method(mol, **params)
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol)
            return mol

    raise ValueError("すべての3D構造生成方法が失敗しました。")
```

この関数は、RDKitを使用して分子の3D構造を生成し最適化します。4つの異なる方法を順番に試み、成功した場合はMMFF力場で構造を最適化します。

#### 2.3 SCADファイルの生成

```python
def mol_to_scad(mol_file, scad_file, max_atoms_per_file=1000):
    # MOLファイルの読み込みと前処理
    mol = Chem.MolFromMolFile(mol_file, removeHs=False)
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)
    mol = add_hydrogens(mol)

    # SCADファイルの生成
    atom_count = mol.GetNumAtoms()
    file_count = (atom_count - 1) // max_atoms_per_file + 1

    for file_index in range(file_count):
        # ファイル名の設定
        current_scad_file = f"{os.path.splitext(scad_file)[0]}_{file_index + 1}.scad" if file_count > 1 else scad_file

        with open(current_scad_file, 'w') as f:
            # SCADファイルのヘッダー部分を書き込み
            write_scad_header(f)

            # 原子と結合の情報を書き込み
            write_atoms_and_bonds(f, mol, file_index * max_atoms_per_file, (file_index + 1) * max_atoms_per_file)

    return True
```

この関数は、MOLファイルを読み込み、3D構造を生成し、SCADファイルを作成します。大きな分子の場合、複数のファイルに分割します。

#### 2.4 OpenSCADまたはFreeCADでの自動表示

```python
def open_scad_file(scad_file):
    openscad_path = r"C:\Program Files (x86)\OpenSCAD\openscad.exe"
    try:
        subprocess.Popen([openscad_path, scad_file])
        print(f"OpenSCADで {scad_file} を開きました。")
    except FileNotFoundError:
        print("OpenSCADが見つかりません。パスを確認してください。")
    except Exception as e:
        print(f"エラーが発生しました: {e}")

def open_freecad_with_scad(scad_file):
    freecad_path = r"C:\Program Files\FreeCAD 0.21\bin\FreeCAD.exe"
    try:
        subprocess.Popen([freecad_path, scad_file])
        print(f"FreeCADで {scad_file} を開きました。")
    except FileNotFoundError:
        print("FreeCADが見つかりません。パスを確認してください。")
    except Exception as e:
        print(f"エラーが発生しました: {e}")
```

これらの関数は、生成されたSCADファイルをOpenSCADまたはFreeCADで自動的に開きます。

### 3. メイン処理の流れ

```python
if __name__ == "__main__":
    mol_file = select_mol_file()
    if mol_file:
        base_name = os.path.splitext(os.path.basename(mol_file))[0]
        scad_file = f"{base_name}.scad"
        
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
```

メイン処理では、MOLファイルの選択、SCADファイルの生成、そしてユーザーの選択に応じたアプリケーションでの表示を行います。

### 4. 追加の注意点

1. OpenSCADとFreeCADのパスは、ユーザーの環境に合わせて変更する必要があります。
2. 大きな分子の場合、複数のSCADファイルが生成されることがあります。
3. 3D構造の生成に失敗した場合、エラーメッセージが表示されます。
4. このプログラムは、RDKit、Tkinter、およびsubprocessモジュールに依存しています。

このプログラムを使用することで、分子構造を簡単に3Dモデル化し、3Dプリンターで出力することができます。
