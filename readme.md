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

1. `main.py` ファイルを開き、入力MOLファイルと出力SCADファイルのパスを設定します：

   ```python
   mol_to_scad('input.mol', 'output.scad')
   ```

2. スクリプトを実行します：

   ```
   python main.py
   ```

3. 生成されたSCADファイルをOpenSCADまたはFreeCADで開いて3Dモデルを表示します。

## カスタマイズ

- `get_atom_color` 関数を編集して、元素の色を変更できます。
- `get_atom_radius` 関数を編集して、元素の大きさを調整できます。
- OpenSCADとFreeCADのパスを環境に合わせて変更してください。

## 注意事項

- 大きな分子の場合、生成されるSCADファイルが非常に大きくなる可能性があります。
- 複雑な分子構造の場合、3D座標の生成に時間がかかることがあります。
- OpenSCADまたはFreeCADが正しくインストールされていることを確認してください。

## ライセンス

このプロジェクトは [MITライセンス](LICENSE) の下で公開されています。

## 貢献

バグ報告や機能リクエストは、GitHubのIssueを通じてお願いします。プルリクエストも歓迎します。
