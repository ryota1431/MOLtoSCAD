# MOL to SCAD コンバーター

このプロジェクトは、MOLファイル形式の分子構造データをOpenSCAD用のSCADファイルに変換するPythonスクリプトです。

## 機能

- MOLファイルから分子構造を読み込む
- 3D座標を生成または最適化
- 原子を球として、結合を円柱としてSCADファイルに出力
- 元素ごとに異なる色と大きさで表現

## 必要条件

- Python 3.6以上
- RDKit

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

2. 仮想環境を有効化した状態で、RDKitをインストールします：
   ```
   pip install rdkit
   ```

## 使用方法

1. `main.py` ファイルを開き、入力MOLファイルと出力SCADファイルのパスを設定します：

   ```python
   mol_to_scad('input.mol', 'output.scad')
   ```

2. スクリプトを実行します：

   ```
   python main.py
   ```

3. 生成されたSCADファイルをOpenSCADで開いて3Dモデルを表示します。

## カスタマイズ

- `get_atom_color` 関数を編集して、元素の色を変更できます。
- `get_atom_radius` 関数を編集して、元素の大きさを調整できます。

## 注意事項

- 大きな分子の場合、生成されるSCADファイルが非常に大きくなる可能性があります。
- 複雑な分子構造の場合、3D座標の生成に時間がかかることがあります。

## ライセンス

このプロジェクトは [MITライセンス](LICENSE) の下で公開されています。

## 貢献

バグ報告や機能リクエストは、GitHubのIssueを通じてお願いします。プルリクエストも歓迎します。
