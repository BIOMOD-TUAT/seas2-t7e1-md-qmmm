# wiki Fig. 3・7・8 用の MD 解析（WT T7E1 ホモ二量体、0A〜5A）

作成：2026-11-13T03:54:53。軌道を読むだけの解析（元のフォルダには何も書いていない）。
各条件 0〜50 ns、100 ps 間隔の 501 フレーム（3A は 500 ns の計算の最初の 50 ns）。3A の 500 ns は別の表（Fig08_3A_500ns_blocks_*.csv）。

## フォルダ
- 01_Analysis_CSV：Fig08_*、Fig07_*、Fig03_*、Threshold_Sensitivity.csv（集計、フレームごと、ヒストグラム）
- 02_Excel、03_Figures：この環境では作っていない（CSV から別に作る）
- 04_3D_Structures/Representative_PDB：代表の構造と representative_frames.csv
- 05_Analysis_Scripts：使ったスクリプト
- 06_Methods_and_QC：input_manifest.csv、atom_mapping、Mn_site_assignment.csv、analysis_definition.json、validation_report.md、この README

## 主張しないこと
- 3A が最も切断されやすい、とは言わない。
- NAC-like の形成率が、切断の速度や活性の順位を表す、とは言わない。
- 条件の間の差は、初期構造や平衡化の履歴の違いを含む、探索的な比較である。各条件は 1 本の軌道。
- Level 3 は、プロトン移動が起きたことを意味しない。3A の 500 ns の時間のブロックは独立な標本ではない。
- Mn は非結合のモデルで、配位の幾何は結晶構造と違うことがある。
