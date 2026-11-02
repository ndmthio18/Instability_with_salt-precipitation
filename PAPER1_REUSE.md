# Paper 1 复用边界与 Phase 0 入口

日期：2026-11-02。机器可读版本与哈希以 `paper1_sources.lock.json` 为准。

## 可复用来源

| 来源（相对 Paper 1 仓库） | 用途 | 边界 |
|---|---|---|
| `data/parameter_database/releases/production_v2/` | 冻结 LFS handoff，原始输入、LSA 标签、模式与既有 DuMuX onset 表 | 2026-09-13 交付，不能视为全部后续运行的完整库 |
| `paper/1/analysis_manifest.json` | 数据文件哈希、总体筛选、协议及配对分析修订 | 历史版本字段保持原样；不同总体不得混用 |
| `paper/1/protocol_contract.md` | Paper 1 方法约定 | 与实际案例 input、manifest 逐一核对 |
| `paper/1/analysis/wp1_paired_v2_raw/` | 后续成对运行的小型诊断与参数记录 | 84 次运行/48 个物理案例；不等于报告总运行库 |
| `paper/1/analysis/wp1_paired_*.csv`、`wp1_paired_results.json` | 已修正的检测类别、成对可检测性和条件时间比较 | 主比较 36 对，时间比较仅双方 resolved 的 27 对 |
| `paper/1/analysis/review_sharpness_window_20261027/` | sharpness 评价窗口审计 | 保留越界问题，不能默认为物理有效的后续输运区间 |
| `code/dumux_kiemle2023b_unsat/CHANGES.md` 及其补丁 | 非饱和模型实现来源 | 新增守恒通量输出前须隔离源码，记录实际编译环境与补丁 |
| `code/parameter_database/`、`code/lsa_unsaturated*` | 映射、LSA、参数与抽样方法 | 先确认协议和接口，再复用；不依赖不断变化的父目录隐式导入 |

外部合作稿仓库：
`D:/3OutputFile/1Research/Article/22-1-surrogate_model/69e50893f5d9b55adfb56494`。
检查时 HEAD 为 `dd824646d5ec6d2d4b0384f043697ea862e9f936`。

- `AWR_onstime/onstime.tex`：合作稿源文件。
- `Reports/report2.tex`：Option A 所依托报告的源文件。
- 本地 `Reports/report2.pdf` 未被 Git 跟踪；不能仅以这个 PDF 文件名作为可复现版本。
- 本次仅记录上述来源位置和源文件哈希，没有修改合作稿或重新核算报告统计。

## 输入冻结方法

1. `paper1_sources.lock.json` 固定 Paper 1 提交号、数据归档哈希、内部文件哈希
   和选定代码/方法文件哈希，外部 Overleaf 提交另记。
2. `scripts/prepare_paper1_inputs.py` 校验本地输入后，只读取并复制 handoff
   白名单文件到被 Git 忽略的本地目录；不解压整库、不执行模拟、不写原始数据。
3. 父仓库后续提交不自动成为 Paper 2 的基线。要升级输入，先记录原因及影响，
   更新锁文件，并为新版本保留独立数据路径。
4. 需要重新取得数据时，在 Paper 1 仓库执行 release README 中的 Git LFS
   获取命令，再运行本仓库脚本。不以人工复制的未知文件替换缺失输入。

## 尚需完成的 Phase 0 审计

- [ ] 从报告统计追溯到唯一 run ID / physical case ID，核对 436/253/293 的单位及筛选。
- [ ] 清点本地与远程二维场、速度、饱和度、离散控制面盐通量和时间覆盖。
- [ ] 区分正常 horizon、达到溶解度阈值及数值失败造成的停止；保留重试来源。
- [ ] 核对 1D/2D 一致性及守恒输出能力，再选 Pilot-0 的明显输运与竞争边界案例。

本次文件校验通过仅说明归档身份与 Paper 1 清单一致，不证明方案可行、实验验证
成立或某组已有运行足以支持盐输运结论。不得混同 `t_crit`、`t_seed`、`t_detect`
及待定义的 `t_eff`，也不得把盐摩尔分数直接作为质量或累计通量。
