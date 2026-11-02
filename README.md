# Paper 2 — effective salt redistribution

本目录已独立为一个 Git 仓库。研究方案见
[Option_A_Research_Plan.md](Option_A_Research_Plan.md)，当前任务是 Phase 0
审计与 Pilot-0 准备，尚未开始新的数值试验。

## 仓库边界

| 内容 | 位置与管理方式 |
|---|---|
| Paper 1、公共基线代码和冻结数据归档 | 父级 `22Surrogate_model` 仓库，既有远程与分支保持不变 |
| Paper 2 研究方案、代码、文稿、分析记录 | 本仓库，分支 `codex/paper2-main` |
| Paper 1 数据版本和文件校验值 | `paper1_sources.lock.json` |
| Paper 2 使用的本地外部数据副本 | `data/external/paper1/`，已忽略，不重复提交 |
| Paper 2 进展 | `PROGRESS.md` |
| 正式合作稿和 11 月 1 日报告源文件 | 独立 Overleaf 仓库，定位见 `PAPER1_REUSE.md` |

本仓库不是父仓库的分支或子模块。父仓库忽略整个 `paper/2/`，因此 Paper 2
的修改不会进入 Paper 1 的提交。本仓库已建立本地历史，尚未配置独立远程；
推送 Paper 1 不会备份 Paper 2，之后需为 Paper 2 配置自己的远程仓库。

## 准备 Paper 1 输入数据

在本目录执行（仅 Python 标准库）：

```powershell
python scripts/prepare_paper1_inputs.py --verify-only
python scripts/prepare_paper1_inputs.py
git status
```

默认读取本目录上两级的 Paper 1 仓库。若本仓库移到别处：

```powershell
python scripts/prepare_paper1_inputs.py --paper1-root 'D:/path/to/22Surrogate_model'
```

脚本核验固定提交存在、选定源文件及 LFS 归档的 SHA-256，再从归档中按白名单
读取 handoff 文件，写到以基线提交号区分的数据目录。已有同名文件必须通过
校验才会复用；不覆盖不匹配文件。数据目录还保存独立 provenance 清单。
源码校验允许 Windows CRLF 与 LF 的纯换行差异；归档和 handoff 数据仍进行
逐字节 SHA-256 校验。

数据准备不代表 Phase 0 已完成：原始二维场、离散盐通量、停止原因和 onset
之后的覆盖长度仍需逐运行审计。方案中的 436 次运行、253 个物理案例和 293
条配对观测是待追溯的报告统计，不能由旧 handoff 直接替代。

## 日常工作

```powershell
git -C 'D:/3OutputFile/1Research/Article/22Surrogate_model/paper/2' status
```

新代码与分析写在本仓库。不要在父仓库现有模型上直接试改 Paper 2 功能；
需要改变数值模型时，先建立可追溯的隔离源码副本并记录上游版本及补丁。
不要把新的输运结果写回 Paper 1 的结果目录。
