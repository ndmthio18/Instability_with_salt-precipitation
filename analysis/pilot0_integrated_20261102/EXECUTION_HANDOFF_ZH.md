# command/3.md 基准工作包完成交接

状态：**M0/M1 completed、paired analysis completed；transport conclusion inconclusive**。
两条sal任务19160/19161均完成600步并达到136423.0637723484 s；没有正在运行的本轮任务，也没有自动重试或新延窗任务。

## 最终证据

- 正式同主机配对通过；实际初始场、库存修正、完整接受步、输出SHA及M0跨主机桥接通过。
- 26帧原始VTU和PVD已归档到本地忽略目录tmp/sal_migration/M0_fields.tar.gz、M1_fields.tar.gz；sal /tmp/zhixin_paper2_pilot0_20261102_r01原文件保留。
- sal_migration/field_archive.json记录压缩包及每个成员SHA；baseline_analysis.json记录完成/完整性及科学计算摘要。
- 六项交付、1803行三深度配对时序、26帧共同模态投影及四类PNG/PDF图已完成。figure_qa.json记录图与数据检查。
- 新M1独立检测约50394.334 s；额外输运差值约1e-17 m，闭合残差量级相当；未知敏感性误差保持unknown，t_eff undefined。

## 条件分支与停止点

Task5：因弱扰动和未分辨差值，不机械执行dt/2/Nz/Nx加密，误差评估未完成而非通过。
Task6：2T配对已形成成本化候选；仅预计约1.4–1.54倍当前sigma，约13.04核时。本轮不执行，4T的M0成本超过12 h限制。restart未验证。
Task7：38852采用来源核验的WP1输入/profile形成匹配草案；旧密度冲突场未使用。未达到执行门槛，不启动。
Task8–9：完成有证据边界的机制解释和下一阶段决策，Phase2 NO-GO。

本轮已无待监测计算；paper-2-pilot-0自动轮询已删除，聊天保留。未来继续工作应从NEXT_STAGE_DECISION_ZH.md读起，不重提旧M0/M1，不重启已取消frodo worker2453288。

## 重现

如需补回本地忽略的原始归档，在确认sal目录仍存在后使用以下命令（不会启动求解器）：

```powershell
python scripts/pilot0_sal_campaign.py collect
python scripts/pilot0_sal_campaign.py archive
python scripts/analyze_matched_pair.py --m0 analysis/pilot0_integrated_20261102/sal_migration/runs/M0 --m1 analysis/pilot0_integrated_20261102/sal_migration/runs/M1 --out analysis/pilot0_integrated_20261102/sal_baseline_pair
python scripts/analyze_pilot0_baseline.py
python scripts/write_pilot0_report.py
```

provenance/extract_onset_2025.py为Paper1提交158395f3d25a的只读副本，SHA冻结。分析不修改Paper1。case38852_design/source的远端manifest与冻结本地副本字节不同但JSON语义相同，初次哈希拒绝及两者SHA均记录在source_hashes.json。
历史队列/迁移/并行记录、冻结计划和原始manifest保持为历史证据，不因最终状态而改写。
