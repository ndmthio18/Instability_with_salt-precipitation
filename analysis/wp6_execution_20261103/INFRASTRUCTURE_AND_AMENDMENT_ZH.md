# 预登记保护、独立实现与执行修订

本轮授权来源为 `command/7.md`。原WP6计划提交314d75d9及其preregistration、门槛、预算和产品全部保持原字节；`amendment_20261103_A01.json`记录新的执行授权：两个微型阈值测试、T桥接、桥接通过后的2T配对，然后停止。不将旧计划中的4T/8T设计解释为本轮启动许可。

## 不可变与恢复机制

`prepare_wp6_plan.py` 无参数或 `--verify` 默认只读检查冻结文件和依赖哈希。历史artifact_manifest记录的旧生成脚本由Git提交314d75d9校验，新脚本变更由本次修订记录，不能为了让验证通过而回写原manifest。

重建只允许 `--recreate-to` 指定Paper2/tmp内不存在的新子目录，显式拒绝原计划目录、tmp根目录以及任何已有目录。更改研究协议须另加带日期的amendment，不允许静默覆盖冻结计划。静态测试检查默认核验不改变计划文件mtime，并验证向冻结目录重建会失败。

`wp6_plan_integrity.py --restore-reference` 从既有正式M0/M1原始场归档恢复忽略缓存，先验压缩包SHA，再验各成员SHA，拒绝路径逃逸、非文件成员或与已存在数据冲突。重算只写新缓存目录，生成NPZ后与WP5原field_audit中的SHA比较，不回写任何WP5报告/CSV。已在全新缓存目录做冷恢复：26帧全部验证，重新生成的NPZ哈希与WP5一致。见 `restoration_checks.json`。

## 独立七面实现

`model/transport_wp6_v1/source`以旧transport_v1复制；除transportdiagnostics.hh以外所有源文件逐字节相同，尤其main.cc/model.hh的物理残差、本构、初始化、边界和冻结检测器未改。补丁见 `provenance/diagnostics.patch`；新ELF SHA为 `112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81`。

七个实际BOX面坐标由冻结计划传入，不存在于该网格上即拒绝；原3面名字保持不变，新增4面，free_domain保持最后一项以维持底部反应诊断语义。元数据schema更新为2。所有区域使用同一真实残差通量核，记录adv/diff/total、水、库存、累计量和预算。

加权在线诊断另写weighted_online.csv：表面长度加权均值/sigma、物理k Fourier、体积加权本征投影、扰动RMS及2/3次表面谐波。M0本征投影无定义，写空及projection_defined=0；不伪造为0。逐步VTU后处理沿用固定A28另作一致性核对，初始场、密度、饱和度也检查。

独立检查按实际VTU张量坐标重建BOX积分权重，对8个记录控制体（7观测面+自由域）的盐/水库存重新积分，并重新对每接受步记录通量累积求和。该检查独立于求解器存储求和，但不是第二套独立物理通量实现；两者不能混称。

## 执行与问题记录

frodo独立构建目录 `/temp/zhixin/paper2_wp6_build_20261103`；执行目录 `/temp/zhixin/paper2_wp6_execution_20261103_r01`。沿用WP5库运行环境复制成独立bundle，仅替换新ELF和路径包装器，并重新记录运行库清单。每腿1 rank、独立物理核0/1、BLAS线程1；启动检查内存、磁盘和负载。预检其他高负载计算位于物理核6/7，不与0/1重叠。

一次暂存脚本因字符串换行转义导致远端Python语法错误，在任何语句执行和任何WP6求解启动前失败；修正为raw-string后暂存成功，记录于 `staging_correction.json`。桥接分析首次JSON输出遇到numpy布尔类型，改为原生bool后重做后处理；没有重启求解或改变验收容差。原始失败信息和修订原因保留，不误计为数值失败。

新ELF完成两个独立微型停止测试：低于初始浓度的测试阈值触发第0步停止；另一测试在第1个接受步首次越过低测试阈值即停止，并保留前后区间。测试阈值不是科学配对阈值；全部正式T/2T运行严格使用0.1。

所有正式运行从t=0开始，未使用restart。原始场压缩包存本地忽略目录与独立远端目录；每条完成后核对包与成员SHA。没有修改或重算Paper1/WP5基准。
