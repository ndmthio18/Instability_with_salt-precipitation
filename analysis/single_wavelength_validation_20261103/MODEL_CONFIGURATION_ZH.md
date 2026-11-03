# 模型配置与执行来源

原 full-domain 数值参照为 Case 36877，600×200、W=0.6000000000000001 m、H=0.2 m、nx=28；R24/R48 分别24×200、48×200，Wr=W/28=0.021428571428571432 m、nx=1。物理波数均293.2153143350473 m⁻¹。原垂向分级、物性、底部压力、蒸发、浓度阈值、幅值1e-6和本征剖面不变。M0为Constant，M1为Eigenprofile。

复用 `model/transport_v1`，未修改或编译C++。ELF SHA256为 `d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc`。输入与剖面SHA逐条存于 `runs/*/plan.json`；原始参照和源码固定于 `reference_hashes.json`。原始600×200结果未覆盖。

所有新运行从零启动，dt=227.3717729539139 s；短窗300步/68211.5318861742 s，完整窗600步/136423.0637723484 s。输出间隔由50改为1以保存每接受步场，不改变积分设置。完整参照仅13帧可重算加权场诊断，短窗共同7帧；守恒预算和原始表面记录均逐步比较。时间仅校正PVD有限有效位到原接受时刻，匹配容差1e-7 s，不进行时间插值。不同网格节点不直接相减。

## 执行修订和失败保留

预登记 `preregistration.json` 保留运行前原文，其中solver_executed=false及sal是当时计划，不代表最终状态。sal无空闲物理核，数组19194在未启动任何求解前取消；证据见 `sal_cancelled_before_start.json`。成功运行迁移至frodo，2条独立任务各一物理核（0/1，不重叠），BLAS线程1，每条12h上限。

首次frodo r01四条启动在MPI初始化处失败，无接受状态；原因是运行库包含orted.real却缺路径适配的orted包装器。保持原ELF/物理输入，添加包装器后在独立r02重新执行，失败目录完整保留。一次r02文件暂存语法错误使管理器文件缺失，没有启动求解器；相应日志/claim另存，未覆盖。见 `launcher_failure_diagnosis.json`、`staging_failure_retained.json`、`failed_initial_launcher/`。这些是启动环境故障，不能混计为数值不收敛或删除。

成功目录 `/temp/zhixin/paper2_wp5_20261103_r02`。先完成4条短窗，登记2条同机全宽10步计时控制，再通过精度、加密和性能门槛后，仅执行R24完整窗M0/M1。总8条成功运行、2420接受步、2428帧场；原始压缩包回传至本地忽略目录 `tmp/single_wavelength_validation_20261103`，全部成员SHA已验证，远端保留。没有2T/3T、幅值修改或其他物理案例。

## 统计约定

独立 `single_wavelength_post.py` 用横向梯形长度权重求均值、方差和标准差；Fourier幅值为2/W乘去均值场与exp(-ikx)的加权积分绝对值。本征投影为去横向均值浓度在A(z)cos(kx)上的体积加权投影，除以该基函数平方积分，不额外乘Sw。盐库存为单位宽度∫phi Sw x dV，水库存为∫phi Sw dV；盐量单位是归一化长度m，并非kg或mol。

通量来自求解器实际守恒SCVF核，未将VTK速度代用为守恒通量。逐步记录局部最大值、浓度阈值和接受步预算。原始sigmatop及冻结onset提取器保持不变；加权sigma只服务网格一致比较。初始库存按积分比较且配对减去初始库存差。节点算术统计的解析采样偏差见 `initial_sampling_bias.csv`。
