# command/14交付索引

已完成方程级诊断、隔离修订、编译和9个有界回归。数值准备GO；新的连续4T未启动，完整4T科学验收仍PENDING。

- [数值失败分析](NUMERICAL_FAILURE_ANALYSIS_ZH.md)：原失败、舍入敏感性、矩阵/差分尺度证据、回归结果。
- [源码审查与准则](SOURCE_REVIEW_AND_CONTRACT_ZH.md)：公式、变量尺度、守恒约束和证据限制。
- [完整4T条件方案与预算](FULL4T_GO_NOGO_AND_COST_ZH.md)：固定输入、GO/NO-GO、原门槛、预计成本。
- [科学解释边界](SCIENTIFIC_BOUNDARY_ZH.md)：早期零模态、再分配、网格和局部阈值的区别。
- [数值诊断图](Figure_Newton_diagnostics.png)：原失败残差、近邻回放、差分尺度和独立守恒预算；同时提供PDF。

复算：恢复archive_receipt所指的忽略缓存后，运行`validate_wp8_v2.py`、`report_wp8_v2.py`和`verify_wp8_v2_delivery.py`。这些脚本不启动模拟。`wp8_v2_regression.py`仅含有界短测名单，已有目录阻止重复prepare/launch；不提供完整4T自动启动入口。
