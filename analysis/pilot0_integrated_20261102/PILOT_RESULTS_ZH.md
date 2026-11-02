# Pilot-0 物理结果：等待新 M1

状态：**in progress，sal M0/M1 各一核并行运行，尚未产生新物理结果**。
执行信息见 `EXECUTION_HANDOFF_ZH.md`、`sal_migration/registration.json`、`sal_migration/latest_status.json`。
本页将在新 M1 完成、实际初始场及正式配对接口通过后更新，不把历史检测时间或旧场当作新输出。

| 科学问题 | 当前答案 |
|---|---|
| 是否观察到相对背景的实际额外盐输运？ | pending：新 M1 未完成，尚无本对证据 |
| 是否高于当前可确定的数值误差？ | pending：尚无配对信号，时间步和网格项 unknown |
| LSA onset 与有效输运是否存在有物理意义的区别？ | 定义已严格区分；本对实证尚 pending，t_eff 未定义 |
| 是否支持 Option A 继续？ | GO 仅限已授权的小规模基准取证；机制结论 inconclusive |
| 是否需要延时、追加案例或改问题？ | 待实际增长、通量、max(x)及成本评价，不能只根据 LSA 阈值时钟决定 |
| 能否进入 Phase 2？ | 当前 NO-GO：尚缺可复现的配对信号与误差评估 |

Figures A–D 均 pending：没有新 M1 数据，不绘制暗示已有结果的图。
后续必须同时呈现净输运方向、额外总通量、adv/diff 补偿及局部最大浓度，避免将向上供应减弱称为净向下迁移。

`ONSET_TRANSPORT_COMPARISON.csv` 当前只有带状态的LSA/历史参考，新的 t_detect 留空。`transport_timeseries.csv` 当前仅表头、零条测量行，等待正式配对输出；不将 pending 填为零输运。
