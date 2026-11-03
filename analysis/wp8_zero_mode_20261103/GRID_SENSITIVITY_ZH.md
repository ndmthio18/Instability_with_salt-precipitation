# R24/R48可观测量审查

| 量 | 窗口 | 分类 | 依据 |
|---|---|---|---|
| canonical增长 | T–3T | Numerically consistent | RMS差0.44–0.91%；仅两档横向网格 |
| canonical增长 | 4T | Grid sensitive | RMS差8.56%；表面与全深投影需各自解释 |
| kE/谱 | T–3T | Numerically consistent | 3T kE差0.00983 rad/m |
| kE/谱 | 4T | Grid sensitive | kE差683.610 rad/m；总能量差约218倍 |
| 均值A | T–3T | Numerically consistent | 3T绝对网格差2.29e-11；小响应非物理有效性 |
| 均值A | 4T | Grid sensitive | 绝对差5.64001e-4；比值约382 |
| M0−候选背景 B | ≤T共有窗 | Below numerical resolution | 最大3.93e-15，且背景严格均匀性失败 |
| 区域库存/累计输运A | 4T约1–5mm | Grid sensitive | 守恒闭合，但网格差远大于可信收敛余量 |
| 区域库存/累计输运A | 4T约10–100mm | Below numerical resolution | 当前dt/z/域宽误差未知，不能用机器闭合当分辨率 |
| 表面均值/最大值 | T–3T | Numerically consistent | 原值接近；不意味着小差量已物理验证 |
| 表面均值/最大值 | 4T | Grid sensitive | M1均值差5.65234e-4；最大值差0.00614977 |
| 相对严格无对流背景的响应 | >T | Unidentifiable because the reference state is unavailable | 缺失长期严格背景；不能用M0替代 |

这些是证据解释分类，不是事后拟合的自动PASS阈值。“Numerically consistent”仅限当前两档横向网格和指定窗口，绝不等于dt、垂向和域宽收敛。“Below numerical resolution”区分机器噪声与在现有验证范围内不能支持物理效应，两者均不填零；原值全部保留。

答案是否定的：晚期水平平均响应并没有比非零谱更稳定。绝对浓度均值的相对误差会被大背景掩盖，应检查A差量的绝对网格差。区域原库存相近同样不能证明新增库存变化稳健。

grid_checkpoints.csv给出原模式、非零能量、物理k、表面均值/最大值、平均响应的T–4T数值和绝对差。2T的既有dt/2证据不能外推为4T晚期误差界；U_dt(晚期)、U_z及U_width均unknown。仅R24/R48两级也不足以作可靠Richardson外推。
