# 成本与来源

| 阶段 | 接受步 | 末时刻/T | 运行+在线检查s | 含归档s | 压缩MiB | 原场MiB | 求解器峰值RSS MiB |
|---|---|---|---|---|---|---|---|
| T | 1200 | 1.0000000 | 43.663 | 47.588 | 13.586 | 61.214 | 31.332 |
| 4T | 4413 | 3.6775000 | 156.705 | 171.292 | 51.704 | 229.202 | 31.637 |

两阶段顺序各占一核，含在线检查和压缩归档总墙钟约218.880 s，合计约0.060800核时；不含本地分析和网络传输等待。没有重试成本，也没有超时。4T是数值失败，不能因为成本低于预算改写其状态。

每次900 s墙钟上限、每端250 MiB**压缩原场归档**上限均满足；此上限并非展开原场大小上限。各阶段原场仍在frodo，压缩包已回传本地忽略的tmp/wp8_strict_background_20261103并逐成员核对SHA，1201+4414=5615帧全部核验。求解器实施2 GiB地址空间硬上限；表中是求解器子进程峰值RSS，监测Python进程峰值内存未单独记录，不能将未测值填0。

源和ELF未变；stage plan保留command/13与原protocol SHA、源码/可执行文件/input SHA、资源快照，以及4T启动所依赖的T gate SHA。原WP8及更早源和结果按旧provenance再核查；Paper1冻结handoff核验通过，未改Paper1。

复现后处理：先按archive_receipt与SHA恢复两个阶段的忽略缓存，然后运行validate_wp8_strict.py、analyze_wp8_strict.py、diagnose_wp8_strict.py、plot_wp8_strict.py、report_wp8_strict.py、verify_wp8_strict_delivery.py。4T验收预期返回失败；这是被保存的科学结果，不应为了让脚本退出0重跑模拟。已有run目录和launch记录会阻止重复提交。
