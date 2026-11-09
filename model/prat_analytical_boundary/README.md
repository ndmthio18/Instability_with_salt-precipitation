# 无盐层网格的有限 Da 生长边界

基线 `fdd4d1f306c09caf3f44cac39f9cc7431aadfaa0`。只生成新的应用；冻结 P1-R2、P2-H、P2-F 源文件均不修改。

## 模型及近似

使用**无内部反应**的界面关系约化式，并非 Prat (2026) 式 (28)：

```
K = beta (1-epsilon) (1-Ceq) kr / (epsilon Deff)
Da_top = K h
v = Ceq j / [rho_cr (1-epsilon) (1-Ceq)] Da_top/(1+Da_top)
Ctop = Ceq + rho_cr v/(beta rho_l kr)
j = gasFilm(Ctop)
```

固定 epsilon=0.43、rho_cr=2165 kg/m3、beta=2、kr=0.0023 m/s、Deff=sqrt(epsilon)*1.3e-9 m2/s。Ceq 由冻结土壤的 xeq=0.1 换算，不把摩尔分数当作质量分数。初始和底部储液槽 x=0.08。二维土壤几何仍为宽0.001 m、高0.20 m的单横向单元柱；盐层没有网格。所有表面库存按单位土壤面积计。

启动 h=0、固体=0、孔液=0。冻结的隐式表面互补条件 `precip>=0, Csoil<=Ceq, precip*(Csoil-Ceq)=0` 计算真实析盐通量，同时按孔隙体积增长渐进填充。首次实际固体高度达到3.5 um时，仅切换状态，继承超出阈值的实际高度和完整库存。不把3.5 um用作网格尺度或人为种子。

激活后每步用已接受 h 计算 Da_top，二分联立气膜与界面动力学（64次），采用显式高度更新；土壤为后向Euler。这是一阶分步耦合，需要独立的 dt/2 检查。新增孔液取 Ceq 和对应土壤密度，原有孔液组成不覆盖。每步生成

```
dSolid = rho_cr*(1-epsilon)*dh
dPoreSalt = rho_l*epsilon*Ceq*dh
dWater = rho_l*epsilon*(1-Ceq)*dh
Fs = (dSolid+dPoreSalt)/dt
Fw = j+dWater/dt
```

DuMux 第一式边界值为 Fs+Fw，第二式为 Fs，均 kg/m2/s、向外正。激活后不再调用预析盐源，也不直接改写土壤库存。外部库存与土壤整步共同接受；拒绝时恢复土壤、变量缓存及边界通量，外部候选库存不提交。最多9次减半试算，持续越界停止。失败不通过裁剪C或增长速率解决。

本闭合**没有**盐层水力阻力、气相扩散阻力、孔隙率演化或内部沉淀。厚度经有限Da动力学和Ctop/水活度影响气膜蒸发；固定孔隙率不能提供堵塞反馈。Sw>=0.99只是供液适用性筛选，不是毛细输液方程。新孔液平衡近似、准稳态浓度关系及连续盐层启动尺度仍需物理验证。

## 有效域与独立预算

同时检查土壤界面和推算盐层顶端相对Ceq偏差<=2%、Sw>=0.99、
`Pe_U = Ftotal*h/(rho_l*epsilon*Deff)<=0.1`、`Pe_move=v*h/Deff<=0.1`。
底部Dirichlet交换从底部BOX库存变化和紧邻内部SCVF通量独立重建，绝不由全域残差反算。组合预算包含土壤水/盐、表面固体、孔液水/盐、底部输入及蒸发输出；每个接受步检查相对残差<=1e-8。激活不搬移土壤质量，库存切换残差<=1e-12。

## 复现

在HPCC的 `/WORK/hulm_work/zhixin/prat_analytical_boundary_20261109` 保留本目录及冻结依赖为 `repo/`，运行：

```
sbatch repo/model/prat_analytical_boundary/build.slurm
sbatch repo/model/prat_analytical_boundary/run.slurm pilot
# 仅当先导 coupling_pass 为 true：
sbatch repo/model/prat_analytical_boundary/run.slurm matrix 0
sbatch repo/model/prat_analytical_boundary/run.slurm matrix 1
sbatch repo/model/prat_analytical_boundary/run.slurm matrix 2
sbatch repo/model/prat_analytical_boundary/run.slurm matrix 3
```

编译脚本核对冻结DUNE/DuMux栈归档SHA256，保存完整编译日志及ELF校验和。不需要SUNDIALS。`run.py` 保存每个实际命令、Slurm job/node、退出码和物理终点。嵌套网格直接读取80层先导的初始节点，40层合并两格，160/320层等分原格；不重新使用固定几何比生成另一套网格。

恒j独立积分参考为 `h-h0+log(h/h0)/K = [Ceq*j/(rho_cr*(1-epsilon)*(1-Ceq))] t`。C++单元测试、Python独立数值实现及积分根同时比较。恒j测试的非零h0是独立边界求积测试初值，生产计算始终从零库存开始。
