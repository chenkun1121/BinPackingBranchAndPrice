# Bin Packing — Branch and Price（Java + CPLEX）

用分支定价（Branch-and-Price）精确求解一维装箱问题（Bin Packing）的最小箱子数。
主问题与定价子问题均用 IBM ILOG CPLEX（22.2）求解。

## 数学模型

### 主问题（集合覆盖，Dantzig-Wolfe 分解后的松弛）

对每个可行装箱模式 s（一组可放进同一箱子的物品），设 0-1 变量 y_s：

```
min  Σ_s y_s
s.t. Σ_{s ∋ i} y_s ≥ 1     ∀ 物品 i        （对偶变量 π_i ≥ 0）
     y_s ≥ 0
```

列（模式）数量是指数级的，故只维护受限列集合 S，其余列由定价子问题按需生成。

### 定价子问题（带冲突约束的 0-1 背包）

列 s 的 reduced cost = `1 − Σ_{i∈s} π_i`。寻找负 reduced cost 列即求解：

```
max  Σ_i π_i x_i
s.t. Σ_i w_i x_i ≤ C
     x_i + x_j ≤ 1     ∀ 分箱分支对 (i,j)
     x_i ∈ {0,1}
```

最优值 > 1 ⟺ 存在改进列；否则列生成收敛，当前 RMP 的 LP 值即节点下界。

### 分支规则（Ryan-Foster）

不能直接对分数变量 y_s 分支（会破坏定价子问题的背包结构）。
改为对物品对 (i,j) 分支：

- **together**：i、j 必须同箱 → 合并为一个"有效物品"（`EffItem`）；
- **separate**：i、j 必须分箱 → 定价问题加冲突约束 `x_i + x_j ≤ 1`。

定理：若 LP 解分数，则必存在物品对 (i,j) 使 `f_ij = Σ_{s ⊇ {i,j}} y_s ∈ (0,1)`，
在该对上分支可保证分支树有限且不遗漏最优解。

## 代码结构

| 文件 | 职责 |
|---|---|
| `Instance.java` | 问题实例（物品重量、容量） |
| `EffItem.java` | 有效物品（together 分支合并后的物品组） |
| `MasterProblem.java` | 受限主问题 RMP（CPLEX LP，动态加列，取对偶） |
| `PricingProblem.java` | 定价子问题（CPLEX 背包 MILP，含冲突约束） |
| `BranchAndPrice.java` | 搜索框架：列生成 + 剪枝 + Ryan-Foster 分支 + FFD 初始上界 |
| `Main.java` | 演示算例入口 |

## 运行

Maven（依赖 `com.ibm:cplex:22.2.0.0` 需已在本地仓库）：

```
mvn compile exec:java -Dexec.mainClass=Main
```

或直接用 javac / java（替换为你的 cplex.jar 路径）：

```
javac -encoding UTF-8 -cp <cplex.jar> -d target/classes src/main/java/*.java
java -cp "target/classes;<cplex.jar>" Main
```

## 演示算例结果

容量 10，物品 `{3×4, 4×2, 8×5}`（总重 60）：

```
连续下界 ceil(总重/容量) = 6
FFD 贪心上界   : 8
根节点 LP 下界 : 7.0
最优箱子数     : 7      （B&P 搜索 3 个节点）
```

最优解：五个 8 各自独占一箱，剩余 `{3×4, 4×2}` 装成 `{4,3,3}` + `{4,3}` 两箱。
FFD 把两个 4 捆绑导致 3 无法跟上，多用了 1 箱 —— 体现列生成全局寻模的优越性。
