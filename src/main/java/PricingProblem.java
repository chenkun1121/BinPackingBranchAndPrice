import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * 定价子问题（Pricing Problem）—— 带冲突约束的 0-1 背包：
 *
 *   max  Σ_e π_e x_e
 *   s.t. Σ_e w_e x_e ≤ C
 *        x_i + x_j ≤ 1    ∀ (i,j) ∈ Forbidden    （Ryan-Foster"分箱"分支）
 *        x_e ∈ {0,1}
 *
 * 任一可行模式的 reduced cost = 1 − Σ_{e∈s} π_e。
 * 当且仅当背包最优值 > 1 时，存在 reduced cost < 0 的改进列，返回该新模式；
 * 否则列生成收敛，当前 RMP 的 LP 值即节点下界。
 *
 * 【求解策略：DP 优先，CPLEX 兜底】
 * 主路径改为**动态规划**，复杂度 O(n·C)，只需整数重量与整数容量：
 *   - dp[w] = 装载总重恰好为 w 时能取得的最大对偶价值；
 *   - 用 lastBuf[e][w] 记录决策，最后自后向前回溯出具体的物品集合。
 *
 * 两种情况下 DP 不适用，自动回退 CPLEX（保证结果始终是最优解）：
 *   1. **存在冲突约束且 DP 解违反了它** —— x_i + x_j ≤ 1 破坏了 DP 的最优子结构，
 *      无法用一维容量的 dp 表达。此时先用 DP 求"无冲突松弛"的最优解：
 *      若该解恰好不违反任何冲突对，它也就是原问题的最优解，直接采用；
 *      否则回退 CPLEX。
 *   2. **规模过大** —— 状态单元数 n×(C+1) 超过 DP_MAX_CELLS 时，O(n·C) 的时间与
 *      内存都不划算（例如 Scholl 算例 C=100000）。
 *
 * 另外重量必须是整数，否则同样回退 CPLEX。
 */
public final class PricingProblem {

    private static final double EPS = 1e-6;

    /**
     * 是否优先使用 DP。置 false 可强制全部走 CPLEX（用于 A/B 对比）。
     * 可用 JVM 参数 -Dbp.pricing.dp=false 覆盖，也可运行时用 setUseDp 切换。
     */
    private static boolean useDp = Boolean.parseBoolean(
            System.getProperty("bp.pricing.dp", "true"));

    /** 运行时切换定价策略（true = 优先 DP，false = 全部走 CPLEX）。 */
    public static void setUseDp(boolean v) {
        useDp = v;
    }

    /** 当前 DP 规模上限，供诊断输出使用。 */
    public static long getDpMaxCells() {
        return DP_MAX_CELLS;
    }

    public static boolean isUseDp() {
        return useDp;
    }

    /** 定价环境复用开关，见 ensureModel 中的实测说明，默认关闭（实测更快）。 */
    private static final boolean REUSE_ENV = false;

    /**
     * DP 规模上限（n × (C+1) 个状态单元）。超过则改用 CPLEX。
     * 可用 JVM 参数 -Dbp.pricing.dp.cells=N 调整。
     *
     * 默认 3e7 是为了覆盖 Scholl 这类大容量算例：n=200、C=100000 时
     * n×(C+1) = 2e7，若上限只有 4e6 就会全程回退 CPLEX，
     * 实测 60 秒连根节点都算不完，而 DP 下根节点约 15 秒。
     */
    private static final long DP_MAX_CELLS = Long.parseLong(
            System.getProperty("bp.pricing.dp.cells", "30000000"));

    /**
     * lastBuf 允许占用的堆内存比例（每个状态单元 2 字节）。
     * 用于兜底，避免大算例直接 OOM。
     */
    private static final double DP_MEMORY_FRACTION = Double.parseDouble(
            System.getProperty("bp.pricing.dp.mem", "0.25"));

    // ---- DP 缓冲区：跨调用复用，避免反复分配大数组 ----
    private static double[] dpBuf;
    /** lastBuf[e][w]：处理完物品 0..e 后，dp[w] 最后一次是被哪个物品更新的（-1 表示不可达）。 */
    private static short[][] lastBuf;

    // ---- 统计 ----
    private static long dpCalls;
    private static long cplexCalls;

    // ---- CPLEX（兜底）----
    private static IloCplex cplex;
    private static IloNumVar[] x;
    private static IloObjective obj;
    private static IloRange capRange;
    private static final List<IloRange> conflictRanges = new ArrayList<>();
    private static int builtN = -1;
    private static String lastStructKey = null;

    private PricingProblem() {}

    /**
     * @param items     当前节点的有效物品
     * @param forbidden 禁止同箱的有效物品下标对
     * @param pi        有效物品层面的对偶价格（长度必须等于 items.size()）
     * @param capacity  箱子容量
     * @return 价值超过 1 的新模式（有效物品下标）；无改进列时返回 null
     */
    public static int[] solve(List<EffItem> items, List<int[]> forbidden,
                              double[] pi, double capacity) throws IloException {
        if (useDp) {
            int[] pat = solveByDp(items, forbidden, pi, capacity);
            if (pat != null) {
                dpCalls++;
                // 空数组 = DP 已判定"无改进列"，可以直接收敛，不必再惊动 CPLEX
                return pat.length == 0 ? null : pat;
            }
            // null = DP 不适用，落到下面的 CPLEX 兜底
        }
        cplexCalls++;
        return solveByCplex(items, forbidden, pi, capacity);
    }

    // ==================================================================
    //                        DP（主路径）
    // ==================================================================

    /** DP 判定"无改进列"时返回此哨兵（区别于 null 的"不适用"）。 */
    private static final int[] NO_IMPROVING_COLUMN = new int[0];

    /**
     * 动态规划求解 0-1 背包。
     *
     * 返回值有三种语义，务必区分：
     *   - 非空数组：找到改进列；
     *   - {@link #NO_IMPROVING_COLUMN}（空数组）：DP **已确定**不存在改进列，
     *     调用方应直接收敛，不要再调用 CPLEX；
     *   - null：DP 不适用（规模过大 / 重量非整数 / 解违反冲突约束），调用方回退 CPLEX。
     */
    private static int[] solveByDp(List<EffItem> items, List<int[]> forbidden,
                                   double[] pi, double capacity) {
        int n = items.size();
        int C = (int) capacity;
        if (capacity != C || C <= 0) return null;                 // 容量非整数
        long cells = (long) n * (C + 1L);
        if (cells > DP_MAX_CELLS) return null;                    // 规模过大
        // lastBuf 用 short 存物品下标，每个状态单元 2 字节；超出堆预算则回退 CPLEX
        long budget = (long) (Runtime.getRuntime().maxMemory() * DP_MEMORY_FRACTION);
        if (cells * 2L > budget) return null;                     // 内存不够

        int[] wt = new int[n];
        for (int e = 0; e < n; e++) {
            double w = items.get(e).weight;
            if (w != Math.floor(w) || w < 0 || w > C) return null; // 非整数重量或单件装不下
            wt[e] = (int) w;
        }

        ensureDpBuffers(n, C);

        // ---- DP：dp[w] = 总重恰为 w 的最大对偶价值，-1 表示该重量不可达 ----
        Arrays.fill(dpBuf, 0, C + 1, -1.0);
        dpBuf[0] = 0.0;
        for (int e = 0; e < n; e++) {
            if (e == 0) {
                Arrays.fill(lastBuf[0], 0, C + 1, (short) -1);
            } else {
                System.arraycopy(lastBuf[e - 1], 0, lastBuf[e], 0, C + 1); // 继承上一阶段
            }
            int we = wt[e];
            double ve = pi[e];
            for (int w = C; w >= we; w--) {          // 倒序保证每件物品只取一次
                double prev = dpBuf[w - we];
                if (prev < 0) continue;              // w-we 不可达
                double cand = prev + ve;
                if (cand > dpBuf[w]) {
                    dpBuf[w] = cand;
                    lastBuf[e][w] = (short) e;       // 标记为"在阶段 e 被更新"
                }
            }
        }

        // ---- 取最优（总重不超过 C 的前提下价值最大）----
        int bestW = 0;
        double best = dpBuf[0];
        for (int w = 1; w <= C; w++) {
            if (dpBuf[w] > best) {
                best = dpBuf[w];
                bestW = w;
            }
        }
        if (best <= 1.0 + EPS) return NO_IMPROVING_COLUMN;   // DP 已确定无改进列

        // ---- 回溯出物品集合 ----
        int[] pat = new int[n];
        int m = 0;
        int w = bestW;
        for (int e = n - 1; e >= 0; e--) {
            if (w >= wt[e] && lastBuf[e][w] == e) {
                pat[m++] = e;
                w -= wt[e];
            }
        }
        if (w != 0) return null;                     // 回溯异常，交给 CPLEX
        int[] result = Arrays.copyOf(pat, m);
        Arrays.sort(result);

        // ---- 冲突约束：DP 求的是无冲突松弛，需检验可行性 ----
        if (forbidden != null && !forbidden.isEmpty()) {
            BitSet b = new BitSet(n);
            for (int e : result) b.set(e);
            for (int[] p : forbidden) {
                if (b.get(p[0]) && b.get(p[1])) return null;   // 违反，回退 CPLEX
            }
        }
        return result;
    }

    /** 按需分配并复用 DP 缓冲区。 */
    private static void ensureDpBuffers(int n, int C) {
        if (dpBuf == null || dpBuf.length < C + 1) {
            dpBuf = new double[C + 1];
        }
        if (lastBuf == null || lastBuf.length < n) {
            lastBuf = new short[n][];
        }
        for (int e = 0; e < n; e++) {
            if (lastBuf[e] == null || lastBuf[e].length < C + 1) {
                lastBuf[e] = new short[C + 1];
            }
        }
    }

    // ==================================================================
    //                        CPLEX（兜底路径）
    // ==================================================================

    private static int[] solveByCplex(List<EffItem> items, List<int[]> forbidden,
                                      double[] pi, double capacity) throws IloException {
        int n = items.size();
        double[] weights = items.stream().mapToDouble(e -> e.weight).toArray();

        // 新建环境时直接用当前数据一次性装配，省掉后续所有改系数的开销；
        // 只有复用环境（REUSE_ENV=true）时才需要增量修改。
        if (!ensureModel(n, capacity, weights, pi, forbidden)) {
            for (int e = 0; e < n; e++) {
                cplex.setLinearCoef(obj, x[e], pi[e]);
            }
            String key = structKey(weights, forbidden);
            if (!key.equals(lastStructKey)) {
                for (int e = 0; e < n; e++) {
                    cplex.setLinearCoef(capRange, x[e], weights[e]);
                }
                rebuildConflicts(forbidden);
                lastStructKey = key;
            }
        }

        if (!cplex.solve() || cplex.getStatus() != IloCplex.Status.Optimal) {
            return null;
        }
        if (cplex.getObjValue() <= 1.0 + EPS) {       // 无负 reduced cost 列
            return null;
        }

        List<Integer> pattern = new ArrayList<>();
        for (int e = 0; e < n; e++) {
            if (cplex.getValue(x[e]) > 0.5) pattern.add(e);
        }
        int[] result = new int[pattern.size()];
        for (int i = 0; i < result.length; i++) result[i] = pattern.get(i);
        return result;
    }

    /**
     * 确保模型按 n 个变量建好；若是新建的，直接用传入的 pi / weights / forbidden
     * 一次性装配完成（避免之后再逐个改系数，实测逐个 setLinearCoef 会慢数倍）。
     *
     * @return true 表示本次是新建且已装配完毕，调用方无需再改任何系数；
     *         false 表示复用了旧环境，调用方需自行更新目标与容量系数。
     */
    private static boolean ensureModel(int n, double capacity, double[] weights,
                                       double[] pi, List<int[]> forbidden) throws IloException {
        boolean fresh = false;
        if (cplex == null || builtN != n || !REUSE_ENV) {
            if (cplex != null) {
                cplex.end();                          // 变量数变了，只能换一个环境
            }
            cplex = new IloCplex();
            cplex.setOut(null);
            cplex.setWarning(null);
            x = cplex.boolVarArray(n);
            obj = cplex.addMaximize(cplex.scalProd(x, pi));
            capRange = cplex.addLe(cplex.scalProd(x, weights), capacity);
            conflictRanges.clear();                   // 旧 range 已随环境销毁
            for (int[] p : forbidden) {
                conflictRanges.add(cplex.addLe(cplex.sum(x[p[0]], x[p[1]]), 1.0));
            }
            lastStructKey = structKey(weights, forbidden);
            builtN = n;
            fresh = true;
        }
        return fresh;
    }

    /** 拆掉旧的冲突约束并按新的 forbidden 重建。 */
    private static void rebuildConflicts(List<int[]> forbidden) throws IloException {
        for (IloRange r : conflictRanges) {
            cplex.remove(r);
        }
        conflictRanges.clear();
        for (int[] p : forbidden) {
            conflictRanges.add(cplex.addLe(cplex.sum(x[p[0]], x[p[1]]), 1.0));
        }
    }

    /** 节点结构签名：有效物品重量序列 + 冲突对列表。 */
    private static String structKey(double[] weights, List<int[]> forbidden) {
        StringBuilder sb = new StringBuilder();
        for (double w : weights) sb.append(w).append(',');
        sb.append('|');
        for (int[] p : forbidden) sb.append(p[0]).append('-').append(p[1]).append(';');
        return sb.toString();
    }

    // ==================================================================
    //                        统计 / 资源释放
    // ==================================================================

    /** 重置 DP / CPLEX 调用计数。 */
    public static void resetStats() {
        dpCalls = 0;
        cplexCalls = 0;
    }

    /** 返回 "DP 次数 / CPLEX 次数" 的统计串。 */
    public static String stats() {
        long total = dpCalls + cplexCalls;
        double ratio = total == 0 ? 0 : 100.0 * dpCalls / total;
        return String.format("DP=%d  CPLEX=%d  (DP 占比 %.1f%%)", dpCalls, cplexCalls, ratio);
    }

    /** 释放 CPLEX 环境（一次求解结束后调用）。 */
    public static void shutdown() {
        if (cplex != null) {
            cplex.end();
            cplex = null;
            x = null;
            obj = null;
            capRange = null;
            conflictRanges.clear();
            builtN = -1;
            lastStructKey = null;
        }
    }
}
