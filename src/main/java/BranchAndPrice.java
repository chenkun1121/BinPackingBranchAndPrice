import ilog.concert.IloException;

import java.util.*;

/**
 * 分支定价（Branch-and-Price）精确求解装箱问题。
 *
 * 框架：
 *   1. FFD 贪心给出初始上界（incumbent）；
 *   2. 每个搜索节点上运行列生成（RMP ↔ Pricing 交替）得到节点 LP 下界；
 *   3. 下界 ≥ 当前上界 → 剪枝；LP 解已整 → 更新上界；
 *   4. 否则按 Ryan-Foster 规则在物品对 (i,j) 上分支：
 *        together：i、j 必须同箱 → 合并为一个有效物品（物品数减一）；
 *        separate：i、j 必须分箱 → 定价问题加冲突约束 x_i + x_j ≤ 1。
 *      该规则保证分支后定价子问题仍保持背包结构（这是不能直接对 y_s 分支的原因），
 *      分支树有限，算法精确。
 */
public class BranchAndPrice {

    private static final double EPS = 1e-6;

    private final double capacity;
    private double incumbent = Double.MAX_VALUE;   // 当前最优箱子数（上界）
    private List<Set<Integer>> bestBins;           // 当前最优装箱方案
    private double rootLpBound = -1;
    private int nodeCount = 0;

    public BranchAndPrice(Instance inst) {
        this.capacity = inst.capacity;
    }

    /** 求解结果。 */
    public static class Result {
        public final int bins;                     // 最优箱子数
        public final List<Set<Integer>> solution;  // 每个箱子包含的原始物品编号
        public final int nodes;                    // 搜索节点数
        public final double rootLpBound;           // 根节点 LP 下界
        public final int ffdBound;                 // FFD 贪心上界

        Result(int bins, List<Set<Integer>> solution, int nodes,
               double rootLpBound, int ffdBound) {
            this.bins = bins;
            this.solution = solution;
            this.nodes = nodes;
            this.rootLpBound = rootLpBound;
            this.ffdBound = ffdBound;
        }
    }

    public Result solve(Instance inst) throws IloException {
        // ---- 初始上界：FFD 贪心 ----
        List<Set<Integer>> ffd = firstFitDecreasing(inst);
        incumbent = ffd.size();
        bestBins = ffd;

        // ---- 根节点：每个原始物品自成一个有效物品，无冲突约束 ----
        List<EffItem> rootItems = new ArrayList<>();
        for (int i = 0; i < inst.n(); i++) {
            rootItems.add(new EffItem(inst.weights[i], new TreeSet<>(Set.of(i))));
        }
        dfs(rootItems, Collections.emptyList(), 0);

        return new Result((int) incumbent, bestBins, nodeCount, rootLpBound, ffd.size());
    }

    /** 深度优先搜索一个节点。 */
    private void dfs(List<EffItem> items, List<int[]> forbidden, int depth) throws IloException {
        nodeCount++;

        // ========== 阶段 1：列生成，求节点 LP 下界 ==========
        double lpObj;
        double[] y;
        List<int[]> cols;
        try (MasterProblem master = new MasterProblem(items.size())) {
            for (int e = 0; e < items.size(); e++) {         // 初始列：单物品模式
                master.addColumn(new int[]{e});
            }
            //添加新列直到没有负 reduced cost 列为止
            while (true) {
                double[] pi = master.solve();
                if (pi == null) return;                      // RMP 不可行 → 剪枝
                int[] pattern = PricingProblem.solve(items, forbidden, pi, capacity);
                if (pattern == null) break;                  // 无负 reduced cost → 收敛
                master.addColumn(pattern);
            }
            lpObj = master.objectiveValue();
            y = master.columnValues();
            cols = new ArrayList<>(master.getColumns());
            //System.out.println("LP objective: " + lpObj);
        }
        if (depth == 0) rootLpBound = lpObj;

        // ========== 阶段 2：剪枝（LP 下界取整后不小于当前上界） ==========
        if (Math.ceil(lpObj - 1e-4) >= incumbent - EPS) return;

        // ========== 阶段 3：整数性判定 / Ryan-Foster 选分支对 ==========
        int[] pair = chooseBranchPair(y, cols);
        if (pair == null) {                                  // LP 解已是整数解
            updateIncumbent(items, cols, y);
            return;
        }

        // ========== 阶段 4：分支 ==========
        // 4a. together：i、j 必须同箱 → 合并为一个有效物品
        mergeAndDfs(items, forbidden, pair[0], pair[1], depth);

        // 4b. separate：i、j 必须分箱 → 定价问题加冲突约束
        List<int[]> sepForbidden = new ArrayList<>(forbidden);
        sepForbidden.add(new int[]{pair[0], pair[1]});
        dfs(items, sepForbidden, depth + 1);
    }

    /**
     * Ryan-Foster 分支对选择：
     * 对每个分数列 s（0 < y_s < 1）中出现的物品对 (i,j)，累加 f_ij = Σ_{s⊇{i,j}} y_s。
     * 若存在 0 < f_ij < 1 则选其分支；若所有 f_ij 均为 0/1，则 LP 解必为整数解
     * （Ryan-Foster 定理）。选取 f_ij 最接近 0.5 的对以平衡分支树。
     */
    private int[] chooseBranchPair(double[] y, List<int[]> cols) {
        Map<Long, Double> f = new HashMap<>();
        for (int s = 0; s < cols.size(); s++) {
            if (y[s] <= EPS || y[s] >= 1.0 - EPS) continue;  // 只看分数列
            //cols.get(s) 是一个模式 s 中包含的有效物品下标数组,也就是一个箱子里放了哪些有效物品
            int[] col = cols.get(s);
            //对于每一个模式 s 中的有效物品对 (i,j)，累加 f_ij = Σ_{s⊇{i,j}} y_s
            for (int a = 0; a < col.length; a++) {
                for (int b = a + 1; b < col.length; b++) {
                    f.merge(pairKey(col[a], col[b]), y[s], Double::sum);
                }
            }
        }
        int[] best = null;
        double bestScore = -1.0;
        // 遍历所有物品对的 f_ij 值，选择最接近 0.5 的对作为分支
        for (Map.Entry<Long, Double> e : f.entrySet()) {
            double v = e.getValue();
            if (v > EPS && v < 1.0 - EPS) {
                double score = Math.min(v, 1.0 - v);         // 越接近 0.5 分数越高
                //更新最优分支对
                if (score > bestScore) {
                    bestScore = score;
                    best = unkey(e.getKey());
                }
            }
        }
        return best;
    }

    /**
     * together 分支：把 i、j 合并为一个有效物品并递归；矛盾/超容量则跳过。
     * 左分支：i、j 必须同箱 → 合并为一个有效物品（物品数减一）；
     * 右分支：i、j 必须分箱 → 在定价问题中加入冲突约束。
     * 在合并时，原始物品集合 origItems 也要合并，冲突对下标也要重映射。
     * */
    private void mergeAndDfs(List<EffItem> items, List<int[]> forbidden,
                             int i, int j, int depth) throws IloException {
        // i、j 已被禁止同箱 → 该分支与已有约束矛盾
        for (int[] p : forbidden) {
            if ((p[0] == i && p[1] == j) || (p[0] == j && p[1] == i)) return;
        }
        double mergedWeight = items.get(i).weight + items.get(j).weight;
        if (mergedWeight > capacity + EPS) return;           // 超容量 → 子节点不可行



        // 合并 i 和 j 的原始物品集合
        Set<Integer> orig = new TreeSet<>(items.get(i).origItems);
        orig.addAll(items.get(j).origItems);

        // 构造新的有效物品列表（i、j 移除，合并物品追加在末尾），并重映射冲突对下标
        List<EffItem> newItems = new ArrayList<>();
        // 用于重映射物品下标，因为合并后物品数量减少了，原来的下标可能不再对应新的列表
        Map<Integer, Integer> remap = new HashMap<>();
        for (int k = 0; k < items.size(); k++) {
            if (k == i || k == j) continue;
            remap.put(k, newItems.size());
            newItems.add(items.get(k));
        }
        // 添加合并后的物品
        int mIdx = newItems.size();
        newItems.add(new EffItem(mergedWeight, orig));
        // 重映射冲突对下标，去掉 i、j，添加合并物品
        Set<Long> keys = new HashSet<>();
        for (int[] p : forbidden) {
            int a = (p[0] == i || p[0] == j) ? mIdx : remap.get(p[0]);
            int b = (p[1] == i || p[1] == j) ? mIdx : remap.get(p[1]);
            if (a != b) keys.add(pairKey(a, b));
        }
        List<int[]> newForbidden = new ArrayList<>();
        for (long k : keys) newForbidden.add(unkey(k));

        dfs(newItems, newForbidden, depth + 1);
    }

    /** LP 解整时把它转成可行装箱方案并尝试更新上界。 */
    private void updateIncumbent(List<EffItem> items, List<int[]> cols, double[] y) {
        List<Set<Integer>> bins = new ArrayList<>();
        //获取每个列中的有效物品对应的原始物品集合，构成箱子
        for (int s = 0; s < cols.size(); s++) {
            if (y[s] < 1.0 - EPS) continue;
            Set<Integer> bin = new TreeSet<>();
            for (int e : cols.get(s)) bin.addAll(items.get(e).origItems);
            bins.add(bin);
        }
        // 去重：同一物品出现在多个箱子时归入第一个箱子，再删掉空箱
        Set<Integer> seen = new HashSet<>();
        List<Set<Integer>> cleaned = new ArrayList<>();
        // 遍历所有箱子，去重并删除空箱
        for (Set<Integer> bin : bins) {
            Set<Integer> nb = new TreeSet<>();
            for (int it : bin) {
                if (seen.add(it)) nb.add(it);
            }
            if (!nb.isEmpty()) cleaned.add(nb);
        }
        if (cleaned.size() < incumbent) {
            incumbent = cleaned.size();
            bestBins = cleaned;
           // System.out.println("New incumbent: " + incumbent + " bins, depth=" + nodeCount);
        }
    }

    /**
     * 这个方法用于将两个整数 a 和 b 组合成一个唯一的长整型键值，
     * 确保无论 a 和 b 的顺序如何，生成的键值都是相同的。它通过将较小的整数乘以 1,000,000 并加上较大的整数来实现这一点，从而保证了唯一性。
     */
    private static long pairKey(int a, int b) {
        return (long) Math.min(a, b) * 1_000_000L + Math.max(a, b);
    }
    /**
     * 这个方法用于将长整型键值解码回原始的两个整数 a 和 b。
     * 它通过将键值除以 1,000,000 来获取较小的整数，并使用取模运算获取较大的整数，从而恢复原始的整数对。
     */
    private static int[] unkey(long key) {
        return new int[]{(int) (key / 1_000_000L), (int) (key % 1_000_000L)};
    }

    /** First-Fit-Decreasing 贪心：给出初始可行上界。 */
    public static List<Set<Integer>> firstFitDecreasing(Instance inst) {
        Integer[] order = new Integer[inst.n()];
        for (int i = 0; i < order.length; i++) order[i] = i;
        // 按重量降序排列物品，排的是索引，便于后续装箱时直接使用原始物品编号
        Arrays.sort(order, (a, b) -> Double.compare(inst.weights[b], inst.weights[a]));

        List<Set<Integer>> bins = new ArrayList<>();
        List<Double> residual = new ArrayList<>();
        //从大到小依次尝试放入已有箱子，若都放不下则新开一个箱子
        for (int idx : order) {
            double w = inst.weights[idx];
            boolean placed = false;
            //遍历已有箱子，尝试放入
            for (int b = 0; b < bins.size(); b++) {
                //能放下就放进去，更新剩余容量
                if (residual.get(b) + EPS >= w) {
                    bins.get(b).add(idx);
                    residual.set(b, residual.get(b) - w);
                    placed = true;
                    break;
                }
            }
            // 若所有已有箱子都放不下，则新开一个箱子
            if (!placed) {
                Set<Integer> nb = new TreeSet<>();
                nb.add(idx);
                bins.add(nb);
                residual.add(inst.capacity - w);
            }
        }
        //返回最终的装箱方案，每个箱子包含的原始物品编号，箱子内部物品标号
        return bins;
    }
}
