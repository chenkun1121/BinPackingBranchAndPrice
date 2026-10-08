import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

/**
 * 受限主问题（Restricted Master Problem，RMP）—— 集合覆盖（Set Covering）松弛：
 *
 *   min  Σ_{s∈S} y_s
 *   s.t. Σ_{s∈S: o∈s} y_s ≥ 1     ∀ 原始物品 o      （对偶变量 π_o ≥ 0）
 *        y_s ≥ 0
 *
 * 【一次建模、跨节点复用】
 * 整个搜索过程只创建 **一个** IloCplex 环境、只建一次目标函数与 N 条覆盖约束
 * （N = 原始物品数，固定不变）。节点之间的差异不再通过"重建模型"表达，而是：
 *
 *   1. 列池全局化 —— 所有节点生成的列都进同一个池，模式一律用**原始物品下标**表示，
 *      因此列的语义在整棵树上一致，可以安全地跨节点复用；
 *   2. 节点切换用变量界 —— {@link #setNode(Predicate)} 把本节点非法的列
 *      变量上界设为 0（屏蔽）、合法的放开为 +∞。变量 UB=0 等价于该列不存在，
 *      CPLEX 预处理会自动剔除，因此无需增删任何约束或变量。
 *
 * 这样每个节点的开销只剩下"改一批 UB"，而不是"重建一个 MILP"。
 */
public class MasterProblem implements AutoCloseable {

    private final IloCplex cplex;
    private final IloObjective objective;
    /** 每个原始物品一条覆盖约束。下标 = 原始物品编号，全树固定。 */
    private final IloRange[] cover;
    /** 全局列池的变量，与 columns 一一对应。 */
    private final List<IloNumVar> vars = new ArrayList<>();
    /** 全局列池：每列覆盖的原始物品下标（已升序排序）。 */
    private final List<int[]> columns = new ArrayList<>();
    /** 列去重：规范化键 -> 池中下标。 */
    private final Map<String, Integer> keyToIndex = new HashMap<>();

    /**
     * @param numOrigItems 原始物品数（注意：不是箱子容量）
     */
    public MasterProblem(int numOrigItems) throws IloException {
        cplex = new IloCplex();
        cplex.setOut(null);                            // 关闭求解日志
        cplex.setWarning(null);
        objective = cplex.addMinimize();
        cover = new IloRange[numOrigItems];
        for (int o = 0; o < numOrigItems; o++) {
            cover[o] = cplex.addRange(1.0, Double.MAX_VALUE); // Σ y_s ≥ 1
        }
    }

    /**
     * 追加一列（一个箱子模式）。
     * @param origPattern 该模式覆盖的**原始物品**下标数组（会被排序，不修改入参）
     * @return 该列在池中的下标；若池中已存在同一模式，直接返回已有下标且不新建变量
     */
    public int addColumn(int[] origPattern) throws IloException {
        int[] sorted = origPattern.clone();
        Arrays.sort(sorted);
        String key = Arrays.toString(sorted);
        Integer existing = keyToIndex.get(key);
        if (existing != null) return existing;

        IloColumn col = cplex.column(objective, 1.0);          // 目标系数 1
        for (int o : sorted) {
            col = col.and(cplex.column(cover[o], 1.0));        // 覆盖约束系数 1
        }
        vars.add(cplex.numVar(col, 0.0, Double.MAX_VALUE, "y" + columns.size()));
        columns.add(sorted);
        int idx = columns.size() - 1;
        keyToIndex.put(key, idx);
        return idx;
    }

    /** 池中是否已存在该模式。 */
    public boolean containsColumn(int[] origPattern) {
        int[] sorted = origPattern.clone();
        Arrays.sort(sorted);
        return keyToIndex.containsKey(Arrays.toString(sorted));
    }

    public int numColumns() {
        return columns.size();
    }

    /** 第 i 列的模式（原始物品下标，升序）。 */
    public int[] getColumn(int i) {
        return columns.get(i);
    }

    /**
     * 切换到某个节点：按该节点的分支约束屏蔽 / 放开列变量。
     *
     * 只修改变量上界，**不动**模型中的任何约束与变量结构 —— 这是"一次建模、
     * 跨节点复用"的技术支点。被屏蔽的列（UB=0）在 LP 中取值恒为 0，
     * 等价于该列不存在，且不会被 {@link #columnValues()} 之后的解读取逻辑采纳。
     *
     * @param isValid 判定某个模式（原始物品下标）在本节点是否合法
     * @return 本节点激活（合法）的列数
     */
    public int setNode(Predicate<int[]> isValid) throws IloException {
        int active = 0;
        for (int i = 0; i < vars.size(); i++) {
            IloNumVar v = vars.get(i);
            if (isValid.test(columns.get(i))) {
                if (v.getUB() <= 0.0) v.setUB(Double.MAX_VALUE);
                active++;
            } else {
                if (v.getUB() > 0.0) v.setUB(0.0);
            }
        }
        return active;
    }

    /** 求解 RMP，返回原始物品覆盖约束的对偶价格 π（不可行返回 null）。 */
    public double[] solve() throws IloException {
        if (!cplex.solve() || cplex.getStatus() != IloCplex.Status.Optimal) {
            return null;
        }
        double[] pi = new double[cover.length];
        for (int o = 0; o < cover.length; o++) {
            pi[o] = cplex.getDual(cover[o]);
        }
        return pi;
    }

    public double objectiveValue() throws IloException {
        return cplex.getObjValue();
    }

    /** 各列当前的 LP 取值（被屏蔽的列恒为 0）。 */
    public double[] columnValues() throws IloException {
        double[] v = new double[vars.size()];
        for (int i = 0; i < vars.size(); i++) {
            v[i] = cplex.getValue(vars.get(i));
        }
        return v;
    }

    /** 当前激活列（UB>0）在列池中的下标。 */
    public int[] activeColumnIndices() throws IloException {
        int[] tmp = new int[vars.size()];
        int k = 0;
        for (int i = 0; i < vars.size(); i++) {
            if (vars.get(i).getUB() > 0.0) tmp[k++] = i;
        }
        return Arrays.copyOf(tmp, k);
    }

    @Override
    public void close() {
        cplex.end();
    }
}
