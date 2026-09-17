import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
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
 * （无冲突约束时可换成动态规划背包，复杂度 O(nC)；这里统一用 CPLEX 精确求解。）
 */
public final class PricingProblem {

    private static final double EPS = 1e-6;

    private PricingProblem() {}

    /**
     * @param items     当前节点的有效物品
     * @param forbidden 禁止同箱的有效物品下标对
     * @param pi        主问题覆盖约束的对偶价格
     * @param capacity  箱子容量
     * @return 价值超过 1 的新模式（装箱方案）；无改进列时返回 null
     */
    public static int[] solve(List<EffItem> items, List<int[]> forbidden,
                              double[] pi, double capacity) throws IloException {
        int n = items.size();
        IloCplex cplex = new IloCplex();
        try {
            cplex.setOut(null);
            IloNumVar[] x = cplex.boolVarArray(n);
            // IloNumExpr obj = cplex.numExpr();
            //重新使用 scalProd 计算目标函数和约束，建立定价问题模型
            cplex.addMaximize(cplex.scalProd(x, pi));
            double[] weights = items.stream().mapToDouble(e -> e.weight).toArray();
            cplex.addLe(cplex.scalProd(x, weights), capacity);

//            IloLinearNumExpr obj = cplex.linearNumExpr();
//            IloLinearNumExpr load = cplex.linearNumExpr();
//            for (int e = 0; e < n; e++) {
//                obj.addTerm(pi[e], x[e]);
//                load.addTerm(items.get(e).weight, x[e]);
//            }
//            cplex.add(cplex.maximize(obj));
//            cplex.addLe(load, capacity);
            for (int[] p : forbidden) {                   // 分箱分支引入的冲突约束
                cplex.addLe(cplex.sum(x[p[0]], x[p[1]]), 1.0);
            }

            if (!cplex.solve() || cplex.getStatus() != IloCplex.Status.Optimal) {
                return null;
            }
            if (cplex.getObjValue() <= 1.0 + EPS) {       // 无负 reduced cost 列
                return null;
            }
            List<Integer> pattern = new ArrayList<>();
            for (int e = 0; e < n; e++) {
                // 检查变量 x[e] 的值是否大于 0.5，如果是，则将该有效物品的下标 e 添加到模式中
                if (cplex.getValue(x[e]) > 0.5) pattern.add(e);
            }
            int[] result = new int[pattern.size()];
            for (int i = 0; i < result.length; i++) result[i] = pattern.get(i);
            return result;
        } finally {
            cplex.end();
        }
    }
}
