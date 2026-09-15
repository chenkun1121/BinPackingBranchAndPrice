import ilog.concert.IloColumn;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.List;

/**
 * 受限主问题（Restricted Master Problem，RMP）—— 集合覆盖（Set Covering）松弛：
 *
 *   min  Σ_{s∈S} y_s
 *   s.t. Σ_{s∈S: e∈s} y_s ≥ 1     ∀ 有效物品 e      （对偶变量 π_e ≥ 0）
 *        y_s ≥ 0
 *
 * 每一列 s 是一个可行的箱子装载模式（哪些有效物品放进同一个箱子）。
 * 初始只有"单物品"列；定价子问题不断向 S 追加 reduced cost < 0 的新列，
 * 直到不存在这样的列，RMP 的 LP 最优值即原问题在该节点的下界。
 */
public class MasterProblem implements AutoCloseable {

    private final IloCplex cplex;
    private final IloObjective objective;
    private final IloRange[] cover;                    // 每个有效物品一条覆盖约束
    private final List<IloNumVar> vars = new ArrayList<>();
    private final List<int[]> columns = new ArrayList<>(); // 每列覆盖的有效物品下标

    public MasterProblem(int numItems) throws IloException {
        cplex = new IloCplex();
        cplex.setOut(null);                            // 关闭求解日志
        objective = cplex.addMinimize();
        cover = new IloRange[numItems];
        for (int e = 0; e < numItems; e++) {
            cover[e] = cplex.addRange(1.0, Double.MAX_VALUE); // Σ y_s ≥ 1
        }
    }

    /** 追加一列（一个箱子模式），pattern 为该模式覆盖的有效物品下标数组。 */
    public void addColumn(int[] pattern) throws IloException {
        IloColumn col = cplex.column(objective, 1.0);        // 目标系数 1
        for (int e : pattern) {
            col = col.and(cplex.column(cover[e], 1.0));      // 覆盖约束系数 1
        }
        vars.add(cplex.numVar(col, 0.0, Double.MAX_VALUE, "y" + columns.size()));
        columns.add(pattern);
    }

    /** 求解 RMP，返回覆盖约束的对偶价格 π_e（不可行返回 null）。 */
    public double[] solve() throws IloException {
        if (!cplex.solve() || cplex.getStatus() != IloCplex.Status.Optimal) {
            return null;
        }
        double[] pi = new double[cover.length];
        for (int e = 0; e < cover.length; e++) {
            pi[e] = cplex.getDual(cover[e]);
        }
        return pi;
    }

    public double objectiveValue() throws IloException { return cplex.getObjValue(); }

    /** 各列当前的 LP 取值（用于分支决策与整数性判定）。 */
    public double[] columnValues() throws IloException {
        double[] v = new double[vars.size()];
        for (int i = 0; i < vars.size(); i++) {
            v[i] = cplex.getValue(vars.get(i));
        }
        return v;
    }

    /** 获取所有列。 */
    public List<int[]> getColumns() { return columns; }

    @Override
    public void close() {
        cplex.end();
    }
}
