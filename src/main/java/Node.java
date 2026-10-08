import ilog.concert.IloException;
import ilog.cplex.IloCplex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * The subproblem at a node in the branch-and-bound tree.
 *
 * 每个节点 = 一组局部约束 + 一份继承列池 + 求解结果：
 *   items   : 有效物品（agg）。together 分支把被绑定的物品合并成一个不可分割的 EffItem，
 *             于是"必须同箱"被编码进问题数据本身（物品数减一）。
 *   cfc     : 冲突对（conflicts）。separate 分支把 (a,b) 加进来，定价子问题据此加
 *             x_a + x_b <= 1。等价于参考文献中的 cfc 矩阵。
 *   columns : 列池。父节点传下来的列（warm start）经相容性过滤 + 下标重映射后继承，
 *             列生成结束后回写为本节点的全部列，供子节点继续复用。
 *
 * 节点 LP 松弛（集合覆盖）：
 *   min  sum_s y_s
 *   s.t. sum_{s ni e} y_s >= 1     for all effective item e     (dual pi_e >= 0)
 *        y_s >= 0
 *   其中 s 为可行模式：sum_{e in s} w_e <= C 且 s 不含任何 cfc 对。
 *
 * 分支沿用 Ryan-Foster：对物品对 (a,b) 分两支——
 *   anb = true  : a、b 同箱（"a+b"），合并为有效物品；
 *   anb = false : a、b 分箱（"a|b"），加入冲突对。
 */
public class Node {

    public static int cnt;                  // 节点计数器
    public static double EPS = 1e-6;

    /** 节点状态。 */
    public enum Status {
        NEW,            // 刚生成，尚未列生成
        SOLVED,         // 列生成收敛，已得 LP 下界
        INFEASIBLE,     // RMP 不可行
        PRUNED,         // 被上界剪枝
        INTEGER,        // LP 解已整，得到可行解
        BRANCHED        // LP 解分数，已分支
    }

    public Node parent;
    public ArrayList<Node> children;
    public int no;                          // index of the node
    public int depth;
    public double lb;                       // lower bound of the subtree at the node
    public int n_lp;                        // 该节点列生成迭代次数

    public ArrayList<EffItem> items;        // aggregate items，表示必须同箱的原始物品集合
    public ArrayList<int[]> cfc;            // items in conflicts，表示必须分箱的有效物品对
    public ArrayList<int[]> columns;        // 列池（warm start + 本节点生成的列）
    public double[] y;                      // 各列的 LP 取值
    public int[] branchPair;                // 本节点选出的 Ryan-Foster 分支对

    public ArrayList<String> last_branch;   // 分支历史: "a+b" / "a|b"
    public Status status;
    public boolean feasible;                // together 分支矛盾或超容量时为 false


    private HashMap<Integer, BitSet> incompatible;//快速查询与i不同箱的有效物品集合
    /** 原始物品总数（全树固定，不随 together 合并而变）。 */
    private final int numOrig;
    /** 原始物品 -> 有效物品下标的映射表，items 定稿后构建。 */
    private int[] origToEff;

    /** 根节点：每个原始物品自成一个有效物品，无冲突约束、空列池。 */
    public Node(Instance inst) throws IloException {
        this.no = Node.cnt++;
        this.parent = null;
        this.depth = 0;
        this.children = new ArrayList<>();
        this.last_branch = new ArrayList<>();
        this.items = new ArrayList<>();
        this.cfc = new ArrayList<>();
        this.columns = new ArrayList<>();
        this.status = Status.NEW;
        this.feasible = true;
        this.numOrig = inst.n();
        for (int i = 0; i < inst.n(); i++) {
            items.add(new EffItem(inst.weights[i], new TreeSet<>(Set.of(i))));
        }
        rebuildOrigToEff();
    }

    /**
     * 分支节点：在父节点约束基础上施加一条分支决策，并继承父节点的列池。
     *
     * @param inst 问题实例（提供箱子容量）
     * @param a,b  分支的有效物品下标
     * @param anb  true = a、b 必须同箱（together）；false = a、b 必须分箱（separate）
     */
    public Node(Instance inst, Node parent, int a, int b, boolean anb) throws IloException {
        this.parent = parent;
        this.no = Node.cnt++;
        this.depth = parent.depth + 1;
        this.lb = parent.lb;
        this.children = new ArrayList<>();
        this.last_branch = new ArrayList<>(parent.last_branch);
        this.items = new ArrayList<>();
        this.cfc = new ArrayList<>();
        this.columns = new ArrayList<>();
        this.status = Status.NEW;
        this.feasible = true;
        this.numOrig = parent.numOrig;
        // 列池不再逐节点复制：由 MasterProblem 全局持有，模式用原始物品下标表示

        if (a > b) {                        // ensure a <= b
            int tmp = b;
            b = a;
            a = tmp;
        }

        if (anb) {
            last_branch.add(0, a + "+" + b);
            buildTogether(inst, parent, a, b);
        } else {
            last_branch.add(0, a + "|" + b);
            buildSeparate(parent, a, b);
        }
        rebuildOrigToEff();                  // items 定稿后再建映射
        parent.children.add(this);
    }


    // ------------------------------------------------------------------
    //                        together 分支 (a+b)
    // ------------------------------------------------------------------
    private void buildTogether(Instance inst, Node parent, int a, int b) throws IloException {
        // 与已有的 separate 决策矛盾
        for (int[] p : parent.cfc) {
            if (p[0] == a && p[1] == b) {
                feasible = false;
                return;
            }
        }
        double w = parent.items.get(a).weight + parent.items.get(b).weight;
        if (w > inst.capacity + EPS) {      // 超容量 → 该子树不可行
            feasible = false;
            return;
        }

        // 合并a、b有效物品的原始物品集合
        TreeSet<Integer> orig = new TreeSet<>(parent.items.get(a).origItems);
        orig.addAll(parent.items.get(b).origItems);

        // 新物品列表：a、b 删除，合并物品追加到末尾
        for (int k = 0; k < parent.items.size(); k++) {
            if (k == a || k == b) continue;
            items.add(parent.items.get(k));
        }
        int m = items.size();
        items.add(new EffItem(w, orig));

        // 下标重映射：a、b 都指向合并物品 m，其他物品按原顺序重新编号
        int[] remap = new int[parent.items.size()];
        for (int k = 0, s = 0; k < remap.length; k++) {
            // 重新映射下标
            if (k == a || k == b) remap[k] = m;
            else remap[k] = s++;
        }

        // 冲突对重映射；两端塌缩成同一下标则该约束消失，其他约束按新下标加入 cfc
        for (int[] p : parent.cfc) {
            int x = remap[p[0]];
            int z = remap[p[1]];
            if (x != z) cfc.add(new int[]{Math.min(x, z), Math.max(x, z)});
        }
        // 列池不再随节点复制 / 重映射：列池由 MasterProblem 全局持有，
        // 模式一律以"原始物品下标"表示，节点切换靠变量上界屏蔽（见 MasterProblem.setNode）。
    }

    // ------------------------------------------------------------------
    //                        separate 分支 (a|b)
    // ------------------------------------------------------------------
    private void buildSeparate(Node parent, int a, int b) {
        items.addAll(parent.items);
        cfc.addAll(parent.cfc);
        cfc.add(new int[]{a, b});
        // 同上：列池不再逐节点复制，交给 MasterProblem.setNode 屏蔽。
    }

    // ------------------------------------------------------------------
    //                              接口
    // ------------------------------------------------------------------

    public boolean isRootNode() {
        return parent == null;
    }

    /** together 分支是否可行（未超容量、未与既有冲突约束矛盾）。 */
    public boolean isFeasible() {
        return feasible;
    }

    // ------------------------------------------------------------------
    //          全局列池下的节点视图（列一律用"原始物品下标"表示）
    // ------------------------------------------------------------------

    /**
     * 判定一个模式（原始物品下标）在本节点是否合法。这是 {@code MasterProblem.setNode}
     * 用来屏蔽 / 放开列的谓词，对应 Ryan-Foster 的两类分支约束：
     *
     *   1) together：每个有效物品必须"整组进箱或整组不进"，不允许只取其中一部分原始物品；
     *   2) separate：模式不得同时包含某个冲突对两端的有效物品。
     */
    public boolean isColumnValid(int[] origPattern) {
        BitSet b = new BitSet();
        for (int o : origPattern) b.set(o);

        // 1) together：整块判定（0 < 命中数 < 组大小 即为非法）
        for (EffItem e : items) {
            int cnt = 0;
            for (int o : e.origItems) {
                if (b.get(o)) cnt++;
            }
            if (cnt > 0 && cnt < e.origItems.size()) return false;
        }

        // 2) separate：冲突对两端不得同时出现在模式中
        for (int[] p : cfc) {
            if (hits(b, items.get(p[0])) && hits(b, items.get(p[1]))) return false;
        }
        return true;
    }

    /** 模式是否命中该有效物品（含其任一原始物品即算命中）。 */
    private static boolean hits(BitSet b, EffItem e) {
        for (int o : e.origItems) {
            if (b.get(o)) return true;
        }
        return false;
    }

    /** 把定价子问题返回的"有效物品下标"模式展开成"原始物品下标"。 */
    public int[] toOrigPattern(int[] localPattern) {
        TreeSet<Integer> orig = new TreeSet<>();
        for (int e : localPattern) {
            orig.addAll(items.get(e).origItems);
        }
        int[] r = new int[orig.size()];
        int k = 0;
        for (int o : orig) r[k++] = o;
        return r;
    }

    /** 第 e 个有效物品对应的原始物品下标数组（用于构造"单组一箱"的兜底列）。 */
    public int[] origPatternOf(int effIndex) {
        Set<Integer> s = items.get(effIndex).origItems;
        int[] r = new int[s.size()];
        int k = 0;
        for (int o : s) r[k++] = o;      // origItems 是 TreeSet，天然升序
        return r;
    }

    /** 重建"原始物品 -> 有效物品下标"映射。items 每次定稿后调用。 */
    private void rebuildOrigToEff() {
        origToEff = new int[numOrig];
        Arrays.fill(origToEff, -1);
        for (int e = 0; e < items.size(); e++) {
            for (int o : items.get(e).origItems) origToEff[o] = e;
        }
    }

    /**
     * 把"原始物品下标"模式映射回本节点的"有效物品下标"模式（去重、升序）。
     * 分支决策必须在有效物品层面进行——together 合并后的组才是一个不可分割的单位，
     * 对组内某个原始物品单独分支会破坏分支的正确性。
     */
    public int[] toEffPattern(int[] origPattern) {
        TreeSet<Integer> eff = new TreeSet<>();
        for (int o : origPattern) {
            int e = origToEff[o];
            if (e >= 0) eff.add(e);
        }
        int[] r = new int[eff.size()];
        int k = 0;
        for (int e : eff) r[k++] = e;
        return r;
    }

    /**
     * 把 RMP 的"原始物品"对偶价格 π 聚合到本节点的"有效物品"层面。
     * together 分支把若干原始物品绑定成一个有效物品，它们必须同进同出，
     * 因此在定价子问题看来该有效物品的价值等于其成员对偶之和：
     *      π_eff(e) = Σ_{o ∈ e.origItems} π_o
     */
    public double[] effDuals(double[] origPi) {
        double[] piEff = new double[items.size()];
        for (int e = 0; e < items.size(); e++) {
            double s = 0;
            for (int o : items.get(e).origItems) s += origPi[o];
            piEff[e] = s;
        }
        return piEff;
    }

    public ArrayList<Node> getChildren() {
        return children;
    }

    /** 与 i 禁止同箱的有效物品集合（BitSet 形式，供定价子问题使用）。 */
    public BitSet getIncompatible(int i) {
        if (incompatible == null) incompatible = new HashMap<>();
        if (!incompatible.containsKey(i)) {
            BitSet bitSet = new BitSet(items.size());
            for (int[] p : cfc) {
                if (p[0] == i) bitSet.set(p[1]);
                else if (p[1] == i) bitSet.set(p[0]);
            }
            incompatible.put(i, bitSet);
        }
        return incompatible.get(i);
    }

    /** 从根到本节点的分支决策串，如 "3+7 -> 2|5"。 */
    public String getBranchPath() {
        List<String> path = new ArrayList<>(last_branch);
        java.util.Collections.reverse(path);
        return String.join(" -> ", path);
    }

    @Override
    public String toString() {
        return "Node{no=" + no + ", depth=" + depth + ", lb="
                + (Double.isNaN(lb) ? "n/a" : String.format("%.4f", lb))
                + ", items=" + items.size() + ", cfc=" + cfc.size()
                + ", cols=" + columns.size() + ", branch=" + getBranchPath()
                + ", status=" + status + "}";
    }

    /** 列的规范化键，用于列池去重，避免 RMP 中出现重复列。 */
    public static String columnKey(int[] col) {
        int[] c = col.clone();
        Arrays.sort(c);
        return Arrays.toString(c);
    }
}
