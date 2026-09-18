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
    /** 根节点：每个原始物品自成一个有效物品，无冲突约束、空列池。 */
    public Node(Instance inst) {
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
        for (int i = 0; i < inst.n(); i++) {
            items.add(new EffItem(inst.weights[i], new TreeSet<>(Set.of(i))));
        }
    }

    /**
     * 分支节点：在父节点约束基础上施加一条分支决策，并继承父节点的列池。
     *
     * @param inst 问题实例（提供箱子容量）
     * @param a,b  分支的有效物品下标
     * @param anb  true = a、b 必须同箱（together）；false = a、b 必须分箱（separate）
     */
    public Node(Instance inst, Node parent, int a, int b, boolean anb) {
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
        parent.children.add(this);
    }


    // ------------------------------------------------------------------
    //                        together 分支 (a+b)
    // ------------------------------------------------------------------
    private void buildTogether(Instance inst, Node parent, int a, int b) {
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

        // 列池重映射：含 a 或 b 但不同时含二者的列违反 together，剔除
        for (int[] col : parent.columns) {
            boolean ha = false, hb = false;
            for (int e : col) {
                if (e == a) ha = true;
                else if (e == b) hb = true;
            }
            if (ha != hb) continue;
            TreeSet<Integer> set = new TreeSet<>();
            for (int e : col) set.add(remap[e]);
            int[] nc = new int[set.size()];
            int s = 0;
            for (int e : set) nc[s++] = e;
            columns.add(nc);
        }
    }

    // ------------------------------------------------------------------
    //                        separate 分支 (a|b)
    // ------------------------------------------------------------------
    private void buildSeparate(Node parent, int a, int b) {
        items.addAll(parent.items);
        cfc.addAll(parent.cfc);
        cfc.add(new int[]{a, b});

        // 同时含 a、b 的列违反新冲突约束，剔除
        for (int[] col : parent.columns) {
            boolean ha = false, hb = false;
            for (int e : col) {
                if (e == a) ha = true;
                else if (e == b) hb = true;
            }
            if (ha && hb) continue;
            columns.add(col.clone());
        }
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
