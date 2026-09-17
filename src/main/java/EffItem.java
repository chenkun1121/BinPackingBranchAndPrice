import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

/**
 * "有效物品"（effective item）。
 *
 * Ryan-Foster 的"同箱"（together）分支会把若干原始物品强制绑定进同一个箱子，
 * 实现方式是把它们合并成一个不可分割的有效物品：整组要么进箱、要么不进。
 * 根节点上每个有效物品恰对应一个原始物品。
 */
public class EffItem {
    public  double weight;          // 组内原始物品重量之和
    public  Set<Integer> origItems; // 组内包含的原始物品编号（用于还原最终方案）
    // 注意：origItems 是一个不可变的 TreeSet，保证元素的有序性和唯一性


    // 用于在分支过程中快速查询与当前物品不同箱的有效物品集合
    public EffItem(double weight, Set<Integer> origItems) {
        this.weight = weight;
        this.origItems = Collections.unmodifiableSet(new TreeSet<>(origItems));
    }

    @Override
    public String toString() { return origItems + "(w=" + weight + ")"; }
}
