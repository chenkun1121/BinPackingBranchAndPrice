import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * 入口：运行演示算例，对比 FFD 贪心上界与 B&P 精确最优解。
 */
public class Main {

    public static void main(String[] args) throws Exception {
        // 演示算例：容量 10，物品 {3×4, 4×2, 8×5}，总重 60。
        // 连续下界 ceil(6)=6，FFD 贪心只得 8 箱，B&P 找到并证明最优解 7 箱：
        //   五个 8 只能独占箱子；剩余 {3×4,4×2} 用 {4,3,3}+{4,3} 两箱装完，
        //   而 FFD 把两个 4 捆在一起导致 3 装不下 —— 体现列生成模式的优越性。

     //   int capacity = 10;
       // double[] weights = {3, 3, 3, 3, 4, 4, 8, 8, 8, 8, 8};

        File file = new File("Scholl/Scholl_3");
        if (!file.exists()) {
            System.out.println("文件夹不存在: " + file.getAbsolutePath());
            return;
        }
        File[] files = file.listFiles();
        if (files == null) {
            System.out.println("无法读取文件夹内容: " + file.getAbsolutePath());
            return;
        }

        long t = System.currentTimeMillis();

        for (File f : files) {
            Instance inst = new Instance();
            inst.initData(f.getAbsolutePath());

            System.out.println("=== Bin Packing — Branch and Price (Java + CPLEX) ===");
            System.out.println("算例: " + f.getName());
            System.out.println(inst);
            System.out.println("连续下界 ceil(总重/容量) = "
                    + (int) Math.ceil(inst.totalWeight() / inst.capacity - 1e-9));

            long t0 = System.currentTimeMillis();
            BranchAndPrice.Result result = new BranchAndPrice(inst).solve();
            long ms = System.currentTimeMillis() - t0;

            System.out.println();
            System.out.println("FFD 贪心上界   : " + result.ffdBound);
            System.out.println("根节点 LP 下界 : " + result.rootLpBound);
            System.out.println("最优箱子数     : " + result.bins);
            System.out.println("B&P 搜索节点数 : " + result.nodes);
            System.out.println("耗时           : " + ms + " ms");
            System.out.println("---------------------------------");
        }
        //System.out.println("最优装箱方案:");
       // printBins(result.solution, inst.weights, inst.capacity);
        System.out.println("总时长" + (System.currentTimeMillis() - t) + " ms");
    }

    private static void printBins(List<Set<Integer>> bins, double[] weights, double capacity) {
        int idx = 1;
        for (Set<Integer> bin : bins) {
            double load = 0;
            for (int it : bin) load += weights[it];
            System.out.printf("  箱 %d: %-16s 装载 %.0f / %.0f%n", idx++, bin, load, capacity);
        }
    }
}
