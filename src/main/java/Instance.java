import ilog.concert.IloNumVar;

import java.io.BufferedReader;
import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

/**
 * 装箱问题（Bin Packing）实例：
 * n 个不可分割物品，重量 w_i；箱子数量不限，每个箱子容量 C；
 * 目标：最小化使用的箱子数。
 */
public class Instance {
    public double[] weights;
    public int capacity;

    public Instance() {

    }

    public void initData(String path) throws Exception {
        try {
            /*
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line = br.readLine();
            line = br.readLine();
            String[] tokens = line.split("\\s+");
            capacity = Integer.parseInt(tokens[0]);
            int n = Integer.parseInt(tokens[1]);
            weights = new double[n];
            for(int i =0; i <3; i++){
                line = br.readLine();
            }
            for(int i = 0; i < n; i++){
                line = br.readLine();
                weights[i] = Double.parseDouble(line);
            }
             */
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line = br.readLine();
            int n = Integer.parseInt(line);
            capacity = Integer.parseInt(br.readLine());
            weights = new double[n];
            for (int i = 0; i < n; i++) {
                weights[i] = Double.parseDouble(br.readLine());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int n() {
        return weights.length;
    }

    public double totalWeight() {
        return Arrays.stream(weights).sum();
    }

    @Override
    public String toString() {
        return "Instance{n=" + n() + ", capacity=" + capacity
                + ", weights=" + Arrays.toString(weights) + "}";
    }
}
