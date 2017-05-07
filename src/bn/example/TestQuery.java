package bn.example;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import dat.EnumVariable;
import dat.Variable;
import dat.Variable.Assignment;
import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.prob.EnumDistrib;

public class TestQuery {
	/**
	 * @param args
	 */

	public static void main(String[] args) {
        String setting = "";
        int max = 0;
        String bn_file = "";
        String data_file = "";
        String query = "";

		if (args.length == 4) {
            bn_file = args[0];
            data_file = args[1];
            query = args[2];
            setting = args[3];
            if (setting.equals("Sample")) {
                System.out.println("Usage: cannot perform sampling without a max value");
                System.exit(1);
            } else if (setting.equals("Infer")) {
                System.out.println("Usage: cannot perform sampling without a max value");
                System.exit(1);
            }
		} else if (args.length == 5) {
            bn_file = args[0];
            data_file = args[1];
            query = args[2];
            setting = args[3];

            if (setting.equals("Sample") || setting.equals("Infer")) {
                max = Integer.parseInt(args[4]);
            } else {
                System.out.println("Usage: not using Sample or Infer so max value will be ignored");
            }
        } else {
            System.out.println("Usage: LoadNTrain <bn-file> <data-file> <query> <setting>");
            System.out.println("Usage: LoadNTrain <bn-file> <data-file> <query> Sample <max>");
            System.exit(1);
        }

		//"micro_data/1413340451503m41-Full-nu-dsd-ut.new" \
		//"micro_data/GenomicLocation.txt" \
		//"VarianceC" "Infer"

        TestQueryCases tq = new TestQueryCases(bn_file, data_file, query, setting, max);

	}
}
