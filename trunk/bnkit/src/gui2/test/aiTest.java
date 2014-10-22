package gui2.test;


import java.util.Collection;
import java.util.List;
import java.util.Map;

import bn.BNet;
import bn.BNode;
import bn.node.CPT;
import bn.Distrib;
import bn.EnumDistrib;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Enumerable;
import bn.JPT;
import bn.Predef;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.alg.QueryResult;
import dat.Variable;
import bn.file.BNBuf;
import bn.file.DataBuf;

public class aiTest {

    public static void main(String[] args) {
//        if (args.length < 1) {
//            System.out.println("Usage: LoadNTrain <bn-file>");
//            System.exit(1);
//        }
//        String bn_file = args[0];
        String bn_file = "cgSimple3.new";
        BNet bn = BNBuf.load(bn_file);

        for (BNode node : bn.getNodes()) {
            if (node.getName().equals("DNase(Open)")) {
                node.setInstance(10.54);
            }
//			if(node.getName().equals("DNase(UWash)")){
//				node.setInstance(10.53);
//			}
//			if(node.getName().equals("Chromatin")){
//				node.setInstance(true);
//			}
            if (node.getName().equals("RepeatSeq")) {
                node.setInstance("two");
            }
//			if(node.getName().equals("Proxy")){
//				node.setInstance(true);
//			}
            if (node.getName().equals("Variance")) {
                node.setInstance(false);
            }
            if (node.getName().equals("Unstable")) {
                node.setInstance(false);
            }

        }

        System.out.println("Variable Elimination------------");
        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn);

        Query q = ve.makeQuery(bn.getNode("DNase(UWash)").getVariable());

        CGTable res = (CGTable) ve.infer(q);
        res.display();
        res.displaySampled();

    }
}
