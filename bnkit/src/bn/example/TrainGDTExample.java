package bn.example;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.node.CPT;
import bn.node.GDT;
import dat.EnumVariable;
import dat.Variable;

import java.util.List;

/**
 * Created by mikael on 1/5/18.
 */
public class TrainGDTExample {

    //static String data_file = "/Users/mikael/simhome/melanoma/h3k27me3_test.tsv";
    static String data_file = "/Users/mikael/simhome/melanoma/h3k27me3_test_signals_only.tsv";
    static String bn_file = "/Users/mikael/simhome/melanoma/bn_h3k27me3_test.xml";

    public static void main(String[] args) {

        // Define variables
        EnumVariable R = Predef.Boolean("Repressor");
        EnumVariable L = Predef.Boolean("Latent factor");
        Variable K = Predef.Real("H3K27me3");
        Variable Y = Predef.Real("Log RNA expression");

        // Define nodes
        CPT r = new CPT(R);
        CPT l = new CPT(L);
        GDT k = new GDT(K, R);    // P(K | R)
        GDT y = new GDT(Y, R, L); // P(Y | R, L)

        // let variances differ
        //k.setTieVariances(0);
        //y.setTieVariances(0);

        // Define BN
        BNet bn = new BNet();
        bn.add(r, l, k, y);
        List<BNode> nodes = bn.getOrdered(); // topologically ordered list of all nodes (i.e. parents always precede children)
        for (BNode node : nodes)
            System.out.println(node);

        // Load data
        Object[][] data = DataBuf.load(data_file, nodes);

        // Train
        EM em = new EM(bn);
        em.train(data, nodes);
        BNBuf.save(bn, bn_file);
    }

}
