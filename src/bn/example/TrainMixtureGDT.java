package bn.example;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.EM;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.node.CPT;
import bn.node.GDT;
import dat.EnumVariable;
import dat.Variable;
import dat.file.TSVFile;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Created by mikael on 1/5/18.
 */
public class TrainMixtureGDT {

    public static void main(String[] args) {

        if (args.length < 2) {
            System.out.println("Usage: TrainMixtureGDT <TSV-data> <BN-file>");
            System.exit(1);
        }
        String data_file = args[0];
        String bn_file = args[1];

        // Define variables
        EnumVariable C = Predef.Boolean("Mixer");
        Variable     Y = Predef.Real("Signal");

        // Define nodes
        CPT c = new CPT(C);
        GDT y = new GDT(Y, C);

        // variance tied?
        y.setTieVariances(2); //  (0 is untied, 1 is use the max, 2 is pooled)

        // Define BN
        BNet bn = new BNet();
        bn.add(c, y);
        List<BNode> nodes = bn.getOrdered(); // topologically ordered list of all nodes (i.e. parents always precede children)
        for (BNode node : nodes)
            System.out.println(node);

        // Load data

            Object[][] data = DataBuf.load(data_file, Arrays.asList(new BNode[] {y}), false);

            // Train
            EM em = new EM(bn);
            em.setEMOption(2);
            em.train(data, Arrays.asList(new BNode[] {y}));
            BNBuf.save(bn, bn_file);
    }

}
