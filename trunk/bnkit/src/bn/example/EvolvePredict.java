package bn.example;

import bn.BNet;
import bn.BNode;
import bn.CPT;
import bn.EnumDistrib;
import bn.EnumVariable;
import bn.Enumerable;
import bn.GDT;
import bn.GaussianDistrib;
import bn.JPT;
import bn.Predef;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.Variable;
import bn.alg.CGVarElim;
import bn.alg.QueryResult;
import bn.file.BNBuf;
import bn.file.DataBuf;

import java.util.List;

/**
 * Created by julianzaugg on 7/04/14.
 */
public class EvolvePredict {

    public static void main(String[] args){
        Variable ESELECT = Predef.Real("E_val");
        EnumVariable SCORE = Predef.Nominal(new String[] {"5-25", "26-50", "51-75", "76-158"}, "Score_group");
//        EnumVariable SCORE = Predef.Nominal(new String[] {"5-51", "52-158"}, "Score_group");

        EnumVariable CUR215 = Predef.AminoAcid("Current_215");
        EnumVariable CUR219 = Predef.AminoAcid("Current_219");
        EnumVariable CUR244 = Predef.AminoAcid("Current_244");
        EnumVariable CUR249 = Predef.AminoAcid("Current_249");
        EnumVariable CUR317 = Predef.AminoAcid("Current_317");
        EnumVariable CUR318 = Predef.AminoAcid("Current_318");
        EnumVariable CUR349 = Predef.AminoAcid("Current_349");
        EnumVariable CUR350 = Predef.AminoAcid("Current_350");

        EnumVariable NEXT215 = Predef.AminoAcid("Next_215");
        EnumVariable NEXT219 = Predef.AminoAcid("Next_219");
        EnumVariable NEXT244 = Predef.AminoAcid("Next_244");
        EnumVariable NEXT249 = Predef.AminoAcid("Next_249");
        EnumVariable NEXT317 = Predef.AminoAcid("Next_317");
        EnumVariable NEXT318 = Predef.AminoAcid("Next_318");
        EnumVariable NEXT349 = Predef.AminoAcid("Next_349");
        EnumVariable NEXT350 = Predef.AminoAcid("Next_350");

        GDT eselect = new GDT(ESELECT, SCORE);
        //Classes are based on approximate quantile groups of the distribution
        //data itself was not normalized and is very slightly skewed.
//        eselect.put(new GaussianDistrib(28.19, 12), "5-51");
//        eselect.put(new GaussianDistrib(80.9, 23), "52-158");
        eselect.put(new GaussianDistrib(18.7, 6.15), "5-25");
        eselect.put(new GaussianDistrib(37.55, 6.9), "26-50");
        eselect.put(new GaussianDistrib(62.24, 7.74), "51-75");
        eselect.put(new GaussianDistrib(100.28, 18.43), "76-158");

        CPT score = new CPT(SCORE);
//        score.put(new EnumDistrib(new Enumerable(new String[] {"5-51", "52-158"}),
//                0.5, 0.5));
        score.put(new EnumDistrib(new Enumerable(new String[] {"5-25", "26-50", "51-75", "76-158"}),
                0.25, 0.25, 0.25, 0.25));


        CPT cur215 = new CPT(CUR215, SCORE);
//        cur215.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur219 = new CPT(CUR219, SCORE);
//        cur219.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur244 = new CPT(CUR244, SCORE);
//        cur244.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur249 = new CPT(CUR249, SCORE);
//        cur249.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur317 = new CPT(CUR317, SCORE);
//        cur317.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur318 = new CPT(CUR318, SCORE);
//        cur318.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur349 = new CPT(CUR349, SCORE);
//        cur349.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT cur350 = new CPT(CUR350, SCORE);
//        cur350.put(new EnumDistrib(Enumerable.aacid), "5-25");

        CPT next215 = new CPT(NEXT215, CUR215);

        CPT next219 = new CPT(NEXT219, CUR219);

        CPT next244 = new CPT(NEXT244, CUR244);

        CPT next249 = new CPT(NEXT249, CUR249);

        CPT next317 = new CPT(NEXT317, CUR317);

        CPT next318 = new CPT(NEXT318, CUR318);

        CPT next349 = new CPT(NEXT349, CUR349);

        CPT next350 = new CPT(NEXT350, CUR350);

        BNet bn = new BNet();
        bn.add(eselect);
        bn.add(score);
        bn.add(cur215,cur219,cur244, cur249, cur317, cur318, cur349, cur350);
        bn.add(next215,next219,next244, next249, next317, next318, next349, next350);

        List<BNode> nodes = bn.getOrdered();

        String data_file = "/Users/mikael/Desktop/feature_vectors.txt";
        Object[][] data = DataBuf.load(data_file, nodes);

        LearningAlg em = new EM(bn);
        EM.EM_MAX_ROUNDS = 10;
        em.train(data, nodes);
        System.out.println("DONE TRAINING");
//        for (BNode node: nodes){
//            node.print();
//        }

        cur215.setInstance('F');
        cur219.setInstance('V');
        cur244.setInstance('F');
        cur249.setInstance('F');
        cur317.setInstance('F');
        cur318.setInstance('I');
        cur349.setInstance('L');
        cur350.setInstance('C');
//        eselect.setInstance(55.0);
//        score.setInstance("5-25");

        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn);

        Query q = ve.makeQuery(new EnumVariable[] {NEXT215});
        QueryResult qr = ve.infer(q);
        JPT jpt = qr.getJPT();
        jpt.display();
    }

}
