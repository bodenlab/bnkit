package bn.example;

import bn.*;
import bn.alg.*;
import bn.file.BNBuf;
import bn.file.DataBuf;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
//        eselect.put(new GaussianDistrib(18.7, 6.15), "5-25");
//        eselect.put(new GaussianDistrib(37.55, 6.9), "26-50");
//        eselect.put(new GaussianDistrib(62.24, 7.74), "51-75");
//        eselect.put(new GaussianDistrib(100.28, 18.43), "76-158");

        CPT score = new CPT(SCORE);
        CPT cur215 = new CPT(CUR215, SCORE);
        CPT cur219 = new CPT(CUR219, SCORE);
        CPT cur244 = new CPT(CUR244, SCORE);
        CPT cur249 = new CPT(CUR249, SCORE);
        CPT cur317 = new CPT(CUR317, SCORE);
        CPT cur318 = new CPT(CUR318, SCORE);
        CPT cur349 = new CPT(CUR349, SCORE);
        CPT cur350 = new CPT(CUR350, SCORE);
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
        Object[][] data = DataBuf.load(data_file, nodes); // load full data file
        LearningAlg em = new EM(bn);

        em.train(data, nodes);


//        cur215.setInstance('F');
//        cur219.setInstance('V');
//        cur244.setInstance('F');
//        cur249.setInstance('F');
//        cur317.setInstance('F');
//        cur318.setInstance('I');
//        cur349.setInstance('L');
//        cur350.setInstance('C');
//        eselect.setInstance(55.0);
//        score.setInstance("5-25");

        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn.getRelevant(new Variable[] {NEXT349}));
        Query q = ve.makeQuery(new EnumVariable[] {NEXT349});
        JPT jpt = ve.infer(q).getJPT();
        jpt.display();

        System.exit(0);

        /**
         * Perform a 'leave one out' cross validation
         */
        Object[][] training_set; // this will hold training set
        for (int i = 0; i < data.length; i++){
            Object[] test_sample = data[i]; // This is our test sample
            training_set = new Object[data.length - 1][];
            int count = 0;
            //Go through and add vectors to training set
            for (int j = 0; j < data.length; j++){
                if (j != i){
                    training_set[count] = data[j];
                    count++;
                }
            }
            // Train network on the training data
            em.train(training_set, nodes);
            //Evaluate the test sample
            System.out.println(Arrays.toString(test_sample));
//            System.exit(0);
        }

//        System.out.println("DONE TRAINING");



        //FBE
//        cur215.setInstance('F');
//        cur219.setInstance('V');
//        cur244.setInstance('F');
//        cur249.setInstance('F');
//        cur317.setInstance('F');
//        cur318.setInstance('I');
//        cur349.setInstance('L');
//        cur350.setInstance('C');
//        eselect.setInstance(55.0);
//        score.setInstance("5-25");


//        CGVarElim ve = new CGVarElim();
//        ve.instantiate(bn.getRelevant(new Variable[] {NEXT349,NEXT350}));
//        Query q = ve.makeQuery(new EnumVariable[] {NEXT349,NEXT350});
//        JPT jpt = ve.infer(q).getJPT();
//        jpt.display();

    }

//    public static void main2(String[] args){
//        String data_file = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary Pathway/Data/ANEH/Data/feature_vectors.txt";
//        Object[][] data = DataBuf.load(data_file, nodes);
////        BNet bn = BNBuf.load(trained_bn);
//
////        for (int i = 0; i < )
//    }

    public static void main1(String[] args){
        String trained_bn = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary Pathway/Code/bnet_zero/data/JZtest.xml";
        BNet bn = BNBuf.load(trained_bn);
        List<BNode> nodes = bn.getOrdered();

        BNode cur215 = bn.getNode("Current_215");
        BNode cur219 = bn.getNode("Current_219");
        BNode cur244 = bn.getNode("Current_244");
        BNode cur249 = bn.getNode("Current_249");
        BNode cur317 = bn.getNode("Current_317");
        BNode cur318 = bn.getNode("Current_318");
        BNode cur349 = bn.getNode("Current_349");
        BNode cur350 = bn.getNode("Current_350");
        BNode score = bn.getNode("Score_group");
        BNode e_val = bn.getNode("E_val");
        cur215.setInstance('F');
        cur219.setInstance('D');
        cur244.setInstance('F');
        cur249.setInstance('Y');
        cur317.setInstance('V');
        cur318.setInstance('F');
        cur349.setInstance('V');
        cur350.setInstance('H');

        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn.getRelevant(new Variable[] {score.getVariable()}));
        Query q = ve.makeQuery(new Variable[] {score.getVariable()});
        QueryResult qr = ve.infer(q);
//        Distrib thing = ((CGVarElim.CGResult)qr).getNonEnumDistrib().get(e_val);
        JPT jpt = qr.getJPT();
        jpt.display();
//        thing.get(e_val);
//        double sum = 0;
//        for (int i = 0; i < 100; i++) {
//            sum += (Double)thing.sample();
//        }
//        System.out.println(sum);
    }
}
