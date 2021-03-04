package bn.ctmc.matrix;

import bn.ctmc.SubstModel;
import dat.Enumerable;

/**
 * Based on Cohen and Pupko (2010)
 * Models for evolutionary loss and gain
 * M1: g = 0.54 l = 6.61 gamma = 0.73
 *
 */
public class GLOOME1 extends SubstModel {

    // true means the "jump" interval is there so "loss"
    // false means absence of "jump" so "gain"
    public static Enumerable bool = new Enumerable(new Boolean[]{true, false});

    public static double[] F = {
            //G  L
            0.5,0.5}; //prior probabilities

//    public static double [][] Q = { // note order is not the same as in Cohen and Pupko, since the meaning of the rows are different
//            {-6.61, 6.61},
//            { 0.54,-0.54},
//    };
    public static double [][] Q = { // note order is not the same as in Cohen and Pupko, since the meaning of the rows are different
            {-0.54, 0.54},
            { 6.61,-6.61},
    };

    public GLOOME1() {
        super (F, Q, bool, false);
    }

    public String getName() { return "GLOOME1"; }

}
