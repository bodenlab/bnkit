package bn.ctmc.matrix;

import bn.ctmc.SubstModel;
import dat.Enumerable;

/**
 * Created by aesseb on 11-Dec-15.
 */
public class gap extends SubstModel {

    public static Character[] S = {'G', 'C'};

    public static double[] F = {
            //G  C
            0.1,0.9}; //prior probabilities

    public static double [][] Q = {
            // G   C
            {-0.54, 0.54},
            {6.61, -6.61}};

    public gap() {
        super (F, Q, new Enumerable(S), false);
    }

    public String getName() { return "gap"; }

}
