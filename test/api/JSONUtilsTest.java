package api;

import json.JSONObject;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertTrue;

class JSONUtilsTest {

    static Object[][] data6 = new Object[][] { // 3 samples, no items, 6 features
            {1,2,3,4,null,6},
            {6,3,1,null,-2,-4},
            {null,null,0,0,0,0}    };
    static Object[][][] data2 = new Object[][][] { // 5 samples, 2 items (e.g. proteins), 3 feature
            { {1,"a",null}, {2,"b",3.14}},
            { {10,"A",-0.2}, {20,"B",3.1}},
            { {-1,"a",1.11}, {-2,"b",3.4}},
            { {2,"a",null}, {3,"b",3.14}},
            { {11,"A",-0.2}, {22,"B",3.1}}
    };

    static JSONUtils.DataSet ds6 = new JSONUtils.DataSet(new String[] {"N1","N2","N3","N4","N5","N6"}); // three samples, each data point has six features
    static JSONUtils.DataSet ds2 = new JSONUtils.DataSet(new String[] {"P1","P2"}, new String[] {"N","C","X"}); // five samples, each data point has 2 items with three features

    @BeforeAll
    static void createDatasets() {
        for (Object[] sample : data6)
            ds6.addFeatureSample(sample);
        for (Object[][] sample : data2)
            ds2.addItemisedSample(sample);
    }

    @Test
    void toJSON() {
        JSONObject js1 = JSONUtils.toJSON(ds6);
        System.out.println(js1);
        JSONUtils.DataSet myds1 = JSONUtils.DataSet.fromJSON(js1);

        JSONObject js2 = JSONUtils.toJSON(ds2);
        System.out.println("Original: \t" + js2);
        JSONUtils.DataSet myds2 = JSONUtils.DataSet.fromJSON(js2);
        System.out.println("Recreated:\t" + JSONUtils.toJSON(myds2));
    }
}