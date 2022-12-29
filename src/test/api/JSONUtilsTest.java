package api;

import dat.file.TSVFile;
import json.JSONObject;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class JSONUtilsTest {

    JSONUtils.DataSet ds1 = new JSONUtils.DataSet();

    void createDataSet() {
        ds1.headers = new String[] {"F1", "F2"};
        ds1.values = new Object[][] {
                {1, "S1"},
                {2, "S2"},
                {null, "S3"},
                {4, null},
                {5, "S5"},
        };
    }

    @Test
    void toJSON() {
        createDataSet();
        JSONObject js1 = JSONUtils.toJSON(ds1);
        System.out.println(js1);
        JSONUtils.DataSet ds2 = JSONUtils.DataSet.fromJSON(js1);
        // TSVFile.print(ds2.values);
        for (int i = 0; i < ds2.values.length; i ++)
            for (int j = 0; j < ds2.values[i].length; j ++)
                if (ds1.values[i][j] != null)
                    assertTrue(ds1.values[i][j].equals(ds2.values[i][j]));
    }
}