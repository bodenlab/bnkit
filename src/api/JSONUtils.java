package api;

import dat.file.TSVFile;
import json.JSONArray;
import json.JSONObject;

public class JSONUtils {

    public static class DataSet {
        String[] headers;
        Object[][] values; // [row][col], but note that the JSON array is first column, then row

        public static DataSet fromJSON(JSONObject json) {
            DataSet ds = new DataSet();
            JSONArray jdata = json.optJSONArray("Data");
            if (jdata != null) {
                JSONArray jheaders = json.optJSONArray("Headers");
                if (jheaders != null) {
                    if (jheaders.length() != jdata.length())
                        throw new JSONUtilsException("Number of headers and data items do not match");
                }
                ds.headers = new String[jdata.length()];
                JSONArray[] jcols = new JSONArray[jdata.length()];
                for (int i = 0; i < jdata.length(); i++) {
                    ds.headers[i] = jheaders == null ? "X" + (i + 1) : jheaders.getString(i);
                    jcols[i] = jdata.getJSONArray(i);
                    if (i > 0)
                        if (jcols[i].length() != jcols[i-1].length()) // check number of data items
                            throw new JSONUtilsException("Invalid data matrix: data items in columns do not match");
                }
                ds.values = new Object[jcols[0].length()][jdata.length()];
                for (int i = 0; i < jcols.length; i++) {
                    for (int j = 0; j < jcols[0].length(); j ++) {
                        ds.values[j][i] = jcols[i].get(j);
                        if (ds.values[j][i] == JSONObject.NULL)
                            ds.values[j][i] = null;
                    }
                }
            } else
                return null;
            return ds;
        }

        public static DataSet fromJSON(JSONArray jdata) {
            DataSet ds = new DataSet();
            if (jdata != null) {
                ds.headers = new String[jdata.length()];
                JSONArray[] jcols = new JSONArray[jdata.length()];
                for (int i = 0; i < jdata.length(); i++) {
                    ds.headers[i] = "X" + (i + 1);
                    jcols[i] = jdata.getJSONArray(i);
                    if (i > 0)
                        if (jcols[i].length() != jcols[i-1].length()) // check number of data items
                            throw new JSONUtilsException("Invalid data matrix: data items in columns do not match");
                }
                ds.values = new Object[jcols[0].length()][jdata.length()];
                for (int i = 0; i < jcols.length; i++) {
                    for (int j = 0; j < jcols[0].length(); j ++) {
                        ds.values[j][i] = jcols[i].get(j);
                        if (ds.values[j][i] == JSONObject.NULL)
                            ds.values[j][i] = null;
                    }
                }
            }
            return ds;
        }
    }

    public static JSONObject toJSON(DataSet dataset) {
        return JSONUtils.toJSON(dataset.headers, dataset.values);
    }

    public static JSONObject toJSON(TSVFile tsv) {
        return toJSON(tsv.getHeaders(), tsv.getRows());
    }

    public static JSONObject toJSON(String[] headers, Object[][] data) {
        JSONObject json = new JSONObject();
        json.put("Headers", headers);
        JSONArray values = new JSONArray();
        for (int i = 0; i < headers.length; i ++) {
            Object[] col = new Object[data.length];
            for (int j = 0; j < data.length; j ++) {
                if (i >= data[j].length)
                    throw new JSONUtilsException("Invalid data matrix");
                col[j] = data[j][i];
            }
            JSONArray arr = new JSONArray(col);
            values.put(arr);
        }
        json.put("Data", values);
        return json;
    }

    public static JSONArray toJSON(Object[][] data) {
        JSONArray jsonArray = new JSONArray();
        int nfeats = 0;
        for (int j = 0; j < data.length; j ++) {
            nfeats = j == 0 ? data[j].length : nfeats;
            if (nfeats != data[j].length)
                throw new JSONUtilsException("Invalid data matrix");
        }
        for (int i = 0; i < nfeats; i ++) {
            Object[] col = new Object[data.length];
            for (int j = 0; j < data.length; j ++)
                col[j] = data[i][j];
            jsonArray.put(new JSONArray(col));
        }
        return jsonArray;
    }


    public static class JSONUtilsException extends RuntimeException {
        public JSONUtilsException(String msg) {
            super(msg);
        }
    }
}
