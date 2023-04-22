package api;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import dat.Variable;
import dat.file.TSVFile;
import json.JSONArray;
import json.JSONObject;

public class JSONUtils {

    /**
     * A data set oriented around
     * - samples (repeated observations about the same items)
     * - items (objects which have features; items group features)
     * - features (fields describing each item)
     * Note that if an item has a single feature, an item is equivalent to a feature;
     * or, features may not be sensibly grouped, in which case one item implicitly encapsulates all features
     * (for backward compatibility these features are referred to as "headers")
     */
    public static class DataSet {
        private String[] features;
        private String[] items = null;

        // store values
        private Object[][][] values;
        // SETTING = ITEMISED;  values[sample][item][feature]
        // SETTING = NONITEMISED; values[0][sample][header]
        private Object[][] flattened = null; // storage of flat version of values (if computed)
        private String[] flattenedHeaders = null; // flat version of items and features (if computed)

        public DataSet(String[] items, String[] features) {
            this.features = features;
            this.items = items;
            this.values = new Object[0][items.length][features.length];
        }

        public DataSet(String[] headers) {
            this.features = headers;
            this.values = new Object[1][0][features.length];
        }

        public DataSet(String[] items, String[] features, Object[][][] values) {
            this(items, features);
            for (Object[][] sample : values)
                addItemisedSample(sample);
        }

        public DataSet(String[] headers, Object[][] values) {
            this(headers);
            for (Object[] sample : values)
                addFeatureSample(sample);
        }

        /**
         * In JSON, values must be one of the following data types:
         * a string
         * a number
         * an object (JSON object)
         * an array
         * a boolean
         * null
         *
         * To correct improperly assigned datatypes, this function will force values to adhere to variables that are recovered from a given BN.
         * For this to work, in the case of non-itemised data, features need to correspond to variable names;
         * in the case of itemised data, items are used as a prefix for features, e.g. item__feature.
         * Another pre-condition is that the datatype is defined in bn.Predef.
         *
         * @param bn Bayesian network with variables (with datatypes)
         */
        public void curateFeatures(BNet bn) {
            if (getNItems() > 0) {
                String[] headers = getFlattenedHeaders();
                Variable[][] vars = new Variable[getNItems()][headers.length];
                for (int i = 0; i < headers.length; i ++) {
                    int j = i / getNFeatures(); // index of item
                    int k = i % getNFeatures(); // index of feature
                    BNode node = bn.getNode(headers[i]);
                    if (node != null)
                        vars[j][k] = node.getVariable();
                }
                for (Object [][] sample : values) {
                    for (int j = 0; j < getNItems(); j ++) {
                        for (int k = 0; k < getNFeatures(); k ++) {
                            if (sample[j][k] != null) {
                                try {
                                    sample[j][k] = Predef.getObject(vars[j][k], (String) sample[j][k]);
                                } catch (ClassCastException e) {
                                    ; // not a string so nothing we can do
                                }
                            }
                        }
                    }
                }
            }
        }

        public String[] getFeatures() {
            return this.features;
        }

        public String[] getItems() {
            return this.items;
        }

        public Object[][][] getItemisedData() {
            return values;
        }

        public Object[][] getNonitemisedData() {
            if (isItemised())
                throw new JSONUtilsException("Data is itemised so cannot be accessed as non-itemised");
            return values[0];
        }

        public static Object[] flattenSample(Object[][] sample) {
            Object[] flatten = null;
            return flatten;
        }

        public Object[][] getFlattenedData() {
            if (getNItems() > 0) {
                if (flattened == null) {
                    flattened = new Object[values.length][getNItems() * getNFeatures()];
                    for (int i = 0; i < values.length; i++) {
                        Object[][] sample = values[i];
                        for (int j = 0; j < getNItems(); j++) {
                            for (int k = 0; k < getNFeatures(); k++)
                                flattened[i][getNFeatures() * j + k] = sample[j][k];
                        }
                    }
                }
                return flattened;
            }
            return getNonitemisedData();
        }

        public String[] getFlattenedHeaders() {
            if (getNItems() > 0) {
                if (flattenedHeaders == null) {
                    flattenedHeaders = new String[getNItems() * getNFeatures()];
                    for (int j = 0; j < getNItems(); j ++) {
                        for (int k = 0; k < getNFeatures(); k ++)
                            flattenedHeaders[getNFeatures() * j + k] = items[j] + "__" + features[k];
                    }
                }
                return flattenedHeaders;
            }
            return getFeatures();
        }

        /**
         * Check if this dataset has (multiple) items. If it has no items or a single item, there are only two dimensions of the dataset.
         * @return true, if the dataset has (multiple) items; false, otherwise.
         */
        public boolean isItemised() {
            if (this.items == null)
                return false;
            else if (this.items.length <= 1)
                return false;
            else
                return true;
        }

        public int getNFeatures() {
            return features == null ? 0 : features.length;
        }

        public int getNHeaders() {
            return getNFeatures();
        }

        public int getNItems() {
            return items == null ? 0 : items.length;
        }

        public int getNSamples() {
            if (getNItems() > 0)
                return values.length;
            else
                return values[0].length;
        }

        /**
         * Add a sample to the (non-itemised) dataset
         * @param sample
         * @return true if successful, else false
         */
        public boolean addFeatureSample(Object[] sample) {
            if (isItemised())
                return false;
            if (getNFeatures() != sample.length) // check number of features in the existing dataset against the new sample
                throw new JSONUtilsException("Invalid dataset: feature count does not match");
            Object[][] tmp = new Object[getNSamples() + 1][]; // extend allocation with ONE sample
            if (getNSamples() >= 0) System.arraycopy(this.values[0], 0, tmp, 0, getNSamples());
            tmp[getNSamples()] = sample;
            this.values[0] = tmp;
            return true;
        }

        /**
         * Add a sample to this itemised dataset
         * @param sample
         * @return true if successful, else false
         */
        public boolean addItemisedSample(Object[][] sample) {
            if (items == null)
                return false;
            if (getNItems() != sample.length) // check number of features in the existing dataset against the new sample
                throw new JSONUtilsException("Invalid dataset: item count does not match");
            if (getNFeatures() != sample[0].length) // check number of features in the existing dataset against the new sample
                throw new JSONUtilsException("Invalid dataset: feature count does not match");
            Object[][][] tmp = new Object[getNSamples() + 1][][]; // extend allocation with ONE sample
            if (getNSamples() >= 0) System.arraycopy(this.values, 0, tmp, 0, getNSamples());
            tmp[getNSamples()] = sample;
            this.values = tmp;
            return true;
        }

        public JSONObject toJSON() {
            if (getNItems() > 0)
                return JSONUtils.toJSON(items, features, values);
            else
                return JSONUtils.toJSON(features, values[0]);
        }

        public static DataSet fromJSON(JSONObject json) {
            DataSet ds = null;
            String[] items = null;
            String[] features = null;
            JSONArray jdata = json.optJSONArray("Data");
            if (jdata != null) {
                JSONArray jitems = json.optJSONArray("Items");
                if (jitems != null) {
                    items = new String[jitems.length()];
                    for (int i = 0; i < jitems.length(); i ++)
                        items[i] = jitems.getString(i);
                } else { // no items, so possibly we have "headers" (equiv to features in the absence of items)
                    JSONArray jheaders = json.optJSONArray("Headers");
                    if (jheaders != null) {
                        features = new String[jheaders.length()];
                        for (int i = 0; i < jheaders.length(); i ++)
                            features[i] = jheaders.getString(i);
                        ds = new DataSet(features);
                    }
                }
                if (features == null) {
                    JSONArray jfeatures = json.optJSONArray("Features");
                    if (jfeatures != null) {
                        features = new String[jfeatures.length()];
                        for (int i = 0; i < jfeatures.length(); i ++)
                            features[i] = jfeatures.getString(i);
                    } else // no headers or features given
                        throw new JSONUtilsException("Features or headers are not given");
                    ds = items == null ? new DataSet(features) : new DataSet(items, features);
                }
                // now lets process the data
                JSONArray[] jsamples = new JSONArray[jdata.length()];
                for (int i = 0; i < jdata.length(); i++) {
                    jsamples[i] = jdata.getJSONArray(i);
                    if (i > 0)
                        if (jsamples[i].length() != jsamples[i-1].length()) // check number of data items
                            throw new JSONUtilsException("Invalid data matrix: data items in columns do not match");
                }
                // unfolding the data set means sample by sample, then either:
                // 1. item by item, then feature-value by feature-value, ... or IF item count is one or none are given...
                // 2. feature-value by feature-value
                boolean OPT1 = true;
                boolean OPT2 = true;
                for (int i = 0; i < jsamples.length; i++) { // sample i
                    Object[] values = new Object[jsamples[i].length()];
                    Object[][] values_opt1 = new Object[jsamples[i].length()][];
                    for (int j = 0; j < jsamples[i].length(); j ++) { // iterate through elements (j)  of sample i
                        // elements can be either option 1, e.g. {1,"A"} in {{1,"A"},{2,"B"},{3,"C"},{4,"D"}}, or option 2, e.g. 1 in {1,2,3,4}
                        if (jsamples[i].isNull(j)) {
                            OPT1 = false;
                            values[j] = null;
                        } else {
                            try {
                                JSONArray jarr = (JSONArray) jsamples[i].get(j); // if cast succeeds it is option 1
                                OPT2 = false;
                                Object[] opt1 = new Object[jarr.length()];
                                for (int k = 0; k < opt1.length; k++) {
                                    if (jarr.isNull(k))
                                        opt1[k] = null;
                                    else
                                        opt1[k] = jarr.get(k);
                                }
                                values_opt1[j] = opt1;
                            } catch (ClassCastException e) {
                                // dataset can still be itemised if only single feature, but typically not itemised and multiple feature values
                                if (ds.isItemised() && ds.getNFeatures() == 1 || !ds.isItemised() && ds.getNFeatures() == values.length) {
                                    OPT1 = false;
                                    values[j] = jsamples[i].get(j);
                                    if (values[j] == JSONObject.NULL)
                                        values[j] = null;
                                } else
                                    throw new JSONUtilsException("Invalid data matrix: data items not consistent with items and/or features given @ sample " + i + " column " + j);
                            }
                        }
                    }
                    if (OPT1)
                        ds.addItemisedSample(values_opt1);
                    else if (OPT2)
                        ds.addFeatureSample(values);
                    else
                        throw new JSONUtilsException("Invalid data matrix: data items are not consistently provided");
                }
            } else
                return null;
            return ds;
        }
    }

    public static JSONObject toJSON(DataSet dataset) {
        if (dataset.getNItems() > 0)
            return JSONUtils.toJSON(dataset.getItems(), dataset.getFeatures(), dataset.getItemisedData());
        else
            return JSONUtils.toJSON(dataset.getFeatures(), dataset.getNonitemisedData());
    }

    public static JSONObject toJSON(TSVFile tsv) {
        return toJSON(tsv.getHeaders(), tsv.getRows());
    }

    /**
     * Transform headers and data matrix into a JSON representation; the matrix
     * has rows that correspond to samples and columns to features (named by headers);
     * the resulting JSON representation joins all samples for each feature/header (i.e. the transpose of the matrix).
     * @param items
     * @param features
     * @param data matrix ([sample][item][feature])
     * @return
     */
    public static JSONObject toJSON(String[] items, String[] features, Object[][][] data) {
        JSONObject json = new JSONObject();
        json.put("Items", new JSONArray(items));
        json.put("Features", new JSONArray(features));
        JSONArray values = new JSONArray();
        for (int i = 0; i < data.length; i ++) {
            JSONArray sample = new JSONArray();
            for (int j = 0; j < data[i].length; j ++) {
                JSONArray col = new JSONArray();
                for (int k = 0; k < data[i][j].length; k ++)
                    col.put(data[i][j][k]);
                sample.put(col);
            }
            values.put(sample);
        }
        json.put("Data", values);
        return json;
    }

    /**
     * Transform headers and data matrix into a JSON representation; the matrix
     * has rows that correspond to samples and columns to features (named by headers);
     * the resulting JSON representation joins all samples for each feature/header (i.e. the transpose of the matrix).
     * @param headers
     * @param data matrix ([sample][feature])
     * @return
     */
    public static JSONObject toJSON(String[] headers, Object[][] data) {
        JSONObject json = new JSONObject();
        json.put("Headers", new JSONArray(headers));
        JSONArray values = new JSONArray();
        for (int i = 0; i < data.length; i ++) {
            JSONArray sample = new JSONArray();
            for (int j = 0; j < data[i].length; j ++)
                sample.put(data[i][j]);
            values.put(sample);
        }
        json.put("Data", values);
        return json;
    }


    public static class JSONUtilsException extends RuntimeException {
        public JSONUtilsException(String msg) {
            super(msg);
        }
    }
}
