package dat.file;

import json.JSONObject;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility class to work with tabulated files, sometimes labeled, sometimes not.
 * Cells can contain arbitrary datatypes, though different applications may have certain requirements
 * and this class is featured to deal with some of them.
 *
 * Created by mikael on 2/08/2016. Cloned from binfkit March 2019.
 * Amended and made specific to bnkit from then on.
 * @author mikael
 */
public class TSVFile {

    /**
     * NULLS identify all string tokens that are interpreted as NULL/missing values (in addition to the empty string) when loading an Object matrix from a TSV file
     */
    public static String[] NULLS = {"null", "NULL", "nil", "NIL", "none", "None", "NONE"};

    List<Object[]> rows = new ArrayList<>();
    int ncols = 0;
    Map<String, Integer> headers = null;
    Map<Integer, Map<Object, int[]>> indexMap = new HashMap<>();

    /**
     * Read TSV from a file.
     * @param filename name of file
     * @param useHeader interpret the first row as headers
     * @param transpose transpose the cells before evaluating values in cells (including headers)
     * @throws IOException
     */
    public TSVFile(String filename, boolean useHeader, boolean transpose) throws IOException {
        this(transpose?Transpose(TSVFile.loadObjects(filename)):TSVFile.loadObjects(filename), useHeader);
    }

    /**
     * Read TSV from a file.
     * @param filename name of file
     * @param useHeader interpret the first row as headers
     * @throws IOException
     */
    public TSVFile(String filename, boolean useHeader) throws IOException {
        this(TSVFile.loadObjects(filename), useHeader);
    }

    /**
     * Read TSV file from standard input
     * @param useHeader interpret the first row as headers
     * @throws IOException
     */
    public TSVFile(boolean useHeader) throws IOException {
        this(TSVFile.loadObjects(new BufferedReader(new InputStreamReader(System.in))));
    }

    /**
     * Construct a TSV data structure from a matrix of objects, which has a header row (first row)
     * @param headers the headers (one for each column)
     * @param objects the values in the matrix, index by row then column
     */
    public TSVFile(String[] headers, Object[][] objects) {
        ncols = headers.length;
        this.headers = new HashMap<>();
        for (int i = 0; i < headers.length; i++)
            this.headers.put(headers[i], i);
        for (int i = 0; i < objects.length; i++)
            if (objects[i] != null)
                rows.add(i, objects[i]);
    }

    /**
     * Construct a TSV data structure from a matrix of objects, which has a header row (first row)
     * @param objects the values in the matrix, index by row then column
     * @param useHeader indicate if the header row is used to label columns or ignored
     */
    public TSVFile(Object[][] objects, boolean useHeader) {
        if (objects.length > 0) {
            ncols = objects[0].length;
            if (useHeader) {
                headers = new HashMap<>();
                for (int i = 0; i < objects[0].length; i++) {
                    try {
                        headers.put(objects[0][i] == null ? null : objects[0][i].toString(), i);
                    } catch (ClassCastException e) {
                        throw new RuntimeException("Failed to index header: not a text string \"" + objects[0][i] + "\"");
                    }
                }
            }
            int index = 0;
            for (int i = (useHeader ? 1 : 0); i < objects.length; i++) { // start at 1 and skip header conditionally
                if (objects[i] != null)
                    rows.add(index ++, objects[i]);
            }
        }
    }

    /**
     * Construct a TSV data structure from a matrix of objects, with all values (no header)
     * @param objects the values in the matrix, index by row then column
     */
    public TSVFile(Object[][] objects) {
        if (objects.length > 0) {
            ncols = objects[0].length;
            for (int i = 0; i < objects.length; i++) { // start at 0, there is no header
                int index = i;
                if (objects[i] != null)
                    rows.add(index, objects[i]);
            }
        }
    }

    /**
     * Return transpose of a specified matrix of objects
     * @param orig original matrix (not modified)
     * @return transposition of specified matrix
     */
    public static Object[][] Transpose(Object[][] orig) {
        int nrows = orig.length;
        int ncols = -1;
        for (int r = 0; r < nrows; r ++) {
            if (ncols != -1 && ncols != orig[r].length)
                throw new RuntimeException("Not a valid matrix");
            ncols = orig[r].length;
        }
        Object[][] transp = new Object[ncols][nrows];
        for (int r = 0; r < nrows; r ++) {
            for (int c = 0; c < ncols; c ++)
                transp[c][r] = orig[r][c];
        }
        return transp;
    }

    /**
     * Concatenates two matrices of objects horizontally.
     * The matrices must have the same number of rows.
     *
     * @param a the first matrix to concatenate
     * @param b the second matrix to concatenate
     * @return a new matrix that is the result of concatenating the two input matrices
     * @throws RuntimeException if the matrices do not have the same number of rows
     */
    public static Object[][] concat(Object[][] a, Object[][] b) {
        if (a.length != b.length)
            throw new RuntimeException("Matrices must have the same number of rows");
        Object[][] c = new Object[a.length][a[0].length + b[0].length];
        for (int i = 0; i < a.length; i ++) {
            for (int j = 0; j < a[0].length; j ++)
                c[i][j] = a[i][j];
            for (int j = 0; j < b[0].length; j ++)
                c[i][j + a[0].length] = b[i][j];
        }
        return c;
    }

    public static TSVFile concat(TSVFile a, TSVFile b) {
        String[] headers = new String[a.ncols + b.ncols];
        String[] aheaders = a.getHeaders();
        String[] bheaders = b.getHeaders();
        for (int i = 0; i < headers.length; i ++)
            headers[i] = i < a.ncols ? aheaders[i] : bheaders[i - a.ncols];
        return new TSVFile(headers, concat(a.getRows(), b.getRows()));
    }
    public static Object[] tokeniseBySep(Object orig, String separator) {
        // tokenise the semi-colon separated fields within the entry
        List<String> tokens = new ArrayList<>();
        StringTokenizer tokenizer = new StringTokenizer(orig.toString(), separator);
        while (tokenizer.hasMoreTokens())
            tokens.add(tokenizer.nextToken());
        Object[] tokarr = new Object[tokens.size()];
        for (int i = 0; i < tokens.size(); i ++) {
            Object y = tokens.get(i);
            try {
                y = Double.parseDouble(tokens.get(i));
            } catch (NumberFormatException e1) {
                try {
                    y = Integer.parseInt(tokens.get(i));
                } catch (NumberFormatException e2) {
                }
            }
            tokarr[i] = y;
        }
        return tokarr;
    }

    public static Object[] tokeniseRange(Object orig) {
        // the observed value is occasionally a "range" with min/max values as a pair separated by hyphen
        Object[] range = tokeniseBySep(orig, "-");
        if (range.length >= 2)
            return range;
        else
            return null;
    }

    public static Object mergeRange(Object[] range) {
        Object[] tokarr = new Object[range.length];
        StringBuilder restring = new StringBuilder();
        Double dsum = null;
        int dcnt = 0;
        Integer isum = null;
        int icnt = 0;
        boolean disqualify = false;
        for (int i = 0; i < range.length; i ++) {
            Object y = range[i];
            restring.append(range[i].toString());
            try {
                y = (Double) range[i];
                if (dsum == null)
                    dsum = (Double) y;
                else
                    dsum += (Double) y;
                dcnt += 1;
            } catch (ClassCastException e1) {
                try {
                    y = (Integer) range[i];
                    if (isum == null)
                        isum = (Integer) y;
                    else
                        isum += (Integer) y;
                    icnt += 1;
                } catch (ClassCastException e2) {
                    disqualify = true;
                }
            }
            restring.append("-");
        }
        if (disqualify)
            return restring.toString();
        if (dsum != null)
            return dsum / dcnt;
        if (isum != null)
            return isum / icnt;
        return restring.toString();
    }

    public static Object[] tokeniseByCount(Object orig) {
        String search = "_count=";
        String line = orig.toString();
        int cidx = line.indexOf(search);
        List<Object> toklst = new ArrayList<>();
        if (cidx < 0) {
            // there is not any "_count" tags...
            cidx = line.length();
        }
        while (cidx >= 0) {
            String token = line.substring(0, cidx);
            Object[] separates = tokeniseBySep(token, ";");
            // ignore everything but the first element
            Object[] range = tokeniseRange(separates[0]);
            if (range == null)
                toklst.add(separates[0]);
            else
                toklst.add(mergeRange(range));
            line = line.substring(Math.min(cidx + search.length(), line.length()));
            int eidx = line.indexOf(";");
            if (eidx == -1) { // last count token reached
                // collect count number?
                break;
            } else { // semi-colon found
                // collect count number?
                line = line.substring(eidx + 1);
            }
            cidx = line.indexOf(search);
        }
        Object[] ret = new Object[toklst.size()];
        toklst.toArray(ret);
        return ret;
    }

    /**
     * Parse the cell as a text string specifying value/values on the BRENDA/Foley format
     * (as dictated by Gabe Foley's ASR curation pipeline format incorporating data from BRENDA)
     * @param orig the original object
     * @return the parsed object with datatype if specified
     * @throws ClassCastException if the parsing fails
    Examples:
        0.13;4-nitrophenyl alpha-D-maltoheptaoside-4,6-O-ethylidene_count=1
        377.0;starch_count=1;619.0;starch_count=2;880.0;starch_count=3
        7.2-7.5_count=1
        34_count=1;55_count=3
        SM00642;SM00632;
     */
//    public static Object tokeniseBRENDA(Object orig) throws ClassCastException {
//    }


    /**
     * Save the present values to a TSV file.
     * Headers are used if they have been set when constructed.
     * @param filename name of file
     * @throws IOException if the IO fails
     */
    public void save(String filename) throws IOException {
        Object[][] ret;
        Object[][] rows = getRows();
        if (headers != null) {
            ret = new Object[rows.length + 1][];
            ret[0] = new String[headers.size()];
            for (Map.Entry<String, Integer> entry : headers.entrySet())
                ret[0][entry.getValue()] = entry.getKey();
            for (int i = 0; i < rows.length; i ++)
                ret[i + 1] = rows[i];
        } else
            ret = rows;
        saveObjects(filename, ret);
    }

    /**
     * Save the present values as a TSV file on the standard output.
     * Headers are used if they have been set when constructed.
     * @throws IOException
     */
    public void save() throws IOException {
        Object[][] ret;
        Object[][] rows = getRows();
        if (headers != null) {
            ret = new Object[rows.length + 1][];
            ret[0] = new String[headers.size()];
            for (Map.Entry<String, Integer> entry : headers.entrySet())
                ret[0][entry.getValue()] = entry.getKey();
            for (int i = 0; i < rows.length; i ++)
                ret[i + 1] = rows[i];
        } else
            ret = rows;
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
        saveObjects(bw, ret);
        bw.close();
    }

    /**
     * Retrieve the column index of a given header
     * @param header name of column
     * @return the column index (starting at 0); returns -1 if not found
     */
    public int getColumn(String header) {
        Integer col = headers.get(header);
        if (col == null)
            return -1;
        return col;
    }

    public String getHeader(int col) {
        for (Map.Entry<String, Integer> entry : headers.entrySet())
            if (entry.getValue().intValue() == col)
                return entry.getKey();
        return null;
    }

    public String[] getHeaders() {
        String[] s = new String[headers.size()];
        for (Map.Entry<String, Integer> entry : headers.entrySet())
            s[entry.getValue().intValue()] = entry.getKey();
        return s;
    }

    public Object getValue(Object[] row, String colname) {
        Integer col = headers.get(colname);
        if (col == null)
            throw new RuntimeException("Invalid column name: " + colname);
        return row[col];
    }

    public Object getValue(int rowidx, String colname) {
        return getValue(rows.get(rowidx), colname);
    }

    public Set<Object> getValues(String colname) {
        Integer col = headers.get(colname);
        if (col == null)
            throw new RuntimeException("Invalid column name: " + colname);
        return getValues(col);
    }

    public Set<Object> getValues(int col) {
        Map<Object, int[]> keyIndex = getIndexMap(col);
        return keyIndex.keySet();
    }

    public Set<Object> getValues() {
        Set<Object> vals = new HashSet<>();
        for (Object[] row : rows)
            for (Object val : row)
                if (val != null)
                    vals.add(val);
        return vals;
    }

    /**
     * Check if all values are either cast-able to double or null.
     * @param col
     * @return
     */
    public static boolean isDouble(Object[] col) {
        try {
            for (int i = 0; i < col.length; i ++) {
                Double y = (Double) col[i];
            }
            return true; // must have all been of type T, because no exception has been thrown
        } catch (ClassCastException e) {
            return false;
        }
    }

    public static boolean isDoubleOrInt(Object[] col) {
        for (int i = 0; i < col.length; i ++) {

            try {
                Double x = (Double) col[i];
            } catch (ClassCastException e) {
                try {
                    Integer y = (Integer) col[i];
                } catch (ClassCastException e2) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Check if all values are either cast-able to double or null.
     * @param rows
     * @return
     */
    public static boolean isDouble(Object[][] rows) {
        try {
            for (int i = 0; i < rows.length; i ++) {
                for (int j = 0; j < rows[i].length; j ++) {
                    Double y = (Double) rows[i][j];
                }
            }
            return true; // must have all been of type T, because no exception has been thrown
        } catch (ClassCastException e) {
            return false;
        }
    }

    /**
     * Check if all values are either cast-able to double, int or null.
     * @param rows
     * @return
     */
    public static boolean isDoubleOrInt(Object[][] rows) {
        for (int i = 0; i < rows.length; i ++) {
            for (int j = 0; j < rows[i].length; j ++) {
                try {
                    Double y = (Double) rows[i][j];
                } catch (ClassCastException e) {
                    try {
                        Integer y = (Integer) rows[i][j];
                    } catch (ClassCastException e2) {
                        return false;
                    }
                }
            }
        }
        return true; // must have all been of type T, because no exception has been thrown
    }

    private Map<Object, int[]> getIndexMap(int col) {
        Map<Object, int[]> keyIndex = indexMap.get(col);
        if (keyIndex == null) { // we need to construct it
            Map<Object, List<Integer>> tmpIndex = new HashMap<>();
            for (int i = 0; i < rows.size(); i ++) {
                List<Integer> tmpIndices = tmpIndex.get(rows.get(i)[col]);
                if (tmpIndices != null)
                    tmpIndices.add(i);
                else {
                    tmpIndices = new ArrayList<>();
                    tmpIndices.add(i);
                    tmpIndex.put(rows.get(i)[col], tmpIndices);
                }
            }
            keyIndex = new HashMap<>();
            for (Map.Entry<Object, List<Integer>> entry : tmpIndex.entrySet()) {
                List<Integer> indices = entry.getValue();
                int[] arr = new int[indices.size()];
                for (int i = 0; i < arr.length; i ++)
                    arr[i] = indices.get(i);
                keyIndex.put(entry.getKey(), arr);
            }
            indexMap.put(col, keyIndex);
        }
        return keyIndex;
    }

    /**
     * Retrieve the list of indices for rows, for which the specified column has a given value.
     * All rows were initially read from a datafile, but only indexed the first time the column has been referred to.
     * @param key value of column
     * @param colname column name
     * @return an array of indices, for rows with that value in the column; null if the value is not found
     */
    public int[] getIndices(Object key, String colname) {
        Integer col = headers.get(colname);
        if (col == null)
            throw new RuntimeException("Invalid column name: " + colname);
        return getIndices(key, col);
    }

    /**
     * Retrieve the list of indices for rows, for which the specified column has a given value.
     * All rows were initially read from a datafile, but only indexed the first time the column has been referred to.
     * @param key value of column
     * @param col column index
     * @return an array of indices, for rows with that value in the column; null if the value is not found
     */
    public int[] getIndices(Object key, int col) {
        Map<Object, int[]> keyIndex = getIndexMap(col);
        return keyIndex.get(key);
    }

    public Object[] getRow(int index) {
        return rows.get(index);
    }

    public Object[][] getRows(int[] indices) {
        Object[][] ret = new Object[indices.length][];
        for (int i = 0; i < indices.length; i ++)
            ret[i] = getRow(indices[i]);
        return ret;
    }

    public Object[][] getRows() {
        Object[][] ret = new Object[rows.size()][];
        for (int i = 0; i < rows.size(); i ++)
            ret[i] = getRow(i);
        return ret;
    }

    /**
     * Retrieve whole column by index
     * @param column_index the index of the column
     * @return values in specified column in order of row
     */
    public Object[] getCol(int column_index) {
        Object[] ret = new Object[rows.size()];
        for (int i = 0; i < rows.size(); i ++)
            ret[i] = getRow(i)[column_index];
        return ret;
    }

    /**
     * Retrieve whole column by index
     * @param column_index the index of the column
     * @param parser parser (e.g. Double, BRENDA)
     * @return values in specified column in order of row, that all conform to type
     */
    public Object[] getCol(int column_index, String parser) {
        Object[] ret = new Object[rows.size()];
        for (int i = 0; i < rows.size(); i ++) {
            ret[i] = getRow(i)[column_index];
            if (ret[i] != null) {
                if (parser != null) {
                    if (parser.equals("BRENDA")) {
                        Object[] tokens = tokeniseByCount(ret[i]);
                        ret[i] = tokens[0];
                    } else if (parser.equals("Double")) {
                        try {
                            ret[i] = Double.parseDouble(ret[i].toString());
                        } catch (NumberFormatException e1) {
                            try {
                                ret[i] = (double) Integer.parseInt(ret[i].toString());
                            } catch (NumberFormatException e2) {
                                ret[i] = null;
                            }
                        }
                    } else {
                        // unknown parser
                    }
                }
            }
        }
        return ret;
    }

    /**
     * Retrieve whole columns by index;
     * this is the matrix "transposed"
     * @param column_indices the indices of the columns
     * @return values in specified columns in order of row
     */
    public Object[][] getCols(int[] column_indices) {
        Object[][] ret = new Object[column_indices.length][rows.size()];
        for (int i = 0; i < rows.size(); i ++) {
            for (int j = 0; j < column_indices.length; j ++)
                ret[j][i] = rows.get(i)[column_indices[j]];
        }
        return ret;
    }

    public static void print(Object[][] rows) {
        for (int i = 0; i < rows.length; i ++) {
            for (int j = 0; j < rows[i].length; j ++) {
                System.out.print(rows[i][j] + "\t");
            }
            System.out.println();
        }
    }

    /**
     * Load a TSV file from standard input and place the contents in a two-dimensional object matrix.
     * @return
     * @throws IOException
     */
    public static Object[][] loadObjects(BufferedReader br) throws IOException {
        Object[][] data = null;
        String line = br.readLine();
        List<Object[]> alldata = new ArrayList<>();
        while (line != null) {
            if (!line.startsWith("#")) {
                String[] tokens = line.split("\t");
                Object[] values = new Object[tokens.length];
                for (int i = 0; i < tokens.length; i++) {
                    try {
                        values[i] = Integer.valueOf(tokens[i]); // value is an int
                    } catch (NumberFormatException e1) {
                        try {
                            values[i] = Double.valueOf(tokens[i]); // value is a double
                        } catch (NumberFormatException e2) {
                            boolean isnull = false;
                            for (String NULL : NULLS) {
                                if (tokens[i].equals(NULL)) {
                                    isnull = true;
                                    break;
                                }
                            }
                            if (isnull || tokens[i].isBlank())
                                values[i] = null;
                            else
                                values[i] = tokens[i]; // value is a string
                        }
                    }
                }
                alldata.add(values);
            }
            line = br.readLine();
        }
        data = new Object[alldata.size()][];
        int ncols = -1; // figure out how many cols we need
        for (int k = 0; k < data.length; k++) {
            if (ncols == -1)
                ncols = alldata.get(k).length;
            else if (ncols != alldata.get(k).length) {
                // could throw error, or just allocate the exact number;
                // but, we ignore to make sure the data is a proper matrix
            }
            data[k] = new Object[ncols];
            for (int j = 0; j < ncols; j++) {
                data[k][j] = (j < alldata.get(k).length ? alldata.get(k)[j] : null);
            }
        }
        br.close();
        return data;
    }


    /**
     * Load a TSV file from a file and place the contents in a two-dimensional object matrix.
     * Values will either be Integer, Double or String.
     * @param filename
     * @return
     * @throws IOException
     */
    public static Object[][] loadObjects(String filename) throws IOException {
        BufferedReader br = null;
        br = new BufferedReader(new FileReader(filename));
        return loadObjects(br);
    }

    /**
     * Save matrix as a TSV file via a given writer handle.
     * For standard output pass this:
     * new BufferedWriter(new OutputStreamWriter(System.out))
     * @param bd buffered writer
     * @param matrix
     * @throws IOException
     */
    public static void saveObjects(BufferedWriter bd, Object[][] matrix) throws IOException {
        for (Object[] row : matrix) {
            if (row == null)
                System.err.println("row is null");
            for (int i = 0; i < row.length; i ++) {
                Object val = row[i] != null ? row[i] : "";
                bd.write(val + (i == row.length - 1 ? "" : "\t"));
            }
            bd.newLine();
        }
        bd.flush();
    }

    /**
     * Save matrix as a TSV file.
     * @param filename name of file
     * @param matrix
     * @throws IOException
     */
    public static void saveObjects(String filename, Object[][] matrix) throws IOException {
        BufferedWriter bd = new BufferedWriter(new FileWriter(filename));
        saveObjects(bd, matrix);
        bd.close();
    }

    public static Integer DEFAULT_SHAPE = 2; // #1: rectangle  #2: circle  #3: star  #4: right pointing triangle  #5: left pointing triangle  #6: checkmark
    public static Integer DEFAULT_SIZE = 3; // #size can be any number. Maximum size in the dataset will be displayed using MAXIMUM_SIZE, while others will be proportionally smaller
    public static Double DEFAULT_POS = 1.0; // #position is a number between 0 and 1 and defines the position of the symbol on the branch (0 is at the start of node branch, 0.5 is in the middle, and 1 is at the end)
    public static Integer DEFAULT_FILL = 1; // #fill can be 1 or 0. If set to 0, only the outline of the symbol will be displayed.

    /**
     * By Sasha Trubetskoy, from https://sashamaps.net/docs/resources/20-colors/ excluding white
     */
    public static String[] NONWHITE_COLORS = {"#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#000000"};

    public static Double getMin(Object[] dblvals) {
        Double min = null;
        for (int i = 0; i < dblvals.length; i ++) {
            if (dblvals[i] == null)
                continue;
            if (min == null)
                min = (Double) dblvals[i];
            else if (min > (Double) dblvals[i])
                min = (Double) dblvals[i];
        }
        return min;
    }
    public static Double getMax(Object[] dblvals) {
        Double max = null;
        for (int i = 0; i < dblvals.length; i ++) {
            if (dblvals[i] == null)
                continue;
            if (max == null)
                max = (Double) dblvals[i];
            else if (max < (Double) dblvals[i])
                max = (Double) dblvals[i];
        }
        return max;
    }

    public static void save2iTOL(String filename, Object[] items, Object[] values, String dataset_label, int nbins, Double setmin, Double setmax) throws IOException {
        save2iTOL(filename, items, values, null, dataset_label, nbins, setmin, setmax);
    }

    public static void save2iTOL(String filename, Object[] items, Object[] values, Object[] confid, String dataset_label, int nbins, Double setmin, Double setmax) throws IOException {
        BufferedWriter bd = new BufferedWriter(new FileWriter(filename));
        bd.write("DATASET_SYMBOL"); bd.newLine();
        bd.write("SEPARATOR SPACE"); bd.newLine();
        bd.write("DATASET_LABEL " + dataset_label); bd.newLine();
        bd.write("COLOR #ffff00"); bd.newLine();
        if (isDouble(values)) { // values are doubles
            Double min = null, max = null;
            if (setmin != null)
                min = setmin;
            else
                min = getMin(values);
            if (setmax != null)
                max = setmax;
            else
                max = getMax(values);
            double binrange = (max - min) / nbins;
            double[] bins = new double[nbins];
            String[] bhex = new String[nbins];
            int hexrange = 255 / (nbins - 1);
            for (int i = 0; i < nbins; i ++) {
                bins[i] = min + binrange / 2 + binrange * i;
                bhex[i] = String.format("#" + "%02x", Math.min(2 * (hexrange * i), 255)) + String.format("%02x", Math.max(2 * (hexrange * i) - 255, 0)) + String.format("%02x", Math.max(0, 255 - hexrange * i));
                // System.out.println(i +"\t" + bhex[i] + "\t" + bins[i]);
            }
            int[] USE_SIZE = new int[values.length];
            Arrays.fill(USE_SIZE, DEFAULT_SIZE);
            if (confid != null) { // there are values for expressing confidence too;
                if (isDouble(confid)) { // confid values are doubles; presently, standard deviation
                    Double cmin = null, cmax = null;
                    cmin = getMin(confid);
                    cmax = getMax(confid);
                    double cbinrange = (cmax - cmin) / 3;
                    double[] cbins = new double[3];
                    Integer[] cbint = new Integer[3]; // the SIZE to use, e.g. 3, 2 and 1
                    for (int i = 0; i < 3; i++) {
                        cbint[i] = Math.min(3, DEFAULT_SIZE) - i;
                        cbins[i] = cmin + cbinrange / 2 + cbinrange * i; // middle value
                    }
                    for (int i = 0; i < confid.length; i++) {
                        if (confid[i] == null)
                            USE_SIZE[i] = DEFAULT_SIZE;
                        else {
                            int cbin = Math.min(3 - 1, (int) (((Double) confid[i] + cbinrange / 2 - cmin) / cbinrange));
                            USE_SIZE[i] = cbint[cbin];
                        }
                    }
                }
            }
            /*
            DATASET_SYMBOL
            SEPARATOR SPACE
            #label is used in the legend table (can be changed later)
            DATASET_LABEL example symbols
            #dataset color (can be changed later)
            COLOR #ffff00
            LEGEND_TITLE Tm
            LEGEND_SHAPES 2 2 2
            LEGEND_COLORS #ff0000 #880088 #0000ff
            LEGEND_LABELS 60 45 30
            #largest symbol will be displayed with this size, others will be proportionally smaller.
            MAXIMUM_SIZE 10
            #ID symbol size,color fill position label
            #symbol should be a number between 1 and 5: #1: rectangle  #2: circle  #3: star  #4: right pointing triangle  #5: left pointing triangle  #6: checkmark
            #size can be any number. Maximum size in the dataset will be displayed using MAXIMUM_SIZE, while others will be proportionally smaller
            #color can be in hexadecimal, RGB or RGBA notation. If RGB or RGBA are used, dataset SEPARATOR cannot be comma.
            #fill can be 1 or 0. If set to 0, only the outline of the symbol will be displayed.
            #position is a number between 0 and 1 and defines the position of the symbol on the branch (for example, position 0 is exactly at the start of node branch, position 0.5 is in the middle, and position 1 is at the end)
            DATA
            ASR01 2 3 #880088 1 1
            */
            bd.write("LEGEND_TITLE " + dataset_label); bd.newLine();
            StringBuffer shapes = new StringBuffer("LEGEND_SHAPES ");
            StringBuffer colors = new StringBuffer("LEGEND_COLORS ");
            StringBuffer labels = new StringBuffer("LEGEND_LABELS ");
            for (int i = 0; i < nbins; i ++) {
                shapes.append(DEFAULT_SHAPE + " ");
                colors.append(bhex[i] + " ");
                // FIXME: figure out suitable number of decimal places,
                //  e.g. range is 0 - 100 and 10 bins, 0 decimals, range is -3 to -2.5 and 10 bins, 2 decimals
                labels.append(String.format("%.2f ", bins[i]));
            }
            bd.write(shapes.toString()); bd.newLine();
            bd.write(colors.toString()); bd.newLine();
            bd.write(labels.toString()); bd.newLine();
            bd.write("MAXIMUM_SIZE 10"); bd.newLine();
            bd.write("DATA"); bd.newLine();
            for (int i = 0; i < values.length; i ++) {
                if (values[i] == null)
                    continue;
                int bin = Math.min (nbins - 1, (int) (( (Double) values[i] + binrange / 2 - min) / binrange));
                if (confid != null) {
                    bd.write(items[i] + " " + DEFAULT_SHAPE + " " + USE_SIZE[i] + " " + bhex[bin] + " " + DEFAULT_FILL + " " + DEFAULT_POS + " " + (confid[i] != null ? String.format("%.3f±%.3f", values[i], confid[i]) : String.format("%.3f", values[i]))); bd.newLine();
                } else {
                    bd.write(items[i] + " " + DEFAULT_SHAPE + " " + DEFAULT_SIZE + " " + bhex[bin] + " " + DEFAULT_FILL + " " + DEFAULT_POS + " " + String.format("%.3f", values[i])); bd.newLine();
                }
            }
        } else { // values are NOT doubles
            Map<Object, String> possible = new HashMap<>();
            Set<Object> valueset = new HashSet<>();

            for (int i = 0; i < values.length; i ++) {
                if (values[i] == null)
                    continue;
                valueset.add(values[i]);
            }
            int j = 0;
            for (Object value : valueset)
                possible.put(value, NONWHITE_COLORS[j ++]);
            bd.write("LEGEND_TITLE " + dataset_label); bd.newLine();
            StringBuffer shapes = new StringBuffer("LEGEND_SHAPES ");
            StringBuffer colors = new StringBuffer("LEGEND_COLORS ");
            StringBuffer labels = new StringBuffer("LEGEND_LABELS ");
            for (Map.Entry<Object, String> entry : possible.entrySet()) {
                shapes.append(DEFAULT_SHAPE + " ");
                colors.append(entry.getValue() + " ");
                labels.append(String.format("%s ", entry.getKey()));
            }
            bd.write(shapes.toString()); bd.newLine();
            bd.write(colors.toString()); bd.newLine();
            bd.write(labels.toString()); bd.newLine();
            bd.write("MAXIMUM_SIZE 10"); bd.newLine();
            bd.write("DATA"); bd.newLine();
            for (int i = 0; i < values.length; i ++) {
                if (values[i] == null)
                    continue;
                bd.write(items[i] + " " + DEFAULT_SHAPE + " " + DEFAULT_SIZE + " " + possible.get(values[i]) + " " + DEFAULT_FILL + " " + DEFAULT_POS); bd.newLine();
            }
        }
        bd.close();
    }

    public static void save2iTOL(String filename, Object[] items, Object[] values, String dataset_label, int nbins) throws IOException {
        save2iTOL(filename, items, values, dataset_label, nbins, null, null);
    }

    public static void save2iTOL(String filename, Object[] items, Object[] values, String dataset_label) throws IOException {
        save2iTOL(filename, items, values, dataset_label, 10, null, null);
    }

    public static void main(String[] args) {
        Object[] items = new String[] {"WT25", "ASR55", "ASR01", "ASR05", "ASR07"};
        Double[] values = new Double[] {25., 60., 46., 41., 38.};
        try {
            save2iTOL("/Users/mikael/simhome/ASR/ReconMode/test_itol.txt", items, values, "Test_iTOL", 7, 25., 60.);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static class Filter {

    }
}
