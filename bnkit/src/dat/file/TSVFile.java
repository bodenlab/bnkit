package dat.file;

import java.io.*;
import java.util.*;

/**
 * Created by mikael on 2/08/2016. Copied from binfkit March 2019.
 */
public class TSVFile {

    List<Object[]> rows = new ArrayList<>();
    int ncols = 0;
    Map<String, Integer> headers = null;
    Map<Integer, Map<Object, int[]>> indexMap = new HashMap<>();

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
                        headers.put((String) objects[0][i], i);
                    } catch (ClassCastException e) {
                        throw new RuntimeException("Failed to index header: not a text string \"" + objects[0][i] + "\"");
                    }
                }
            }
            int index = 0;
            for (int i = (useHeader ? 1 : 0); i < objects.length; i++) { // start at 1 and skip header conditionally
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
                rows.add(index, objects[i]);
            }
        }
    }

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
     * @return the column index (starting at 0)
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
                            values[i] = tokens[i]; // value is a string
                        }
                    }
                }
                alldata.add(values);
            }
            line = br.readLine();
        }
        data = new Object[alldata.size()][];
        for (int k = 0; k < data.length; k++) {
            data[k] = new Object[alldata.get(k).length];
            for (int j = 0; j < data[k].length; j++) {
                data[k][j] = alldata.get(k)[j];
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
                Object val = row[i];
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

}
