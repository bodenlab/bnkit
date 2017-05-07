package bn.node;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * This is just a very simple class for loading/storing a tab delimited matrix file and
 * accessing values within this matrix. It primary purpose is to store pseudocounts
 * for CPTPseudo. It is likely that this will be removed at a later and replaced
 * by a more universal and useful class/s for applying prior probabilities.
 *
 * @author julian
 */
public class PseudoMatrix {

    private ArrayList<ArrayList<Double>> matrix;
    public ArrayList<Object> valuemap = new ArrayList<>();

    /**
     * Create a PseudoMatrix
     */
    public PseudoMatrix() {
        this.matrix = new ArrayList<>();
    }

    /**
     * Create a PseudoMatrix
     * @param values
     */
    public PseudoMatrix(ArrayList<ArrayList<Double>> values) {
    	matrix = values;
    }

    public ArrayList<ArrayList<Double>> getValues(){
        return this.matrix;
    }

    public Double getValue(int index1, int index2){
        return this.matrix.get(index1).get(index2);
    }

    public Double getValue(int index1, Object key){
        if (valuemap == null){
            return null;
        }
        Double value = null;
        for (int i = 0; i < valuemap.size(); i++){
            if ((valuemap.get(i).toString()).equals(key.toString())){
                value = getValue(index1, i);
            }
        }
        return value;
    }

    public Double getValue(Object key){
        if (valuemap == null){
            return null;
        }
        Double value = null;
        for (int i = 0; i < valuemap.size(); i++){
            if (valuemap.get(i).equals(key)){
                value = getValue(0, i);
            }
        }
        return value;
    }

    public void normalize(){
        for (int row = 0; row < matrix.size(); row ++){
            double total = Sum(matrix.get(row));
            ArrayList<Double> newrow = new ArrayList<Double>();
            if (total == 0){
                return;
            } else {
                for (int col = 0; col < matrix.get(row).size(); col++) {
                    newrow.add(matrix.get(row).get(col)/total);
                }
            }
            matrix.set(row, newrow);
        }
    }

    public double Sum(List<Double> values){
        double sum = 0;
        for (int i = 0; i < values.size(); i++) {
            sum += values.get(i);
        }
        return sum;
    }

    public void scaleValues(double scalar){
        ArrayList<ArrayList<Double>> newmatrix = new ArrayList<ArrayList<Double>>();
        for (int i =0; i < this.matrix.size(); i++){
            ArrayList<Double> row = new ArrayList<Double>();
            for (int j =0; j < this.matrix.get(i).size(); j++){
                Double curval = this.matrix.get(i).get(j);
                row.add(curval * scalar);
            }
            newmatrix.add(row);
        }
        this.matrix = newmatrix;
    }

    public void scaleDiagonal(double scalar){
        ArrayList<ArrayList<Double>> newmatrix = new ArrayList<ArrayList<Double>>();
        for (int i =0; i < this.matrix.size(); i++){
            ArrayList<Double> row = new ArrayList<Double>();
            for (int j =0; j < this.matrix.get(i).size(); j++){
                Double curval = this.matrix.get(i).get(j);
                if (i == j){
                    row.add(curval * scalar);
                } else{
                    row.add(curval);
                }
            }
            newmatrix.add(row);
        }
        this.matrix = newmatrix;
    }

    public void addValues(ArrayList<Double> values) {
        this.matrix.add(values);
    }

    public void setMapping(List<Object> amap){
        valuemap.clear();
        for (int i = 0; i < amap.size(); i++)
            valuemap.add(amap.get(i));
        if (valuemap.size() != this.matrix.get(0).size()){
            System.err.println("Matrix map does not match the length of matrix entries");
        }
    }

    public void addValues(ArrayList<Double> values, int repeat) {
        for (int i = 0; i < repeat; i++)
            this.matrix.add(values);
    }

    /**
     *
     * @return
     */
    public ArrayList<Double> getMarginalized(){
        ArrayList<Double> out = new ArrayList<>();
        for (int i = 0; i < this.matrix.size(); i++){
            double sum = 0;
            for (double j : this.matrix.get(i))
                sum += j;
            out.add(sum);
        }
        return out;
    }

    public void load(String filename){
        File file = new File(filename);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line;
            ArrayList<ArrayList<Double>> lines = new ArrayList<>();
//            this.matrix = new ArrayList<>();
            try {

                while ((line = reader.readLine()) != null) {
                    String[] entries = line.split("\t");
                    ArrayList<Double> values = new ArrayList<>();
                    for (int i = 0; i < entries.length; i++){ // for each entry on line
                        values.add(Double.parseDouble(entries[i]));
                    }
                    this.matrix.add(values);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (FileNotFoundException e){
            e.printStackTrace();
        }
    }

    /**
     * Load a tab delimited matrix from a file
     * @param filename -> filename
     * @param repeat -> how many times ou want to repeat add each line
     *               (useful for single vector matrices that need to be applied for multiple parents in a BN)
     */
    public void load(String filename, int repeat){
        File file = new File(filename);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line;
            ArrayList<ArrayList<Double>> lines = new ArrayList<>();
            this.matrix = new ArrayList<>();
            try {

                while ((line = reader.readLine()) != null) {
                    String[] entries = line.split("\t");
                    ArrayList<Double> values = new ArrayList<>();
                    for (int i = 0; i < entries.length; i++){ // for each entry on line
                        values.add(Double.parseDouble(entries[i]));
                    }
                    int cnt = 0;
                    while (cnt < repeat){
                        this.matrix.add(values);
                        cnt++;
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (FileNotFoundException e){
            e.printStackTrace();
        }
    }


}
