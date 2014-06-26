package bn;

import java.io.*;
import java.util.ArrayList;

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
    public PseudoMatrix(ArrayList<Double> values){
        this.matrix = new ArrayList(values);
    }

    public ArrayList<ArrayList<Double>> getValues(){
        return this.matrix;
    }

    public Double getValue(int index1, int index2){
        return this.matrix.get(index1).get(index2);
//        if (index2<index1) {
//            return this.matrix.get(index1).get(index2);
//        }else {
//            return this.matrix.get(index2).get(index1);
//        }
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

    public void addValues(ArrayList<Double> values) {
        this.matrix.add(values);
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
            this.matrix = new ArrayList<>();
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
