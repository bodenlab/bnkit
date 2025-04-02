package bn.example;

import bn.BNet;
import bn.alg.MAP;
import bn.file.BNBuf;
import bn.node.CPT;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.file.TSVFile;

import java.io.IOException;
import java.util.*;

public class TF2Gene {

    /**
     * Loads a TSV file, extracts data for specified variables, and converts it to a structured format.
     *
     * @param filename the path to the TSV file to be loaded
     * @param vars the list of variables to be extracted from the TSV file
     * @return a TSVFile object containing data for the specified variables
     * @throws IOException if an error occurs during file reading
     */
    public static TSVFile load(String filename, List<EnumVariable> vars) throws IOException {
        Object[][] values = TSVFile.loadObjects(filename);
        Object[][] ret = new Object[values.length - 1][vars.size()];
        Map<String, Integer> var2idx = new HashMap<>();
        for (int i = 0; i < vars.size(); i ++) {
            String name = vars.get(i).getName(); // a variable we want to have values for
            for (int j = 1; j < values[0].length; j ++) {
                if (name.equals(values[0][j])) { // we found it in the data file
                    var2idx.put(name, j);
                    break;
                }
            }
        }
        int valid_vars_cnt = 0;
        List<String> valid_vars = new ArrayList<>();
        for (int i = 0; i < vars.size(); i ++) {            // check each variable...
            String name = vars.get(i).getName();
            Integer index = var2idx.get(name);              // we should find the variable in big-file column "index"
            if (index != null) {
                for (int row = 1; row < values.length; row ++)
                    ret[row - 1][valid_vars_cnt] = values[row][index].toString(); // note that we convert everything to text, including Integer "0"
                valid_vars.add(name);
                valid_vars_cnt ++;
            } else
                System.out.println("Variable " + name + " not found in file.");
        }
        String[] valid_vars_arr = new String[valid_vars.size()];
        for (int i = 0; i < valid_vars.size(); i ++)
            valid_vars_arr[i] = valid_vars.get(i);
        return new TSVFile(valid_vars_arr, ret);
    }


    public static void main(String[] args) {


        if (true) {
            Object[] expvals = new Object[]{"0", "low", "medium", "high"};
            Enumerable expvals_enum = new Enumerable(expvals);

            String TFGENEDATA = "TF_target_gene.tsv";
            String GENEDATA = "training_target_discrete.tsv";
            String TFDATA = "training_TF_discrete.tsv";

            String GENEDATA_TST = "testing_target_discrete.tsv";
            String TFDATA_TST = "testing_TF_discrete.tsv";

            TSVFile tf2gene = null;
            try {
                tf2gene = new TSVFile(TSVFile.loadObjects(TFGENEDATA), true);
                System.out.println(tf2gene.getRows().length + " rows loaded.");
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(1);
            }

            Map<String, Set<String>> gene2tf = new HashMap<>();
            Map<String, EnumVariable> tf2var = new HashMap<>();
            Map<EnumVariable, Set<EnumVariable>> gene2tfvars = new HashMap<>();

            BNet bn = new BNet();
            int i = 0;
            String prev = "";
            for (Object[] row : tf2gene.getRows()) {
                String tf = (String) row[0];
                if (!prev.equals(tf)) {
                    EnumVariable tfvar = new EnumVariable(expvals_enum, tf);
                    tf2var.put(tf, tfvar);
                    bn.add(new CPT(tfvar));
                    prev = tf;
                }
                String gene = (String) row[1];
                if (!gene2tf.containsKey(gene)) {
                    gene2tf.put(gene, new java.util.HashSet<>());
                }
                gene2tf.get(gene).add(tf);
                i++;
            }
            System.out.println(i + " rows processed. " + gene2tf.size() + " genes found.");

            Map<Integer, Integer> numbers = new HashMap<>();
            for (String gene : gene2tf.keySet()) {
                Set<String> tfs = gene2tf.get(gene);
                int n = tfs.size();
                if (!numbers.containsKey(n)) {
                    numbers.put(n, 0);
                }
                numbers.put(n, numbers.get(n) + 1);
            }
            for (int n : numbers.keySet()) {
                System.out.println(numbers.get(n) + " genes have " + n + " TFs.");
            }

            Map<String, EnumVariable> gene2var = new HashMap<>();
            for (String gene : gene2tf.keySet()) {
                Set<String> tfs = gene2tf.get(gene);
                int n = tfs.size();
                if (n > 0) {
                    EnumVariable genevar = new EnumVariable(expvals_enum, gene);
                    gene2var.put(gene, genevar);
                }
            }

            // load training data, which in turn will identify what genes that need not be represented in the BN
            List<EnumVariable> tfvars = new ArrayList<>(tf2var.values());
            String[] tfnames = new String[tfvars.size()];
            for (int j = 0; j < tfvars.size(); j++)
                tfnames[j] = tfvars.get(j).getName();
            List<EnumVariable> genevars = new ArrayList<>(gene2var.values());
            String[] genenames = new String[genevars.size()];
            for (int j = 0; j < genevars.size(); j++)
                genenames[j] = genevars.get(j).getName();
            TSVFile tfdata = null;
            TSVFile genedata = null;
            try {
                tfdata = load(TFDATA, tfvars);
                System.out.println(tfdata.getRows().length + " rows included.");
                tfdata.save("training_TF_amended.tsv");
                genedata = load(GENEDATA, genevars);
                System.out.println(genedata.getRows().length + " rows included.");
                genedata.save("training_target_amended.tsv");
            } catch (IOException e) {
                System.err.println("Error loading or saving TF data: " + e.getMessage());
            }

            // build the rest of the BN, namely the TF->Gene relations, including only genes that were found in training data
            int DOWNSAMPLE = 2; // accept this many TFs per gene, else downsample to this number
            for (String gene : gene2tf.keySet()) {
                if (genedata.getColumn(gene) >= 0) { // gene is in the data
                    EnumVariable genevar = gene2var.get(gene);
                    Set<String> tfs = gene2tf.get(gene);
                    Set<String> tfs_in_data = new HashSet<>();
                    for (String tf : tfs) {
                        if (tfdata.getColumn(tf) >= 0)
                            tfs_in_data.add(tf);
                    }
                    int n = tfs_in_data.size();
                    if (n > 0) {
                        Set<String> valid_tfs = tfs_in_data;
                        if (n > DOWNSAMPLE) { // too many to handle, so down-sample
                            List<String> tfs_list = new ArrayList<>(tfs_in_data);
                            Collections.shuffle(tfs_list);
                            valid_tfs = new HashSet<>(tfs_list.subList(0, DOWNSAMPLE));
                        }
                        EnumVariable[] tfarr = new EnumVariable[valid_tfs.size()];
                        int t = 0;
                        for (String tf : valid_tfs) {
                            tfarr[t++] = tf2var.get(tf);
                        }
                        CPT cpt = new CPT(genevar, tfarr);
                        //cpt.setPseudoCounts(0.1);
                        bn.add(cpt);
                    }
                }
            }
            System.out.println(bn.getNodes().size() + " nodes added to BN.");

            BNBuf.save(bn, "tf2gene.bnet");
            TSVFile alltrn = TSVFile.concat(tfdata, genedata);
        /*
        try {
            alltrn.save("training_all_amended.tsv");
        } catch (IOException e) {
            System.err.println("Error saving training data: " + e.getMessage());
        }
         */

            Variable[] allvars = new Variable[alltrn.getHeaders().length];
            Variable[] genevars_arr = new Variable[genevars.size()];
            int gcnt = 0;
            for (int j = 0; j < alltrn.getHeaders().length; j++) {
                String name = alltrn.getHeaders()[j];
                if (tf2var.containsKey(name))
                    allvars[j] = tf2var.get(name);
                else if (gene2var.containsKey(name)) {
                    allvars[j] = gene2var.get(name);
                    genevars_arr[gcnt++] = allvars[j];
                } else {
                    System.err.println("Variable " + name + " not found in BN.");
                    System.exit(2);
                }
            }

            MAP map = new MAP(bn);
            map.train(alltrn.getRows(), allvars, System.currentTimeMillis());
            BNBuf.save(bn, "tf2gene_trained.bnet");

            try {
                //TSVFile tfdata_tst = load(TFDATA_TST, tfvars);
                TSVFile genedata_tst = load(GENEDATA_TST, genevars);
                //System.out.println("Test: " + tfdata_tst.getRows().length + " rows included.");
                System.out.println("Test: " + genedata_tst.getRows().length + " rows included.");
                //map.test(genedata_tst.getRows(), genevars_arr);
                map.test(genedata_tst.getRows(), genevars_arr, new EnumVariable[] {tfvars.get(0)});
            } catch (IOException e) {
                System.err.println("Error loading TF test data: " + e.getMessage());
            }
        }
    }



}
