package bn.file;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.node.CPT;
import dat.EnumVariable;
import dat.Variable;
import dat.file.TSVFile;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 * Class to extract a table from a BN, and save it to a TSV file.
 */
public class BN2TSV {


    public static void usage() {
        usage(null);
    }
    public static void usage(String err) {
        if (err != null)
            System.err.println(err);
        System.err.println("BN2TSV [-f <bnetfile> -v <variable> -tsv <output-tsvfile> -transpose -header]");
        System.err.println("\t<bnetfile> is the XML file with the BN");
        System.err.println("\t<output-tsvfile> is the TSV file that is generated");
        System.err.println("\t<variable> is the name of the variable for which a table will be generated");
        System.err.println("\t-transpose means that table is \"turned on its side\" and therefore rows are indexed by the values the node assigns.");
        System.err.println("\t-header means that the table is decorated with headers");
        System.err.println("\t-keepnull means that null rows or columns will not be removed (like they will by default)");
    }

    public static void main(String[] args) {

        String bnetfile = null;
        String tsvfile = null;
        String queryvar = null;
        boolean TRANSPOSE = false;
        boolean HEADERS = false;
        boolean KEEPNULL = false;

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-f") && a < args.length - 1)
                bnetfile = args[++a];
            else if (args[a].startsWith("-tsv") && a < args.length - 1)
                tsvfile = args[++a];
            else if (args[a].startsWith("-v") && a < args.length - 1)
                queryvar = args[++a];
            else if (args[a].startsWith("-transp"))
                TRANSPOSE = true;
            else if (args[a].startsWith("-head"))
                HEADERS = true;
            else if (args[a].startsWith("-keepn"))
                KEEPNULL = true;
        }
        if (bnetfile == null || queryvar == null || tsvfile == null) {
            usage();
            System.exit(0);
        }
        BNet bn = BNBuf.load(bnetfile);
        if (bn == null) {
            usage("Bayesian network could not be loaded from " + bnetfile);
            System.exit(1);
        }
        BNode node = bn.getNode(queryvar);
        if (node == null) {
            usage("Bayesian network node " + queryvar + " is not available");
            System.err.println("These nodes are found in " + bnetfile);
            for (BNode n : bn.getNodes()) {
                try {
                    CPT cpt = (CPT) n;
                    System.err.println("\t" + n + " [valid]");
                } catch (ClassCastException e) {
                    System.err.println("\t" + n + " [invalid]");
                }
            }
            System.exit(2);
        }
        try {
            CPT cpt = (CPT)node;
            Set<Integer> nulls = new HashSet<>();
            Object[][] rows = null;
            if (!TRANSPOSE) {
                EnumVariable v = cpt.getVariable();
                if (cpt.isPrior()) {
                    rows = new Object[1 + (HEADERS?1:0)][v.size()];
                    int rowcnt = 0;
                    if (HEADERS) {
                        for (int j = 0; j < v.getDomain().getValues().length; j ++)
                            rows[rowcnt][j] = v.getDomain().getValues()[j].toString();
                        rowcnt += 1;
                    }
                    try {
                        Distrib val = cpt.getDistrib();
                        if (val != null) {
                            for (int j = 0; j < v.getDomain().getValues().length; j ++)
                                rows[rowcnt][j] = val.get(v.getDomain().getValues()[j]);
                        } else {
                            for (int j = 0; j < v.getDomain().getValues().length; j++)
                                rows[rowcnt][j] = "null";
                        }
                    } catch (ClassCastException e) {
                        System.err.println("Invalid distribution");
                        System.exit(3);
                    }
                } else {
                    rows = new Object[cpt.getTable().getSize() + (HEADERS ? 1 : 0)][cpt.getParents().size() + v.size()];
                    int rowcnt = 0;
                    if (HEADERS) {
                        for (int j = 0; j < cpt.getParents().size(); j++)
                            rows[rowcnt][j] = cpt.getParents().get(j).getName();
                        for (int j = 0; j < v.getDomain().getValues().length; j++)
                            rows[rowcnt][j + cpt.getParents().size()] = v.getDomain().getValues()[j].toString();
                        rowcnt += 1;
                    }
                    for (int i = 0; i < cpt.getTable().getSize(); i++) {
                        Object[] key = cpt.getTable().getKey(i);
                        try {
                            Distrib val = (Distrib) cpt.getTable().getValue(i);
                            for (int j = 0; j < key.length; j++)
                                rows[rowcnt][j] = key[j].toString();
                            if (val != null) {
                                for (int j = 0; j < v.getDomain().getValues().length; j++)
                                    rows[rowcnt][j + cpt.getParents().size()] = val.get(v.getDomain().getValues()[j]);
                            } else {
                                for (int j = 0; j < v.getDomain().getValues().length; j++) {
                                    rows[rowcnt][j + cpt.getParents().size()] = "null";
                                    nulls.add(rowcnt);
                                }
                            }
                        } catch (ClassCastException e) {
                            System.err.println("Invalid distribution");
                            System.exit(3);
                        }
                        rowcnt += 1;
                    }
                }
                try {
                    if (!KEEPNULL) {
                        Object[][] nrows = new Object[rows.length - nulls.size()][rows[0].length];
                        int r = 0;
                        for (int rr = 0; rr < rows.length; rr++) {
                            if (!nulls.contains(rr)) {
                                for (int cc = 0; cc < rows[rr].length; cc++)
                                    nrows[r][cc] = rows[rr][cc];
                                r += 1;
                            }
                        }
                        TSVFile.saveObjects(tsvfile, nrows);
                    } else
                        TSVFile.saveObjects(tsvfile, rows);
                } catch (IOException e) {
                    usage("Failed to save to TSV-file " + tsvfile);
                    System.exit(4);
                }
            } else {
                EnumVariable v = cpt.getVariable();
                if (cpt.isPrior()) {
                    rows = new Object[v.size()][1 + 1];
                    int rowcnt = 0;
                    Object[] values = v.getDomain().getValues();
                    for (int i = 0; i < values.length; i++) {
                        rows[rowcnt][0] = values[i].toString();
                        try {
                            Distrib val = cpt.getDistrib();
                            if (val != null)
                                rows[rowcnt][1] = val.get(values[i]);
                            else {
                                rows[rowcnt][1] = "null";
                                nulls.add(rowcnt);
                            }
                        } catch (ClassCastException e) {
                            System.err.println("Invalid distribution");
                            System.exit(3);
                        }
                        rowcnt += 1;
                    }
                } else {
                    rows = new Object[v.size() + (HEADERS ? cpt.getTable().getKey(0).length : 0)][cpt.getTable().getSize() + 1];
                    int rowcnt = 0;
                    if (HEADERS) {
                        for (int k = 0; k < cpt.getTable().getSize(); k++) {
                            Object key[] = cpt.getTable().getKey(k);
                            for (int rk = 0; rk < key.length; rk++) {
                                rows[rowcnt + rk][0] = "";
                                rows[rowcnt + rk][k + 1] = key[rk];
                            }
                        }
                        rows[rowcnt + cpt.getTable().getKey(0).length - 1][0] = v.getName();
                        rowcnt += cpt.getTable().getKey(0).length;
                    }
                    Object[] values = v.getDomain().getValues();
                    for (int i = 0; i < values.length; i++) {
                        rows[rowcnt][0] = values[i].toString();
                        for (int k = 0; k < cpt.getTable().getSize(); k++) {
                            try {
                                Distrib val = (Distrib) cpt.getTable().getValue(k);
                                if (val != null)
                                    rows[rowcnt][k + 1] = val.get(values[i]);
                                else {
                                    rows[rowcnt][k + 1] = "null";
                                    nulls.add(k + 1);
                                }
                            } catch (ClassCastException e) {
                                System.err.println("Invalid distribution");
                                System.exit(3);
                            }
                        }
                        rowcnt += 1;
                    }
                }
                try {
                    if (!KEEPNULL) {
                        Object[][] nrows = new Object[rows.length][rows[0].length - nulls.size()];
                        for (int rr = 0; rr < rows.length; rr++) {
                            int c = 0;
                            for (int cc = 0; cc < rows[rr].length; cc++) {
                                if (!nulls.contains(cc))
                                    nrows[rr][c++] = rows[rr][cc];
                            }
                        }
                        TSVFile.saveObjects(tsvfile, nrows);
                    } else
                        TSVFile.saveObjects(tsvfile, rows);
                } catch (IOException e) {
                    usage("Failed to save to TSV-file " + tsvfile);
                    System.exit(4);
                }
            }
        } catch (ClassCastException e) {
            usage("Node for variable " + node + " is not a CPT [invalid]");
        }

    }
}
