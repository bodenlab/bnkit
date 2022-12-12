package bn.example;

import bn.BNet;
import bn.Distrib;
import bn.Predef;
import bn.alg.*;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Variable;
import dat.file.TSVFile;

import java.io.IOException;
import java.util.*;

public class TFCascade {

    public static void usage() {
        usage(null);
    }
    public static void usage(String error) {
        if (error != null) {
            System.err.println(error);
        }
        System.out.println("Usage: TFCascade <command> [ <options> | -f <data-file> | -b <bnet-file> | -p <pred-file> ]");
        System.out.println("Command is one of");
        System.out.println("\tdefine\trequires <data-file> to define variables, and saves untrained BN to <bnet-file>");
        System.out.println("\ttrain\trequires <data-file> to train parameters, and loads an untrained from and saves trained BN to <bnet-file>");
        System.out.println("\tpredict\trequires <data-file> to make predictions of <query-var>, which are output on standard output or as a <pred-file> (for the latter, each row is a distribution over <query-var>)");
        System.out.println("Options include");
        System.out.println("\t-g <num, default 4> (generation 3 or 4)");
        System.out.println("\t-mc <num, default 3> (number of categorical values that Celltype latent can take; effective only for define-command; 1 means Boolean, 0 means no celltype latent)");
        System.out.println("\t-mt1 <num, default 3> (number of categorical values that TF1 latent can take; effective only for define-command; 1 means Boolean, 0 means no TF1 latent)");
        System.out.println("\t-r <rounds, default 100> (max number of EM rounds; effective only for train-command)");
        System.out.println("\t-nudge <count> (add negative count to pre-tuned BN; effective only for train-command)");
        System.out.println("\t-exit (exit after preset of counts, but before training; effective only for train-command)");
        System.out.println("\t-v <query-var, default Bnd> (must be either CObs, CLat, CAct, TF1, TF2 or Bnd [default]; effective only for predict-command)");
        System.out.println("\t-even (use a uniform prior; effective only for predict-command)");
        System.out.println("\t-tabulate (use naive tabulation of data; effective only for predict-command)");
        System.out.println("<data-file> is a tab-separated file in the order celltype, TF source, TF target and [optionally] binding status (otherwise assumed true)");
        System.out.println("<pred-file> is a tab-separated file with a header identifying values defining the distribution");
        System.out.println("<bnet-file> is a bnkit XML file describing the BN");
        System.out.println("This is generation 4 of TFCascade");
    }

    public static void main(String[] args) {

        String command  = null;
        String datafile = null;
        String predfile = null;
        String bnetfile = null;
        String queryvar = null;
        String[] vars_in_file = new String[] {"CAct", "TF1", "TF2", "Bnd"};
        int ROUNDS = 100;
        int MIXTURE_CT = 3;
        int MIXTURE_T1 = 3;
        int GEN = 4;
        Double NUDGE = null;
        boolean EVENPRIOR = false;
        boolean PRESET_EXIT = false;

        for (int a = 0; a < args.length; a ++) {
            if (a == 0)
                command = args[a];
            else {
                if (args[a].startsWith("-f") && a < args.length - 1)
                    datafile = args[++a];
                else if (args[a].startsWith("-b") && a < args.length - 1)
                    bnetfile = args[++a];
                else if (args[a].startsWith("-p") && a < args.length - 1)
                    predfile = args[++a];
                else if (args[a].startsWith("-v") && a < args.length - 1)
                    queryvar = args[++a];
                else if (args[a].startsWith("-g") && a < args.length - 1)
                    GEN = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-r") && a < args.length - 1)
                    ROUNDS = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-nudge") && a < args.length - 1)
                    NUDGE = Double.parseDouble(args[++a]);
                else if (args[a].startsWith("-even"))
                    EVENPRIOR = true;
                else if (args[a].startsWith("-exit"))
                    PRESET_EXIT = true;
                else if (args[a].startsWith("-mc") && a < args.length - 1)
                    MIXTURE_CT = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-mt") && a < args.length - 1)
                    MIXTURE_T1 = Integer.parseInt(args[++a]);
            }
        }
        if (command != null)
            System.out.println("RUN: " + command);
        else {
            usage("Command not given");
            System.exit(1);
        }

        if (datafile != null)
            System.out.println("DATA FILE: " + datafile);
        if (predfile != null)
            System.out.println("PRED FILE: " + predfile);
        if (bnetfile != null)
            System.out.println("BNET FILE: " + bnetfile);
        if (queryvar != null)
            System.out.println("QUERY VAR: " + queryvar);

        if (command.equalsIgnoreCase("define") && datafile != null) {
            TSVFile tsv = null;
            String[] tf_names = null;
            String[] ct_names = null;
            try {
                tsv = new TSVFile(datafile, false);
                Set<Object> tfs2 = tsv.getValues(2);
                Set<Object> tfs1 = tsv.getValues(1);
                Set<Object> ctypes = tsv.getValues(0);
                Set<String> tfs = new HashSet<>();
                for (Object obj : tfs1)
                    tfs.add((String) obj);
                for (Object obj : tfs2)
                    tfs.add((String) obj);
                tf_names = new String[tfs.size()];
                ct_names = new String[ctypes.size()];
                int i = 0;
                for (String obj : tfs) {
                    tf_names[i] = obj;
                    i ++;
                }
                i = 0;
                for (Object obj : ctypes) {
                    ct_names[i] = obj.toString();
                    i ++;
                }
                Arrays.sort(tf_names);
                Arrays.sort(ct_names);
            } catch (IOException e) {
                System.err.println(e);
                System.exit(1);
            }

            System.out.println("Defining variables and tables");
            // Define variables
            EnumVariable T1 = Predef.Nominal(tf_names, "TF1");
            System.out.println("Defining " + T1 + " with " + T1.getDomain().size() + " values");
            EnumVariable T2 = Predef.Nominal(tf_names, "TF2");
            System.out.println("Defining " + T2 + " with " + T2.getDomain().size() + " values");
            EnumVariable CObs = Predef.Nominal(ct_names, "CObs");
            System.out.println("Defining " + CObs + " with " + CObs.getDomain().size() + " values");
            EnumVariable CAct = Predef.Nominal(ct_names, "CAct");
            System.out.println("Defining " + CAct + " with " + CAct.getDomain().size() + " values");
            EnumVariable CLat = null;
            if (MIXTURE_CT > 1)
                CLat = Predef.Number(MIXTURE_CT, "CLat");
            else if (MIXTURE_CT == 1)
                CLat = Predef.Boolean("CLat");
            else
                CLat = null;
            if (CLat != null)
                System.out.println("Defining " + CLat + " with " + CLat.getDomain().size() + " values");
            EnumVariable T1Lat = null;
            if (MIXTURE_T1 > 1)
                T1Lat = Predef.Number(MIXTURE_T1, "T1Lat");
            else if (MIXTURE_T1 == 1)
                T1Lat = Predef.Boolean("T1Lat");
            else
                T1Lat = null;
            if (T1Lat != null)
                System.out.println("Defining " + T1Lat + " with " + T1Lat.getDomain().size() + " values");
            EnumVariable Bnd = Predef.Boolean("Bnd");

            // Define nodes (connecting the variables into an acyclic graph, i.e. the structure)
            BNet bn = new BNet();
            CPT bnd = new CPT(Bnd);
            CPT t1 = new CPT(T1);
            CPT t2 = null;
            CPT t1lat = null;
            CPT clat = null;
            CPT cact = null;
            if (GEN == 4) {
                if (T1Lat == null) {
                    t2 = new CPT(T2, T1, Bnd);
                    if (CLat != null) {
                        clat = new CPT(CLat, T2, Bnd);
                        cact = new CPT(CAct, CLat);
                        bn.add(t1, t2, bnd, cact, clat);
                    } else {
                        cact = new CPT(CAct, T2, Bnd);
                        bn.add(t1, t2, bnd, cact);
                    }
                } else {
                    t1lat = new CPT(T1Lat, T1);
                    t2 = new CPT(T2, T1Lat, Bnd);
                    if (CLat != null) {
                        clat = new CPT(CLat, T2, Bnd);
                        cact = new CPT(CAct, CLat);
                        bn.add(t1, t2, bnd, cact, clat, t1lat);
                    } else {
                        cact = new CPT(CAct, T2, Bnd);
                        bn.add(t1, t2, bnd, cact, t1lat);
                    }
                }
            } else if (GEN == 3) {
                if (T1Lat == null) {
                    t2 = new CPT(T2, T1, Bnd);
                    cact = new CPT(CAct);
                    if (CLat != null) {
                        clat = new CPT(CLat, CAct);
                        t2 = new CPT(T2, CLat, T1, Bnd);
                        bn.add(t1, t2, bnd, cact, clat);
                    } else {
                        t2 = new CPT(T2, CAct, T1, Bnd);
                        bn.add(t1, t2, bnd, cact);
                    }
                } else {
                    t1lat = new CPT(T1Lat, T1);
                    t2 = new CPT(T2, T1Lat, Bnd);
                    cact = new CPT(CAct);
                    if (CLat != null) {
                        clat = new CPT(CLat, CAct);
                        t2 = new CPT(T2, CLat, T1, Bnd);
                        bn.add(t1, t2, bnd, cact, clat);
                    } else {
                        t2 = new CPT(T2, CAct, T1, Bnd);
                        bn.add(t1, t2, bnd, cact);
                    }
                }
            } else
                System.err.println("Gen " + GEN + " is not recognised");
            System.out.println("This is TFCascade gen " + GEN);
            if (bnetfile != null)
                BNBuf.save(bn, bnetfile);
        } else

        if ((command.equalsIgnoreCase("train") || command.equalsIgnoreCase("predict")) && datafile != null) {

            BNet bn = BNBuf.load(bnetfile);
            System.out.println("Loaded variables and tables");
            // Nodes, and their variables
            CPT t1 = (CPT)bn.getNode("TF1");
            CPT t2 = (CPT)bn.getNode("TF2");
            CPT cact = (CPT)bn.getNode("CAct");
            CPT clat = (CPT)bn.getNode("CLat");
            CPT bnd = (CPT)bn.getNode("Bnd");

            EnumVariable T1 = (EnumVariable)t1.getVariable();
            System.out.println("Loaded " + T1 + " with " + T1.getDomain().size() + " values");
            EnumVariable T2 = (EnumVariable)t2.getVariable();;
            System.out.println("Loaded " + T2 + " with " + T2.getDomain().size() + " values");
            EnumVariable CAct = (EnumVariable)cact.getVariable();
            System.out.println("Loaded " + CAct + " with " + CAct.getDomain().size() + " values");
            EnumVariable CLat = null;
            if (clat != null) {
                CLat = (EnumVariable) clat.getVariable();
                System.out.println("Loaded " + CLat + " with " + CLat.getDomain().size() + " values");
            }
            EnumVariable Bnd = (EnumVariable)bnd.getVariable();

            Object[] cts = CAct.getDomain().getValues();
            Map<Object, Map<Object, Set<Object>>> posmap = new HashMap<>();
            Map<Object, Map<Object, Set<Object>>> negmap = new HashMap<>();
            for (Object ct : cts)
                posmap.put(ct, new HashMap<>());
            for (Object ct : cts)
                negmap.put(ct, new HashMap<>());
            Object[] tfs = t1.getVariable().getDomain().getValues();
            if (command.equalsIgnoreCase("train")) {
                Object[][] obs_preset = DataBuf.load(datafile, new Variable[]{CAct, T1, T2}, false);
                System.out.println(obs_preset.length + " samples read");
                // Stage 2: parameterise the latent variable, and class variables
                for (int i = 0; i < obs_preset.length; i++) {
                    Map<Object, Set<Object>> postfs = posmap.get(obs_preset[i][0]);
                    Map<Object, Set<Object>> negtfs = negmap.get(obs_preset[i][0]);
                    if (!postfs.containsKey(obs_preset[i][1])) {
                        Set<Object> tf2s = new HashSet<>();
                        postfs.put(obs_preset[i][1], tf2s);
                    }
                    postfs.get(obs_preset[i][1]).add(obs_preset[i][2]);
                    if (!negtfs.containsKey(obs_preset[i][1])) // eventually there will be negatives for this tf1
                        negtfs.put(obs_preset[i][1], new HashSet<>());
                }
                // Check how many negatives
                int negcnt = 0;
                Object[] tf_names = T2.getDomain().getValues();
                for (Map.Entry<Object, Map<Object, Set<Object>>> ct_tf1 : posmap.entrySet()) {
                    Object ct = ct_tf1.getKey();
                    for (Map.Entry<Object, Set<Object>> tf1_tf2s : ct_tf1.getValue().entrySet()) {
                        Object tf1 = tf1_tf2s.getKey();
                        for (Object tf2 : tf_names) {
                            if (!tf1_tf2s.getValue().contains(tf2)) {
                                negmap.get(ct).get(tf1).add(tf2);
                                negcnt += 1;
                            }
                        }
                    }
                }
                Object[][] obs_train = new Object[obs_preset.length + negcnt][4]; // four variables
                for (int i = 0; i < obs_preset.length; i++) {
                    obs_train[i][0] = obs_preset[i][0]; // CAct
                    obs_train[i][1] = obs_preset[i][1]; // T1
                    obs_train[i][2] = obs_preset[i][2]; // T2
                    obs_train[i][3] = true; // Bnd
                }
                int cnt = 0;
                for (Map.Entry<Object, Map<Object, Set<Object>>> ct_tf1 : negmap.entrySet()) {
                    Object ct = ct_tf1.getKey();
                    for (Map.Entry<Object, Set<Object>> tf1_tf2s : ct_tf1.getValue().entrySet()) {
                        Object tf1 = tf1_tf2s.getKey();
                        for (Object tf2 : tf1_tf2s.getValue()) {
                            obs_train[obs_preset.length + cnt][0] = ct;
                            obs_train[obs_preset.length + cnt][1] = tf1;
                            obs_train[obs_preset.length + cnt][2] = tf2;
                            obs_train[obs_preset.length + cnt][3] = false;
                            cnt += 1;
                        }
                    }
                }
                System.out.println(obs_train.length + " samples total (with negatives)");
                LearningAlg em = new EM(bn);
                ((EM) em).setMaxRounds(ROUNDS);
                ((EM) em).setEMOption(2);
                em.train(obs_train, new Variable[]{CAct, T1, T2, Bnd}, System.currentTimeMillis());
                BNBuf.save(bn, bnetfile + ".trained");
            } else if (command.equalsIgnoreCase("predict")) { // predict
                if (EVENPRIOR) {
                    t1.put(EnumDistrib.uniform(T1.getDomain()));
                    t2.put(EnumDistrib.uniform(T2.getDomain()));
                    cact.put(EnumDistrib.uniform(CAct.getDomain()));
                }
                TSVFile input = null;
                try {
                    input = new TSVFile(datafile, false);
                    if (queryvar == null)
                        queryvar = "Bnd";
                    Object[][] rows = input.getRows();
                    Object[][] output = new Object[rows.length + 1][];
                    for (int k = 0; k < rows.length; k ++) {
                        VarElim ve = new VarElim();
                        ve.instantiate(bn);
                        EnumVariable qvar = null;
                        Integer overwrite = null;
                        Object overwrite_value = null;
                        if (queryvar.equalsIgnoreCase("TF1")) {
                            qvar = T1;
                            overwrite = 1;
                            overwrite_value = rows[k][overwrite.intValue()].toString();
                        } else if (queryvar.equalsIgnoreCase("TF2")) {
                            qvar = T2;
                            overwrite = 2;
                            overwrite_value = rows[k][overwrite.intValue()].toString();
                        } else if (queryvar.equalsIgnoreCase("CAct")) {
                            qvar = CAct;
                            overwrite = 0;
                            overwrite_value = rows[k][overwrite.intValue()].toString();
                        } else if (queryvar.equalsIgnoreCase("CObs") || queryvar.equalsIgnoreCase("CAct"))
                            qvar = CAct;
                        else if (queryvar.equalsIgnoreCase("CLat"))
                            qvar = CLat;
                        else {
                            qvar = Bnd;
                            if (rows[k].length > 3) {
                                overwrite = 3;
                                overwrite_value = Boolean.parseBoolean(rows[k][overwrite.intValue()].toString());
                            }
                        }
                        if (k == 0) {
                            output[k] = new Object[qvar.getDomain().size() + rows[k].length];
                            for (int kk = 0; kk < vars_in_file.length; kk ++)
                                output[k][kk] = vars_in_file[kk];
                        }
                        output[k + 1] = new Object[qvar.getDomain().size() + rows[k].length];
                        for (int kk = 0; kk < rows[k].length; kk ++) {
                            output[k + 1][kk] = rows[k][kk];
                            if (kk == 0 && !queryvar.equalsIgnoreCase("CAct")) {
                                cact.setInstance(rows[k][kk].toString());
                            } else if (kk == 1 && !queryvar.equalsIgnoreCase("TF1")) {
                                t1.setInstance(rows[k][kk].toString());
                            } else if (kk == 2 && !queryvar.equalsIgnoreCase("TF2")) {
                                t2.setInstance(rows[k][kk].toString());
                            } else if (kk == 3 && !queryvar.equalsIgnoreCase("Bnd")) {
                                bnd.setInstance(Boolean.parseBoolean(rows[k][kk].toString()));
                            }
                        }
                        if (rows[k].length <= 3 && !queryvar.equalsIgnoreCase("Bnd")) {
                            bnd.setInstance(true);
                        }
                        Query q = ve.makeQuery(qvar);
                        CGTable r = (CGTable) ve.infer(q);
                        Distrib d = r.query(qvar);
                        for (int i = 0; i < qvar.getDomain().size(); i++) {
                            Object instance = qvar.getDomain().get(i);
                            if (k == 0)
                                output[k][i + rows[k].length] = instance;
                            Double p = d.get(instance);
                            output[k + 1][i + rows[k].length] = p;
                            if (overwrite != null)
                                if (overwrite_value.toString().equals(instance.toString()))
                                    output[k + 1][overwrite.intValue()] = p;
                        }
                    }
                    if (predfile != null)
                        TSVFile.saveObjects(predfile, output);
                    else {
                        for (int i = 0; i < output.length; i ++) {
                            for (int j = 0; j < Math.min(output[i].length, rows[0].length); j ++) {
                                try {
                                    System.out.printf("%8.6f\t", output[i][j]);
                                } catch (java.util.IllegalFormatConversionException e) {
                                    System.out.printf("%-8s\t", output[i][j]);
                                }
                            }
                            System.out.println();
                        }
                    }
                } catch (IOException e) {
                    usage(e.getMessage());
                    System.exit(2);
                }
            }
        }
        System.out.println("Done.");
    }
}
