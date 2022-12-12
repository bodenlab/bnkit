package bn.example;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.Predef;
import bn.alg.*;
import bn.file.BNBuf;
import bn.file.DataBuf;
import bn.node.CPT;
import bn.node.GDT;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.file.TSVFile;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by mikael on 1/5/18.
 */
public class MergeMixtures {

    public static void main(String[] args) {

        boolean TRAIN = true;

        if (args.length < 1) {
            System.out.println("Usage: MergeMixtures <command>");
            System.out.println("command is ");
            System.out.println("\ttrain");
            System.out.println("\tquery");
            System.out.println("For train the following arguments can follow");
            System.out.println("\t-bx1 <BN-file-Mixture-Y1> -bx2 <BN-file-Mixture-Y2> -mx1 <TSV-file-with-Y1s> -mx2 <TSV-file-with-Y2s> -bm <BN-file-merged>");
            System.out.println("For query the following arguments can follow");
            System.out.println("\t-bm <BN-file> -in <query-TSV-file> -v <query-var> -out <inferred-TSV-file");
            System.exit(1);
        }

        String command  = null;
        String tsvmix1 = null;
        String tsvmix2 = null;
        String inputfile = null;
        String outputfile = null;
        String bnetmix1 = null;
        String bnetmix2 = null;
        String bnetmerge = null;
        String queryvar = null;
        int EM_OPTION = 1;
        int ROUNDS = -1;
        int MIXTURE_CT = 5;
        int MIXTURE_TF = 10;
        int GEN = 1;

        for (int a = 0; a < args.length; a ++) {
            if (a == 0)
                command = args[a];
            else {
                if (args[a].startsWith("-mx1") && a < args.length - 1)
                    tsvmix1 = args[++a];
                else if (args[a].startsWith("-mx2") && a < args.length - 1)
                    tsvmix2 = args[++a];
                else if (args[a].startsWith("-bx1") && a < args.length - 1)
                    bnetmix1 = args[++a];
                else if (args[a].startsWith("-bx2") && a < args.length - 1)
                    bnetmix2 = args[++a];
                else if (args[a].startsWith("-bm") && a < args.length - 1)
                    bnetmerge = args[++a];
                else if (args[a].startsWith("-em") && a < args.length - 1)
                    EM_OPTION = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-in") && a < args.length - 1)
                    inputfile = args[++a];
                else if (args[a].startsWith("-out") && a < args.length - 1)
                    outputfile = args[++a];
                else if (args[a].startsWith("-v") && a < args.length - 1)
                    queryvar = args[++a];
                else if (args[a].startsWith("-g") && a < args.length - 1)
                    GEN = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-r") && a < args.length - 1)
                    ROUNDS = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-mc") && a < args.length - 1)
                    MIXTURE_CT = Integer.parseInt(args[++a]);
                else if (args[a].startsWith("-mt") && a < args.length - 1)
                    MIXTURE_TF = Integer.parseInt(args[++a]);
            }
        }

        BNet bn = null;
        CPT ct = null, ct_lat = null, tf_lat = null, tf = null, mix_hm = null, mix_exp = null;
        EnumVariable CObs = null, TF = null, CT_lat = null, TF_lat = null, Mix_Exp = null, Mix_HM = null;
        Variable Y1 = null, Y2 = null;

        if (bnetmerge != null && bnetmix1 == null && bnetmix2 == null) { // load BN
            bn = BNBuf.load(bnetmerge); // RPKM
            System.out.println("Loaded variables and tables");
            // Nodes, and their variables
            GDT y1 = (GDT) bn.getNode("Exp");
            if (y1 != null) {
                y1.setTrainable(false);
                Y1 = y1.getVariable();
                System.out.println("Loaded " + Y1 + " defined with domain " + Y1.getDomain());
            }
            GDT y2 = (GDT) bn.getNode("H3K27me3");
            if (y2 != null) {
                y2.setTrainable(false);
                Y2 = y2.getVariable();
                System.out.println("Loaded " + Y2 + " defined with domain " + Y2.getDomain());
            }
            ct = (CPT) bn.getNode("CObs");
            if (ct != null) {
                CObs = ct.getVariable();
                System.out.println("Loaded " + CObs + " defined with domain " + CObs.getDomain());
            }
            tf = (CPT) bn.getNode("TF");
            if (tf != null) {
                TF = tf.getVariable();
                System.out.println("Loaded " + TF + " defined with domain " + TF.getDomain());
            }
            ct_lat = (CPT) bn.getNode("CT_lat");
            if (ct_lat != null) {
                CT_lat = ct_lat.getVariable();
                System.out.println("Loaded " + CT_lat + " defined with domain " + CT_lat.getDomain());
            }
            tf_lat = (CPT) bn.getNode("TF_lat");
            if (tf_lat != null) {
                TF_lat = tf_lat.getVariable();
                System.out.println("Loaded " + TF_lat + " defined with domain " + TF_lat.getDomain());
            }
            mix_exp = (CPT) bn.getNode("Mix_Exp");
            if (mix_exp != null) {
                Mix_Exp = mix_exp.getVariable();
                System.out.println("Loaded " + Mix_Exp + " defined with domain " + Mix_Exp.getDomain());
            }
            mix_hm = (CPT) bn.getNode("Mix_HM");
            if (mix_hm != null) {
                Mix_HM = mix_hm.getVariable();
                System.out.println("Loaded " + Mix_HM + " defined with domain " + Mix_HM.getDomain());
            }
        }

        if (command.equalsIgnoreCase("TRAIN") && tsvmix1 != null && tsvmix2 != null) {

            TSVFile tsv1 = null;
            TSVFile tsv2 = null;
            String[] tf_names = null;
            String[] ct_names = null;
            Set<Object> ctypes1 = null;
            try {
                tsv1 = new TSVFile(tsvmix1, false); // RPKM
                tsv2 = new TSVFile(tsvmix2, false); // H3K27me3
                Set<Object> tfs1 = tsv1.getValues(1);
                Set<Object> tfs2 = tsv2.getValues(1);
                ctypes1 = tsv1.getValues(0);
                Set<Object> ctypes2 = tsv2.getValues(0);
                Set<String> tfs = new HashSet<>();
                Set<String> cts = new HashSet<>();
                for (Object obj : tfs1)
                    tfs.add((String) obj);
                for (Object obj : tfs2)
                    tfs.add((String) obj);
                for (Object obj : ctypes1)
                    cts.add((String) obj);
                for (Object obj : ctypes2)
                    cts.add((String) obj);
                tf_names = new String[tfs.size()];
                ct_names = new String[cts.size()];
                int i = 0;
                for (String obj : tfs) {
                    tf_names[i] = obj;
                    i++;
                }
                i = 0;
                for (Object obj : cts) {
                    ct_names[i] = obj.toString();
                    i++;
                }
                Arrays.sort(tf_names);
                Arrays.sort(ct_names);
            } catch (IOException e) {
                System.err.println(e);
                System.exit(1);
            }

            if (bnetmix1 != null && bnetmix2 != null) { // load and mix two models

                BNet bn_mix1 = BNBuf.load(bnetmix1); // RPKM
                System.out.println("Loaded variables and tables");
                // Nodes, and their variables
                GDT y1 = (GDT) bn_mix1.getNode("Signal");
                y1.setTrainable(false);
                Y1 = y1.getVariable();
                System.out.println("Loaded " + Y1 + " defined with domain " + Y1.getDomain());

                String mix2_file = args[2];
                BNet bn_mix2 = BNBuf.load(bnetmix2); // H3K27me3
                System.out.println("Loaded variables and tables");
                // Nodes, and their variables
                GDT y2 = (GDT) bn_mix2.getNode("Signal");
                y2.setTrainable(false);
                Y2 = y2.getVariable();
                System.out.println("Loaded " + Y2 + " defined with domain " + Y1.getDomain());

                System.out.println("Defining variables and tables: model generation " + GEN);
                TF = Predef.Nominal(tf_names, "TF");
                System.out.println("Defining " + TF + " with " + TF.getDomain().size() + " values");
                CObs = Predef.Nominal(ct_names, "CObs");
                System.out.println("Defining " + CObs + " with " + CObs.getDomain().size() + " values");
                Y1 = Predef.Real("Exp"); // re-define
                Y2 = Predef.Real("H3K27me3"); // re-define
                tf = new CPT(TF);
                ct = new CPT(CObs);
                Mix_HM = Predef.Boolean("Mix_HM");
                Mix_Exp = Predef.Boolean("Mix_Exp");
                // Define merged BN
                GDT yy1 = new GDT(Y1, Mix_Exp);
                GDT yy2 = new GDT(Y2, Mix_HM);
                for (Object v : Enumerable.bool.getValues()) {
                    yy1.put(new Object[]{v}, y1.getDistrib(new Object[]{v}));
                    yy2.put(new Object[]{v}, y2.getDistrib(new Object[]{v}));
                }
                yy1.setTrainable(false);
                yy2.setTrainable(false);
                bn = new BNet();
                if (MIXTURE_TF >= 2 && MIXTURE_CT >= 2) {
                    TF_lat = Predef.Number(MIXTURE_TF, "TF_lat");
                    tf_lat = new CPT(TF_lat, TF);
                    CT_lat = Predef.Number(MIXTURE_CT, "CT_lat");
                    ct_lat = new CPT(CT_lat, CObs);
                    mix_exp = new CPT(Mix_Exp, CT_lat, TF_lat);
                    mix_hm = new CPT(Mix_HM, CT_lat, TF_lat);
                    bn.add(ct, tf, ct_lat, tf_lat, mix_exp, mix_hm, yy1, yy2);
                } else if (MIXTURE_TF < 2 && MIXTURE_CT >= 2) {  // no latent TF node
                    CT_lat = Predef.Number(MIXTURE_CT, "CT_lat");
                    ct_lat = new CPT(CT_lat, CObs);
                    mix_exp = new CPT(Mix_Exp, CT_lat, TF);
                    mix_hm = new CPT(Mix_HM, CT_lat, TF);
                    bn.add(ct, tf, ct_lat, mix_exp, mix_hm, yy1, yy2);
                } else if (MIXTURE_CT < 2 && MIXTURE_TF >= 2) { // no latent cell type node
                    ;
                } else { // no latent nodes at all
                    mix_exp = new CPT(Mix_Exp, CObs, TF);
                    mix_hm = new CPT(Mix_HM, CObs, TF);
                    bn.add(ct, tf, mix_exp, mix_hm, yy1, yy2);
                }
                if (bnetmerge != null)
                    BNBuf.save(bn, bnetmerge + ".preset");

            }

            Object[][] rows1 = tsv1.getRows();
            Object[][] rows2 = tsv2.getRows();
            Object[][] all = new Object[rows2.length][4];
            try {
                for (int i = 0; i < rows2.length; i++) {
                    all[i][0] = rows2[i][0];
                    all[i][1] = rows2[i][1];
                    all[i][3] = Double.parseDouble(rows2[i][2].toString());
                    if (ctypes1.contains(rows2[i][0])) {
                        for (Object[] row1 : rows1) {
                            if (row1[0].toString().equalsIgnoreCase(rows2[i][0].toString()) && row1[1].toString().equalsIgnoreCase(rows2[i][1].toString()))
                                all[i][2] = Double.parseDouble(row1[2].toString());
                        }
                    }
                }
            } catch (RuntimeException e) {
                System.err.println(rows1[2] + "\t" +rows2[2]);
                e.printStackTrace();
            }
            LearningAlg em = new EM(bn);
            if (ROUNDS != -1)
                ((EM) em).setMaxRounds(ROUNDS);
            ((EM) em).setEMOption(EM_OPTION);
            em.train(all, new Variable[]{CObs, TF, Y1, Y2}, System.currentTimeMillis());
            if (bnetmerge != null)
                BNBuf.save(bn, bnetmerge);

        } else if (command.equalsIgnoreCase("QUERY")){ // QUERY

            List<BNode> nodes = bn.getOrdered();
            Object[][] qdata = DataBuf.load(inputfile, nodes);

            BNode qnode = bn.getNode(queryvar);
            if (qnode != null) {
                VarElim ve = new VarElim();
                ve.instantiate(bn);
                Variable qvar = qnode.getVariable();
                boolean qIsEnum = true;
                Object[] cols = null;
                try {
                    EnumVariable evar = (EnumVariable) qvar;
                    cols = new Object[evar.size()];
                    for (int i = 0; i < evar.size(); i ++)
                        cols[i] = evar.getDomain().get(i);
                } catch (ClassCastException e) { // query variable is non-enumerable
                    qIsEnum = false;
                    cols = new String[1];
                    cols[0] = qvar.getName();
                }

                int rowcnt = 0;
                Object[][] result = new Object[qdata.length + 1][];
                result[0] = new Object[nodes.size() + cols.length];
                for (int i = 0; i < nodes.size(); i ++)
                    result[0][i] = nodes.get(i).getName();
                for (int i = 0; i < cols.length; i ++)
                    result[0][i + nodes.size()] = (String)cols[i].toString();
                for (Object[] q : qdata) {
                    if (q.length != nodes.size()) {
                        System.err.println("Invalid query @ row " + (rowcnt + 1));
                        continue;
                    }
                    for (int col = 0; col < q.length; col++) {
                        BNode node = nodes.get(col);
                        if (q[col] != null && !node.getVariable().equals(qvar))
                            node.setInstance(q[col]);
                    }
                    Query query = ve.makeQuery(qvar);
                    CGTable r = (CGTable) ve.infer(query);
                    Distrib d = r.query(qvar);
                    Object[] save = new Object[q.length + cols.length];
                    for (int i = 0; i < q.length; i ++)
                        save[i] = q[i];
                    if (qIsEnum) {
                        for (int i = 0; i < cols.length; i ++)
                            save[q.length + i] = d.get(cols[i]);
                    } else { // qvar is non-enumerable
                        int nSamples = 10000;
                        double y = 0;
                        for (int j = 0; j < nSamples; j ++)
                            y += ((Double)d.sample());
                        save[q.length] = y / nSamples;
                    }
                    result[rowcnt + 1] = save;
                    rowcnt ++;
                }
                try {
                    TSVFile.saveObjects(outputfile, result);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                System.err.println("Invalid query variable");
            }
        } else {
            System.err.println("Invalid input parameters");
        }
    }

}
