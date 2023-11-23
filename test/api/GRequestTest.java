package api;

import asr.Prediction;
import bn.node.CPT;
import bn.node.GDT;
import bn.prob.GaussianDistrib;
import dat.*;
import dat.file.FastaReader;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.*;
import dat.pog.POGraph;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.*;

import static org.junit.jupiter.api.Assertions.assertTrue;

class GRequestTest {

    Tree tree1 = Newick.parse("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2");

    private Socket socket = null;
    private BufferedReader server_input = null;
    private PrintWriter server_output = null;

    @BeforeEach
    void setup_socket() {
       // establish a connection
        try {
            socket = new Socket("127.0.0.1", 4072);
            // stream FROM server
            server_input = new BufferedReader(new InputStreamReader(socket.getInputStream()));
            // stream TO server
            server_output = new PrintWriter(socket.getOutputStream(), true);
        } catch (UnknownHostException u) {
            System.err.println("Problem host: " + u);
        } catch (IOException i) {
            System.err.println("Problem I/O [note asr.GServer must be running before these tests]: " + i);
        }
    }

    @AfterEach
    void shutdown_socket() {
        try {
            server_input.close();
            server_output.close();
            socket.close();
        } catch (UnknownHostException u) {
            System.err.println("Problem host: " + u);
        } catch (IOException i) {
            System.err.println("Problem I/O:" + i);
        }
    }

    /**
     * Load an alignment from the project data directory
     * @param filename
     * @return
     */
    EnumSeq.Alignment loadAln(String filename) {
        try {
            FastaReader r = new FastaReader("data/"+filename, Enumerable.aacid, '-');
            EnumSeq.Gappy[] eseqs = r.loadGappy();
            return new EnumSeq.Alignment(eseqs);
        } catch (IOException e) {
            return null;
        }
    }

    /**
     * Load a phylogenetic treefrom the project data directory
     * @param filename
     * @return
     */
    Tree loadNwk(String filename) {
        try {
            Tree t = Newick.load("data/"+filename);
            return t;
        } catch (IOException e) {
            return null;
        }
    }

    @Test
    void fromJSON_queue() {
        int[] jobtimes = new int[] {4000, 3000, 2000, 1000};
        int[] jobnumbers = new int[jobtimes.length];
        try {
            for (int i = 0; i < jobtimes.length; i ++) {
                JSONObject jreq1 = new JSONObject();
                jreq1.put("Command", "Fake");
                jreq1.put("Auth", "Guest");
                JSONObject params = new JSONObject();
                params.put("Sleep", jobtimes[i]);
                jreq1.put("Params", params);
                System.out.println("My server-request: " + jreq1);
                server_output.println(jreq1);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                jobnumbers[i] = GMessage.fromJSON2Job(jresponse);
            }

            Thread.sleep(500);
            for (int i = 0; i < jobnumbers.length; i ++) {
                JSONObject jreq2 = new JSONObject();
                jreq2.put("Job", jobnumbers[i]);
                jreq2.put("Command", "Place");
                System.out.println("My server-request: " + jreq2);
                server_output.println(jreq2);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
                jreq2.put("Command", "Status");
                System.out.println("My server-request: " + jreq2);
                server_output.println(jreq2);
                jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
            }

            Thread.sleep(Arrays.stream(jobtimes).sum()); // wait until all jobs have been finished
            for (int i = 0; i < jobnumbers.length; i ++) {
                JSONObject jreq2 = new JSONObject();
                jreq2.put("Job", jobnumbers[i]);
                jreq2.put("Command", "Output");
                System.out.println("My server-request: " + jreq2);
                server_output.println(jreq2);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void fromJSON_request_recon() {
        try {
            Tree tree = loadNwk("66.nwk");
            EnumSeq.Alignment aln = loadAln("66.aln");
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Recon");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Alignment", aln.toJSON());
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            // System.out.println("Server responded: " + jresponse);
            Thread.sleep(2000); // waiting 5 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            int idx108 = tree.getIndex("sequence108");
            int parent_108 = tree.getParent(idx108);
            Object label = tree.getLabel(parent_108).toString();
            JSONObject jresult = jresponse.getJSONObject("Result");
            Map<Object, POGraph> ancestors = Prediction.fromJSONJustAncestors(jresult);
            POGraph pog = ancestors.get(label);
            // System.out.println(pog);
            int[] supported = pog.getMostSupported();
            String[] names = aln.getNames();
            int idx = -1;
            for (int i = 0; i < names.length; i ++) {
                if (names[i].equals("sequence108")) {
                    idx = i;
                    break;
                }
            }
            for (int i = 0; i < supported.length; i ++) {
                //System.out.println(pog.getNode(supported[i]).toJSON().get("Value") + "\t" + aln.getEnumSeq(idx).getStripped()[i]);
                assertTrue(pog.getNode(supported[i]).toJSON().get("Value").toString().toCharArray()[0] == (Character) aln.getEnumSeq(idx).getStripped()[i]);
            }
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    @Test
    void fromJSON_request_train() {
        IdxTree tree = IdxTree.fromJSON(new JSONObject("{\"Parents\":[-1,0,1,2,2,4,4,1,0,8,9,9,11,11,8,14,14,16,16],\"Labels\":[\"0\",\"1\",\"2\",\"S001\",\"3\",\"S002\",\"S003\",\"S004\",\"4\",\"5\",\"S005\",\"6\",\"S006\",\"S007\",\"7\",\"S008\",\"8\",\"S009\",\"S010\"],\"Distances\":[0,0.14,0.03,0.14,0.08,0.16,0.10,0.12,0.06,0.06,0.28,0.13,0.12,0.14,0.11,0.20,0.07,0.12,0.19],\"Branchpoints\":19}\n"));
        String[] headers = {   "S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"};
        // The values assigned to the tree above, tabulated as per headers; they were intended to group nicely per clade
        Double[][] rows1 = {{   3.63,  3.81,  2.89,  3.81,  2.54,  2.76,  3.79,  3.70,  1.94,  3.97}};

        try {
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Train");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(headers,rows1));
            params.put("States", new JSONArray(new Character[] {'A','B'}));
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void fromJSON_request_infer() {
        IdxTree tree = IdxTree.fromJSON(new JSONObject("{\"Parents\":[-1,0,1,2,2,4,4,1,0,8,9,9,11,11,8,14,14,16,16],\"Labels\":[\"0\",\"1\",\"2\",\"S001\",\"3\",\"S002\",\"S003\",\"S004\",\"4\",\"5\",\"S005\",\"6\",\"S006\",\"S007\",\"7\",\"S008\",\"8\",\"S009\",\"S010\"],\"Distances\":[0,0.14,0.03,0.14,0.08,0.16,0.10,0.12,0.06,0.06,0.28,0.13,0.12,0.14,0.11,0.20,0.07,0.12,0.19],\"Branchpoints\":19}\n"));
        String[] headers = {   "S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"};
        // The values assigned to the tree above, tabulated as per headers; they were intended to group nicely per clade
        Double[][] rows = {
                {   3.63,  3.81,  2.89,  3.81,  2.54,  2.76,  3.79,  3.70,  1.94,  3.97},
                {   3.33,  3.21,  2.93,  3.51,  2.59,  2.96,  3.49,  3.40,  2.24,  3.44}};

        try {
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Infer");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(headers,rows));
            params.put("States", new JSONArray(new Character[] {'A','B'}));
            params.put("Distrib", new JSONObject("{\"Condition\":[[\"A\"],[\"B\"]],\"Pr\":[[3.784926135903969,0.056738891699391655],[2.532458829359575,0.056738891699391655]],\"Index\":[0,1],\"Domain\":\"dat.Continuous@3cf71bc7\"}"));
            params.put("Leaves-only", true);
            params.put("Inference", "Marginal");
            params.put("Ancestor", 0);
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(2000); // waiting 2 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void fromJSON_retrieve() {
        try {
            Tree tree = loadNwk("66.nwk");
            EnumSeq.Alignment aln = loadAln("66.aln");
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Pogit");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Alignment", aln.toJSON());
            jreq1.put("Params", params);
            server_output.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            System.out.println(GMessage.fromJSON2Result(jresponse));
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void fromJSON_cancel() {
        int[] jobtimes = new int[] {4000, 3000, 2000, 1000};
        int[] jobnumbers = new int[jobtimes.length];

        try {
            for (int i = 0; i < jobtimes.length; i++) {
                JSONObject jreq1 = new JSONObject();
                jreq1.put("Command", "Fake");
                jreq1.put("Auth", "Guest");
                JSONObject params = new JSONObject();
                params.put("Sleep", jobtimes[i]);
                jreq1.put("Params", params);
                System.out.println("My server-request: " + jreq1);
                server_output.println(jreq1);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                jobnumbers[i] = GMessage.fromJSON2Job(jresponse);
            }

            Thread.sleep(500);
            // should succeed to cancel three (last) jobs, not the first, which should've started
            boolean[] success = new boolean[] {false, true, true, true};
            for (int i = jobnumbers.length - 1; i >= 0; i--) {
                JSONObject jreq2 = new JSONObject();
                jreq2.put("Job", jobnumbers[i]);
                jreq2.put("Command", "Cancel");
                System.out.println("My server-request: " + jreq2);
                server_output.println(jreq2);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
                assertTrue(jresponse.getBoolean("Cancel") == success[i]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    /*
    Test tree from Seb: Glycosyl hydrolase 13 family (stored in the ~/data folder
    Properties with good coverage include
    pH optimum: 56%
    TO: 62%
    Tn: 31%
    TS: 49%
     */
    JSONUtils.DataSet getProperty(String propname) throws IOException {
        TSVFile tsv = new TSVFile("data/3_2_1_1_annotations.txt", true);
        int col = tsv.getColumn(propname);
        if (col == -1)
            return null;
        Object[] vals = tsv.getCol(col);
        Object[] curated = new Object[vals.length];
        for (int i = 0; i < vals.length; i ++) {
            if (vals == null)
                continue;
            String sval = vals[i].toString();
            StringTokenizer stok = new StringTokenizer(sval, ";");
            String[] subsvals = new String[stok.countTokens()];
            boolean none = false;
            int j = 0;
            while (stok.hasMoreTokens()) {
                String tok = stok.nextToken().trim();
                if (tok.equals("None") || tok.equals("none") || tok.equals("Null") || tok.equals("null")) {
                    none = true;
                    break;
                }
                int stophere = tok.indexOf('_');
                if (stophere == -1)
                    subsvals[j] = tok;
                else
                    subsvals[j] = tok.substring(0, stophere);
                j += 1;
            }
            if (none)
                continue;
            Double[] subsvald = new Double[subsvals.length];
            double sum = 0;
            for (j = 0; j < subsvald.length; j ++) {
                try {
                    // could be a range
                    int stophere = subsvals[j].indexOf('-');
                    if (stophere == -1) // not a range
                        subsvald[j] = Double.parseDouble(subsvals[j]);
                    else { // range indeed, so pick middle point
                        String part2 = subsvals[j].substring(stophere + 1);
                        String part1 = subsvals[j].substring(0, stophere);
                        subsvald[j] = (Double.parseDouble(part2) + Double.parseDouble(part1)) / 2.0;
                    }
                } catch (NumberFormatException e) {
                    // not a number
                    throw new RuntimeException("Invalid format: value \"" + subsvals[j] + "\" on row " + (j + 1) + " is not a number, or range of numbers");
                }
                sum += subsvald[j];
            }
            curated[i] = sum / (double)subsvald.length;
        }
        String[] features = new String[vals.length];
        Object[][] values = new Object[1][vals.length];
        int entrycol = tsv.getColumn("Entry");
        for (int i = 0; i < vals.length; i ++) {
            features[i] = tsv.getCol(entrycol)[i].toString();
            values[0][i] = curated[i];
        }
        JSONUtils.DataSet ds = new JSONUtils.DataSet(features, values);
        return ds;
    }

    @Test
    void request_Sebs_tree() {
        try {
            Tree tree = Newick.load("data/3_2_1_1_filt.nwk");
            tree.save("data/glycosol_hydrolase.nwk", "ancestor");
            JSONUtils.DataSet ds = getProperty("BRENDA_PHO_DATA");
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Train");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(ds));
            params.put("States", new JSONArray(new Character[] {'A','B','C'}));
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult = jresponse.getJSONObject("Result");
            JSONObject jdistrib = jresult.getJSONObject("Distrib");

            JSONObject jreq2 = new JSONObject();
            jreq2.put("Command", "Infer");
            jreq2.put("Auth", "Guest");
            params.put("Distrib", jdistrib);
            params.put("Inference", "Marginal");
            for (int idx : tree.getAncestors()) {
                Integer ancid = (Integer) tree.getLabel(idx);
                params.put("Ancestor", ancid);
                params.put("Leaves-only", false); // switch to true if you want to see the parents distribution
                jreq2.put("Params", params);
                server_output.println(jreq2);
                System.out.println(jreq2);
                jresponse = new JSONObject(server_input.readLine());
                job = GMessage.fromJSON2Job(jresponse);
                System.out.println("Server responded: " + jresponse);
                Thread.sleep(500); // waiting 0.5 secs to make sure the job has finished
                JSONObject jreq2b = new JSONObject();
                jreq2b.put("Job", job);
                jreq2b.put("Command", "Output"); // request the output/result
                server_output.println(jreq2b);
                jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
            }
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    JSONUtils.DataSet getDataset(EnumSeq.Alignment aln, int col) {
        JSONUtils.DataSet ds = new JSONUtils.DataSet(aln.getNames(), new Object[][] {aln.getColumn(col)});
        return ds;
    }

    @Test
    void request_Train_1xmotif() {
        try {
            Tree tree = Newick.load("data/3_2_1_1_filt.nwk");
            tree.save("data/glycosol_hydrolase.nwk", "ancestor");
            EnumSeq.Alignment aln = loadAln("3_2_1_1_filt.aln");
            JSONUtils.DataSet ds = getDataset(aln, 709);
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Train");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(ds));
            params.put("States", new JSONArray(new Character[] {'A','B','C'}));
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult = jresponse.getJSONObject("Result");
            JSONObject jdistrib = jresult.getJSONObject("Distrib");

            JSONObject jreq2 = new JSONObject();
            jreq2.put("Command", "Infer");
            jreq2.put("Auth", "Guest");
            params.put("Distrib", jdistrib);
            params.put("Inference", "Marginal");
            for (int idx : tree.getAncestors()) {
                Integer ancid = (Integer) tree.getLabel(idx);
                params.put("Ancestor", ancid);
                params.put("Leaves-only", false); // switch to true if you want to see the parents distribution
                jreq2.put("Params", params);
                server_output.println(jreq2);
                System.out.println(jreq2);
                jresponse = new JSONObject(server_input.readLine());
                job = GMessage.fromJSON2Job(jresponse);
                System.out.println("Server responded: " + jresponse);
                Thread.sleep(1000); // waiting 1 sec to make sure the job has finished
                JSONObject jreq2b = new JSONObject();
                jreq2b.put("Job", job);
                jreq2b.put("Command", "Output"); // request the output/result
                server_output.println(jreq2b);
                jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
            }
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void request_Train_1xmotif_joint() {
        try {
            Tree tree = Newick.load("data/3_2_1_1_filt.nwk");
            tree.save("data/glycosol_hydrolase.nwk", "ancestor");
            EnumSeq.Alignment aln = loadAln("3_2_1_1_filt.aln");
            JSONUtils.DataSet ds = getDataset(aln, 709);
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Train");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(ds));
            params.put("States", new JSONArray(new Character[] {'a','b','c','d'}));
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult = jresponse.getJSONObject("Result");
            JSONObject jdistrib = jresult.getJSONObject("Distrib");

            JSONObject jreq2 = new JSONObject();
            jreq2.put("Command", "Infer");
            jreq2.put("Auth", "Guest");
            params.put("Distrib", jdistrib);
            params.put("Inference", "Joint");
            params.put("Leaves-only", true); // switch to true if you want to see the parents distribution
            jreq2.put("Params", params);
            server_output.println(jreq2);
            System.out.println(jreq2);
            jresponse = new JSONObject(server_input.readLine());
            job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            JSONObject jreq2b = new JSONObject();
            jreq2b.put("Job", job);
            jreq2b.put("Command", "Output"); // request the output/result
            server_output.println(jreq2b);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult2 = jresponse.getJSONObject("Result");
            JSONUtils.DataSet dspred = JSONUtils.DataSet.fromJSON(jresult2.getJSONObject("Predict"));
            TreeInstance[] multi = TreeInstance.createFromDataset(tree, dspred.getFeatures(), dspred.getNonitemisedData());
            for (TreeInstance ti : multi)
                ti.save("data/3_2_1_1_pred.nwk");
            TSVFile.print(new Object[][] {dspred.getFeatures()});
            TSVFile.print(dspred.getNonitemisedData());
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void request_TrainModes_motif_joint() {
        IdxTree mblTree = null;
        EnumSeq.Alignment mblAln = null;
        Tree tree;
        try {
            tree = Tree.load("data/master_relabel.nwk", "newick");
            mblTree = (IdxTree)tree;
            FastaReader r = new FastaReader("data/master_t10.aln", Enumerable.aacid, '-');
            EnumSeq.Gappy[] eseqs = r.loadGappy();
            mblAln = new EnumSeq.Alignment(eseqs);

            /* Metal binding at (positions in alignment, starting at column 1)
            Alpha site: 184, 186, 306
            Beta site: 188, 189, 567
            Bridging residue: 332
             */
            int[] sites = new int[] {183, 185, 305};//,  187, 188, 566};
            long SEED = 3;
            JSONUtils.DataSet ds = new JSONUtils.DataSet(mblAln.getNames(), new String[] {"Pos184", "Pos186", "Pos306"}, new Object[][][] {mblAln.getColumns(sites)});
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "TrainModes");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(ds));
            PhyloPlate.Modes template = new PhyloPlate.Modes(new Enumerable[] {new Enumerable(new Object[] {'A','B','C'})}); //,new Enumerable(new Object[] {'a','b','c'})});
            PhyloPlate phybn = new PhyloPlate(mblTree, template);
            PhyloPlate.Plate plate1 = phybn.getPlate(0);
            CPT[] masters = new CPT[sites.length];
            for (int s = 0; s < sites.length; s ++) {
                CPT pos = new CPT(new EnumVariable(Enumerable.aacid, "Pos" + (sites[s] + 1)), plate1.getParents(new int[] {0})); //,1}));
                plate1.addNode(pos);
            }
            params.put("Distrib", plate1.toJSON());
            params.put("Seed", SEED);
            params.put("Gamma", 1.0);
            params.put("Rate", 1.0);
            params.put("Rounds", 10); // training rounds
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);

            Thread.sleep(20000); // waiting 20 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult = jresponse.getJSONObject("Result");
            JSONObject jdistrib = jresult.getJSONObject("Distrib");

            JSONObject jreq2 = new JSONObject();
            jreq2.put("Command", "InferModes");
            jreq2.put("Auth", "Guest");
            params.put("Distrib", jdistrib);
            // params.put("Inference", "Marginal");
            params.put("Inference", "Marginal");
            params.put("Leaves-only", false);
            params.put("Latent", true);
            params.put("Queries", new JSONArray(new Object[] {0, 1, 2, "Q704V1"}));
            jreq2.put("Params", params);
            server_output.println(jreq2);
            System.out.println(jreq2);
            jresponse = new JSONObject(server_input.readLine());
            job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(5000); // waiting 5 secs to make sure the job has finished
            JSONObject jreq2b = new JSONObject();
            jreq2b.put("Job", job);
            jreq2b.put("Command", "Output"); // request the output/result
            server_output.println(jreq2b);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult2 = jresponse.getJSONObject("Result");
            // JSONUtils.DataSet dspred = JSONUtils.DataSet.fromJSON(jresult2.getJSONObject("Predict"));
 /*
            TreeInstance[] multi = TreeInstance.createFromDataset(tree, dspred.getFeatures(), dspred.getNonitemisedData());
            for (TreeInstance ti : multi)
                ti.save("data/3_2_1_1_pred.nwk");
            TSVFile.print(new Object[][] {dspred.getFeatures()});
            TSVFile.print(dspred.getNonitemisedData());
 */
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    @Test
    void request_TrainModes_real() {
        long SEED = 2;
        IdxTree tree = Newick.parse("(D:0.3,((A:0.1,B:0.1):0.1,C:0.2):0.1)");
        try {
            JSONUtils.DataSet ds = new JSONUtils.DataSet(tree.getNames(), new String[] {"Prop1"}, new Object[][][] {
                    // sample 1:
                    new Object[][] {
                            // item 1 (D):
                            {/* prop 1 */ 1.1},
                            // item 2 (A):
                            {/* prop 1 */ 3.2},
                            // item 3 (B)
                            {/* prop 1 */ null},
                            // item 4 (C)
                            {/* prop 1 */ 3.3}},
                    // sample 2:
                    new Object[][] {
                            // item 1 (D):
                            {/* prop 1 */ 1.5},
                            // item 2 (A):
                            {/* prop 1 */ null},
                            // item 3 (B)
                            {/* prop 1 */ 3.0},
                            // item 4 (C)
                            {/* prop 1 */ 3.4}}});
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "TrainModes");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            System.out.println(tree);
            params.put("Dataset", JSONUtils.toJSON(ds));
            PhyloPlate.Modes template = new PhyloPlate.Modes(new Enumerable[] {new Enumerable(new Object[] {'a','b'})});
            PhyloPlate phybn = new PhyloPlate(tree, template);
            PhyloPlate.Plate plate1 = phybn.getPlate(0);
            Variable var = new Variable(new Continuous(), "Prop1");
            GDT gdt = new GDT(var, plate1.getParents(new int[] {0})); //,1}));
            plate1.addNode(gdt);
            params.put("Distrib", plate1.toJSON());
            params.put("Seed", SEED);
            params.put("Gamma", 1.0);
            params.put("Rate", 1.0);
            params.put("Rounds", 10); // training rounds
            jreq1.put("Params", params);
            server_output.println(jreq1);
            System.out.println(jreq1);
            JSONObject jresponse = new JSONObject(server_input.readLine());
            int job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);

            Thread.sleep(2000); // waiting 2 secs to make sure the job has finished
            jreq1 = new JSONObject();
            jreq1.put("Job", job);
            jreq1.put("Command", "Output"); // request the output/result
            server_output.println(jreq1);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult = jresponse.getJSONObject("Result");
            JSONObject jdistrib = jresult.getJSONObject("Distrib");

            JSONObject jreq2 = new JSONObject();
            jreq2.put("Command", "InferModes");
            jreq2.put("Auth", "Guest");
            params.put("Distrib", jdistrib);
            params.put("Inference", "Marginal");
            params.put("Leaves-only", false);
            params.put("Latent", true);
            params.put("Queries", new JSONArray(new Object[] {0, 1, 2, "B"}));
            jreq2.put("Params", params);
            server_output.println(jreq2);
            System.out.println(jreq2);
            jresponse = new JSONObject(server_input.readLine());
            job = GMessage.fromJSON2Job(jresponse);
            System.out.println("Server responded: " + jresponse);
            Thread.sleep(1000); // waiting 1 sec to make sure the job has finished
            JSONObject jreq2b = new JSONObject();
            jreq2b.put("Job", job);
            jreq2b.put("Command", "Output"); // request the output/result
            server_output.println(jreq2b);
            jresponse = new JSONObject(server_input.readLine());
            System.out.println("Server responded: " + jresponse);
            JSONObject jresult2 = jresponse.getJSONObject("Result");
        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    @Test
    void request_InferModes_Lewis() {
        String folder = "/Users/mikael/simhome/ASR/ReconMode/";
        try {
            Set<String> extants = new HashSet<>();
            Set<String> ancestors = new HashSet<>();
            TSVFile names = new TSVFile(folder + "ERED_names.tsv", true);
            Map<String, String> namemap = new HashMap<>();
            for (Object[] row : names.getRows())
                namemap.put((String)row[1], (String)row[4]);
            TSVFile wt = new TSVFile(folder + "s2-2.tsv", true);
            Map<String, Object[]> propmap = new HashMap<>();
            for (Object[] row : wt.getRows()) {
                propmap.put((String) row[0], row);
                extants.add((String) row[0]);
            }
            TSVFile anc = new TSVFile(folder + "s2-1.tsv", true);
            for (Object[] row : anc.getRows()) {
                propmap.put((String) row[0], row);
                ancestors.add((String) row[0]);
            }
            String[] extants_arr = new String[extants.size()];
            String[] ancestors_arr = new String[ancestors.size()];
            extants.toArray(extants_arr);
            ancestors.toArray(ancestors_arr);
            IdxTree tree = Newick.load(folder + "lewis_tree.nwk");

            int NPROTEINS = 20;
            int NCOMPONENTS = 3;
            int EM_ROUNDS = 25;
            int[] percent_ancestors = new int[] {0, 20, 40, 60, 80, 100};
            Map<Integer, JSONObject> jobs = new HashMap<>();
            for (int percent : percent_ancestors) {
                int NANCESTORS = NPROTEINS * percent / 100;
                int NEXTANTS = NPROTEINS - NANCESTORS;
                for (int SEED = 0; SEED < 5; SEED += 1) {
                    Random rand = new Random(SEED);
                    Set<String> select = new HashSet<>();
                    while (select.size() < NEXTANTS)
                        select.add(extants_arr[rand.nextInt(extants_arr.length)]);
                    while (select.size() < NPROTEINS)
                        select.add(ancestors_arr[rand.nextInt(ancestors_arr.length)]);
                    String[] select_arr = new String[select.size()];
                    select.toArray(select_arr);
                    Object[][][] data = new Object[1][select_arr.length][1];
                    for (int i = 0; i < select_arr.length; i ++) {
                        data[0][i][0] = propmap.get(select_arr[i])[1] ;
                        if (data[0][i][0] != null)
                            data[0][i][0] = (((double)((Integer)data[0][i][0]) + i/100.0));
                    }
                    Set<Object> queryset = new HashSet<>();
                    for (Object candidate : extants)
                        if (!select.contains(candidate))
                            queryset.add(candidate);
                    for (Object candidate : ancestors)
                        if (!select.contains(candidate))
                            queryset.add(candidate);
                    String[] unselect = new String[queryset.size()];
                    queryset.toArray(unselect);
                    JSONUtils.DataSet ds = new JSONUtils.DataSet(select_arr, new String[] {"Tm"}, data);
                    GDT gdt = GDT.trainGMM(ds, NCOMPONENTS, EM_ROUNDS, SEED); // this could come from a PhyloPlate BN as well
                    JSONObject jreq1 = new JSONObject();
                    jreq1.put("Command", "InferModes");
                    jreq1.put("Auth", "Guest");
                    JSONObject params = new JSONObject();
                    params.put("Tree", tree.toJSON());
                    params.put("Leaves-only", false);
                    params.put("Dataset", JSONUtils.toJSON(ds));
                    PhyloPlate.Modes template = new PhyloPlate.Modes(new Enumerable[] {gdt.getParents().get(0).getDomain()});
                    PhyloPlate phybn = new PhyloPlate(tree, template);
                    PhyloPlate.Plate plate1 = phybn.getPlate(0);
                    plate1.addNode(gdt);
                    params.put("Distrib", plate1.toJSON());
                    params.put("Gamma", 1.0);
                    params.put("Rate", 1.0);
                    params.put("Inference", "Marginal");
                    params.put("Leaves-only", false);
                    params.put("Latent", false);
                    params.put("Queries", new JSONArray(unselect));
                    jreq1.put("Params", params);
                    server_output.println(jreq1);
                    System.out.println(jreq1);
                    JSONObject jresponse = new JSONObject(server_input.readLine());
                    int job = GMessage.fromJSON2Job(jresponse);
                    System.out.println("Server responded: " + jresponse + " for job " + job);
                    jobs.put(job, params);
                }
            }
            // done submitting jobs... wait for them to finish
            Map<Integer, JSONObject> results = new HashMap<>();
            for (Map.Entry<Integer, JSONObject> entry : jobs.entrySet()) {
                Integer job = entry.getKey();
                JSONObject jreq2 = new JSONObject();
                jreq2.put("Job", job);
                jreq2.put("Command", "Output");
                System.out.println("My server-request: " + jreq2);
                server_output.println(jreq2);
                JSONObject jresponse = new JSONObject(server_input.readLine());
                System.out.println("Server responded: " + jresponse);
                JSONObject result = GMessage.fromJSON2Result(jresponse);
                while (result == null) {
                    Thread.sleep(200); // waiting 0.2 secs before checking again job has finished
                    server_output.println(jreq2);
                    jresponse = new JSONObject(server_input.readLine());
                    System.out.println("Server responded: " + jresponse);
                    result = GMessage.fromJSON2Result(jresponse);
                }
                results.put(job, result);
            }
            // all results have been collected, now sort them into spreads
            Object[][] out = new Object[propmap.size() + 1][jobs.size() + 2];
            Map<String, Map<Integer, Object>> outmap = new HashMap<>();
            for (Map.Entry<Integer, JSONObject> entry : jobs.entrySet()) {
                Integer job = entry.getKey();
                JSONObject params = entry.getValue();
                JSONUtils.DataSet inp_ds = JSONUtils.DataSet.fromJSON(params.getJSONObject("Dataset"));
                String[] inp_names = inp_ds.getItems();
                for (int i = 0; i < inp_names.length; i ++) {
                    if (!outmap.containsKey(inp_names[i]))
                        outmap.put(inp_names[i], new HashMap<>());
                    Map<Integer, Object> valmap = outmap.get(inp_names[i]);
                    valmap.put(job, inp_ds.getItemisedData()[0][i][0] != null ? -(Double)inp_ds.getItemisedData()[0][i][0] : null);
                }
                JSONObject result = results.get(job);
                JSONUtils.DataSet res_ds = JSONUtils.DataSet.fromJSON(result.getJSONObject("Predict"));
                String[] out_names = res_ds.getItems();
                for (int i = 0; i < out_names.length; i ++) {
                    if (!outmap.containsKey(out_names[i]))
                        outmap.put(out_names[i], new HashMap<>());
                    Map<Integer, Object> valmap = outmap.get(out_names[i]);
                    valmap.put(job, res_ds.getItemisedData()[0][i][0]);
                }
            }
            Object[][] trn = new Object[out.length][out[0].length];
            Object[][] tst = new Object[out.length][out[0].length];
            List<String> alphabetical = new ArrayList<>();
            alphabetical.addAll(outmap.keySet());
            Collections.sort(alphabetical);
            List<Integer> jobsascend = new ArrayList<>();
            jobsascend.addAll(jobs.keySet());
            Collections.sort(jobsascend);
            int row = 1;
            for (String name : alphabetical) {
                int col = 2;
                out[row][0] = name;
                tst[row][0] = name;
                trn[row][0] = name;
                out[row][1] = propmap.get(name)[1] == null ? null : Double.parseDouble(propmap.get(name)[1].toString()) ;
                tst[row][1] = propmap.get(name)[1] == null ? null : Double.parseDouble(propmap.get(name)[1].toString()) ;
                trn[row][1] = propmap.get(name)[1] == null ? null : Double.parseDouble(propmap.get(name)[1].toString()) ;
                for (int job : jobsascend) {
                    if (row == 1) {
                        out[0][col] = job;
                        tst[0][col] = job;
                        trn[0][col] = job;
                    }
                    Map<Integer, Object> myjobs = outmap.get(name);
                    Object y = myjobs.get(job);
                    if (y != null) {
                        if ((Double) y < 0)
                            trn[row][col] = -(Double) y;
                        else
                            tst[row][col] = y;
                    }
                    out[row][col ++] = y;
                }
                row ++;
            }
            TSVFile tsv_tst = new TSVFile(tst, true);
            TSVFile tsv_trn = new TSVFile(trn, true);
            for (int i = 1; i <= jobs.size(); i ++) {
                TSVFile.DEFAULT_SHAPE = 2;
                TSVFile.DEFAULT_FILL = 1;
                TSVFile.DEFAULT_SIZE = 3;
                TSVFile.save2iTOL(folder + "test_" + i + ".itol", tsv_tst.getCol(0), tsv_tst.getCol(i+1), "Tm(test)", 12, 35., 65.);
                TSVFile.DEFAULT_SHAPE = 4;
                TSVFile.save2iTOL(folder + "trn_" + i + ".itol", tsv_trn.getCol(0), tsv_trn.getCol(i+1), "Tm(trn)", 12, 35., 65.);
                TSVFile.DEFAULT_SHAPE = 3;
                TSVFile.DEFAULT_FILL = 1;
                TSVFile.DEFAULT_SIZE = 5;
                TSVFile.save2iTOL(folder + "exp_" + i + ".itol", tsv_trn.getCol(0), tsv_trn.getCol(1), "Tm(exp)", 12, 35., 65.);
                Object[] match = new Boolean[tsv_tst.getCol(0).length];
                Double[] myexp = new Double[match.length];
                Double[] mytst = new Double[match.length];
                int cnt_match = 0;
                for (int j = 0; j < match.length; j ++) {
                    myexp[j] = (Double) tsv_trn.getCol(1)[j];
                    mytst[j] = (Double) tsv_tst.getCol(i + 1)[j];
                    if (mytst[j] != null && myexp[j] != null) {
                        match[j] = (Math.abs(myexp[j] - mytst[j]) / myexp[j] < 0.05);
                        if ((Boolean)match[j])
                            cnt_match++;
                    }
                }
                System.out.println("Job index " + i + " has " + cnt_match + " matches");
                TSVFile.save2iTOL(folder + "mat_" + i + ".itol", tsv_trn.getCol(0), match, "Match", 12, 35., 65.);

            }
            TSVFile.saveObjects(folder + "results.tsv", out);

        } catch (JSONException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}