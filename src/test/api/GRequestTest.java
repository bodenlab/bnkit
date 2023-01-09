package api;

import asr.Prediction;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.FastaReader;
import dat.file.Newick;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.pog.POGraph;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Arrays;
import java.util.Map;

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
        Double[][] rows1 = {{   3.63,  3.81,  2.89,  3.81,  2.54,  2.76,  3.79,  3.70,  1.94,  3.97}};

        try {
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Infer");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Tree", tree.toJSON());
            params.put("Dataset", JSONUtils.toJSON(headers,rows1));
            params.put("States", new JSONArray(new Character[] {'A','B'}));
            params.put("Distrib", new JSONObject("{\"Condition\":[[\"A\"],[\"B\"]],\"Pr\":[[3.784926135903969,0.056738891699391655],[2.532458829359575,0.056738891699391655]],\"Index\":[0,1],\"Domain\":\"dat.Continuous@3cf71bc7\"}"));
            params.put("Inference", "Marginal");
            params.put("Ancestor", 5);
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

}