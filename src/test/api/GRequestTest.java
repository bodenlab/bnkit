package api;

import asr.Prediction;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.FastaReader;
import dat.file.Newick;
import dat.phylo.Tree;
import dat.pog.POGraph;
import json.JSONException;
import json.JSONObject;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.*;
import java.net.Socket;
import java.net.UnknownHostException;

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
        try {
            JSONObject jreq1 = new JSONObject();
            jreq1.put("Command", "Fake");
            jreq1.put("Auth", "Guest");
            JSONObject params = new JSONObject();
            params.put("Sleep", 10000);
            jreq1.put("Params", params);
            System.out.println("My server-request: " + jreq1);
            server_output.println(jreq1);
            String response = server_input.readLine();
            System.out.println("Server responded: " + response);
            response = server_input.readLine();
            System.out.println("Server responded: " + response);
            Thread.sleep(1000);
            JSONObject jreq2 = new JSONObject();
            jreq2.put("Job", 1);
            jreq2.put("Command", "Retrieve");
            System.out.println("My server-request: " + jreq2);
            server_output.println(jreq2);
            response = server_input.readLine();
            System.out.println("Server responded: " + response);
            Thread.sleep(5000);
            JSONObject jreq3 = new JSONObject();
            jreq3.put("Job", 1);
            jreq3.put("Command", "Status");
            System.out.println("My server-request: " + jreq3);
            server_output.println(jreq3);
            response = server_input.readLine();
            System.out.println("Server responded: " + response);
            Thread.sleep(5000);
            jreq3 = new JSONObject();
            jreq3.put("Job", 1);
            jreq3.put("Command", "Output");
            System.out.println("My server-request: " + jreq3);
            server_output.println(jreq3);
            response = server_input.readLine();
            System.out.println("Server responded: " + response);
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
            int idx108 = tree.getIndex("sequence108");
            int parent_108 = tree.getParent(idx108);
            Object label = tree.getLabel(parent_108);
            JSONObject jresult = jresponse.getJSONObject("Result");
            Prediction pred = Prediction.fromJSON(jresult.getJSONObject("Prediction"));
            POGraph pog = pred.getAncestor(label);
            System.out.println(pog);
            int[] supported = pog.getMostSupported();
            for (int i = 0; i < supported.length; i ++)
                System.out.println(pog.getNode(supported[i]).toJSON());
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

    }

    @Test
    void fromJSON_cancel() {

    }

}