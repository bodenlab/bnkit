package api;

import dat.file.Newick;
import dat.phylo.Tree;
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
            System.out.println("Problem host: " + u);
        } catch (IOException i) {
            System.out.println("Problem I/O:" + i);
        }
    }

    @AfterEach
    void shutdown_socket() {
        try {
            server_input.close();
            server_output.close();
            socket.close();
        } catch (UnknownHostException u) {
            System.out.println("Problem host: " + u);
        } catch (IOException i) {
            System.out.println("Problem I/O:" + i);
        }
    }

    @Test
    void fromJSON_queue() {

    }

    @Test
    void fromJSON_request() {
        JSONObject jreq1 = new JSONObject();
        jreq1.put("Command", "Label-tree");
        jreq1.put("Auth", "Guest");
        JSONObject params = new JSONObject();
        params.put("Tree", tree1.toJSON());
        jreq1.put("Params", params);
        GRequest req = GRequest.fromJSON(jreq1);
        System.out.println(req.toString() + "\t" + jreq1.toString());
        server_output.println(jreq1);
        String response = "";
        try {
            response = server_input.readLine();
            System.out.println("Server responded: " + response);
        } catch (IOException e) {
            e.printStackTrace();
        }
        //System.out.println("|" + response.substring(response.indexOf('\t') + 1) + "|");
        //assertTrue(new JSONObject(response.substring(response.indexOf('\t') + 1)).toString().equals(jreq1.toString()));
    }

    @Test
    void fromJSON_retrieve() {

    }

    @Test
    void fromJSON_cancel() {

    }

}