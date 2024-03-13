package asr;

import api.CommandCentral;
import api.GMessage;
import api.GRequest;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class GCommand {

    /**
     * Create a GRASP Command Line instance
     */
    public GCommand(JSONObject json) {
        CommandCentral central = new CommandCentral();
        try {
            GRequest greq = central.createRequest(json);
            if (greq != null) {
                greq.runnow();
                if (greq.isFailed())
                    System.err.println("Failed to complete request because: " + greq.getError());
                JSONObject result = greq.getResult();               // request to retrieve output of completed job
                if (result != null) {                               // if available...
                    System.out.println(result.toString());
                } else {
                    System.out.println(GMessage.errorToJSON(2));    // job not available; pass error
                }
            } else { // command not recognised
                System.err.println("Failed to identify command in request");
            }
        } catch (CommandCentral.GRequestRuntimeException syntax) {
            System.out.println(GMessage.errorToJSON(4, syntax.getMessage()));
        }
    }

    /**
     * Print usage instructions without error
     */
    public static void usage() {
        usage(0, null);
    }

    /**
     * Print usage instructions with (optional) error
     */
    public static void usage(int error, String message) {
        System.out.println("Usage:");
        System.out.println("asr.GCommand [<command>]");
        System.out.println("\tcommand to be executed, with parameters provided on the standard input in JSON format");
        if (error > 0) {
            System.err.println("Error " + error + ": " + message);
        }
        System.exit(error);
    }

    /**
     * Commandline front-end of GRASP JSON requests.
     * Program operates by accepting a single job request via JSON arguments.
     *
     * @param args see usage instructions
     */
    public static void main(String args[]) {
        String command = null;
        String auth = "Command-line user";
        String json = null;
        JSONObject jobj = null;
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {

            }
        }
        BufferedReader br = null;
        try { // read stdin
            br = new BufferedReader(new InputStreamReader(System.in));
            StringBuilder sb = new StringBuilder();
            String input = br.readLine();
            while (input != null) {
                sb.append(input);
                input = br.readLine();
            }
            json = sb.toString();
            jobj = new JSONObject(json);
        } catch (IOException e) {
            System.err.println("Error in standard input: " + e.getMessage());
        } catch (JSONException e) {
            System.err.println("Error in JSON: " + e.getMessage());
        }
        if (jobj != null) {
            new GCommand(jobj);
        }
    }
}

