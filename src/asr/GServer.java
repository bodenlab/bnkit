package asr;

import api.CommandCentral;
import api.GRequest;
import api.GMessage;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.net.*;
import java.io.*;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

public class GServer {

    static int MAX_CLIENTS = 20; // maximum number of clients that can connect to server

    private ServerSocket serverSocket;
    protected Map<ClientHandler, Date> clients = new HashMap<>();
    private final GRequest.JobQueue queue;
    private final CommandCentral commandCentral;

    /**
     * Open the socket on the server using specified port
     * @param port
     * @throws IOException
     */
    public void start(int port) throws IOException {
        serverSocket = new ServerSocket(port);
        while (true) { // this runs forever
            Socket client = serverSocket.accept(); // accept new client at the socket
            ClientHandler handler = new ClientHandler(client);
            // before we start a thread, make sure we are not exceeding maximum number of clients
            if (clients.size() < MAX_CLIENTS) {
                clients.put(handler, new Date());
                handler.start(); // spawn-off client thread
            } else {
                // deny service to client
                handler.deny();
            }
        }
    }

    /**
     * Terminate the socket running on the server
     * @throws IOException
     */
    public void stop() throws IOException {
        serverSocket.close();
    }

    /**
     * Create a server instance
     */
    public GServer(String directory) {
        // the command central is where actual jobs are managed; for now, we just create the instance
        commandCentral = new CommandCentral(this);
        // as requests come in (via the command central) they may be queued; here we create the job queue
        // TODO: check that the directory exists and/or is writeable

        queue = new GRequest.JobQueue(new File(directory)); // the directory is where completed jobs are stored
        // the job queue is a separate thread (which monitors and dispatches jobs as they pile up), which we start here
        queue.start();
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
        System.out.println("asr.GServer [-p <port>] [-d <dir>] [-c <max-clients>] [-h]");
        System.out.println("\twhere <port> is the port number that this server is listening to (default is 4072)");
        System.out.println("\t<dir> is the directory in which completed jobs are stored temporarily (default is /tmp)");
        System.out.println("\t<max-clients> is the maximum number of clients that will be allowed (default is 20)");
        System.out.println("\t-h prints this help screen");
        if (error > 0) {
            System.err.println("Error " + error + ": " + message);
        }
        System.exit(error);
    }

    /**
     * Main commandline front-end of GRASP server.
     * Program operates by accepting requests via sockets.
     * @param args see usage instructions
     */
    public static void main(String args[]) {
        Integer SERVERPORT = 4072; // default UQ's post-code
        String  DIRECTORY = "/tmp/";

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("p") && args.length > a + 1) {
                    try {
                        SERVERPORT = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(1, args[a] + "port must be an integer 0 - 65535 (default is 4072)");
                    }
                } else if (arg.equalsIgnoreCase("c") && args.length > a + 1) {
                    try {
                        MAX_CLIENTS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(3, args[a] + "max-clients must be a positive integer 1 or above (default is 20)");
                    }
                    if (MAX_CLIENTS < 1)
                        usage(3, args[a] + "max-clients must be a positive integer 1 or above (default is 20)");
                } else if (arg.equalsIgnoreCase("d") && args.length > a + 1) {
                    DIRECTORY = args[++a];
                } else if (arg.equalsIgnoreCase("h")) {
                    usage();
                }
            }
        }
        GServer server = new GServer(DIRECTORY);
        System.out.println("Server initialises a socket open for requests");
        try {
            server.start(SERVERPORT);
            System.out.println("Server proceeds to terminate its socket");
            server.stop();
        } catch (IOException e) {
            usage(2, e.getMessage());
        }
    }

    /**
     * Inner class that is responsible to handle the communications with a specific client.
     * If client is not denied service, the server class (above) will "start" the "run" method (hence a separate thread).
     */
    private class ClientHandler extends Thread {
        private Socket clientSocket;
        private PrintWriter out;   // the stream to which messages to the client are written
        private BufferedReader in; // the stream from which messages from the client is read

        /**
         * Start to handle a client on a socket.
         * @param socket the reference to the client
         */
        public ClientHandler(Socket socket) {
            this.clientSocket = socket;
            System.out.println("Client connected: \"" + clientSocket.toString() + "\"");
            try {
                OutputStream out_stream = clientSocket.getOutputStream();
                InputStream in_stream = clientSocket.getInputStream();
                out = new PrintWriter(out_stream, true);
                in = new BufferedReader(new InputStreamReader(in_stream));
            } catch (IOException e) {
                System.err.println("In thread; client socket thread failed: \"" + clientSocket.toString() + "\"");
            }
        }

        public void deny() {
            System.out.println("[" + clientSocket.toString() + "]: Denied client, server's closing their socket and associated thread");
            out.println(GMessage.errorToJSON(5, "Denied, since max number of clients (" + MAX_CLIENTS + ") reached on server"));
            try {
                in.close();
                out.close();
                clientSocket.close();
            } catch (IOException e) {
                System.err.println("Closing; client socket thread failed: \"" + clientSocket.toString() + "\"");
            }
        }

        /**
         * This is the method that is called if server-to-client socket thread is "started".
         * It receives messages from the client, which are parsed and made into "requests"
         * The first type of requests are simply concerned with the server, e.g. queries about queue.
         * The method has access to the command central instance, which is used to manage
         * jobs that are initiated by the client (the second, important type of requests).
         */
        public void run() { // client-specific thread running on server
            try {
                System.out.println("[" + clientSocket.toString() + "] client connected");
                String inputLine = in.readLine();
                while (inputLine != null) {
                    System.out.println("[" + clientSocket.toString() + "] message (50 chars max): " + inputLine.substring(0, Math.min(inputLine.length(), 50)) + ((inputLine.length() > 50) ? "..." : ""));
                    try {
                        // client message is parsed as JSON, and an exception is generated if it fails (+ error passed to client)
                        JSONObject json = new JSONObject(inputLine);
                        // Distinguish between Server/Job related commands (first type, when there's mention of "Job")
                        int job = GMessage.fromJSON2Job(json);
                        if (job > 0) { // this command has job number, so needs to be handled by server-aware code
                            GRequest request = queue.getRequest(job);               // retrieve from current job queue
                            if (request != null) {                                  // found job...
                                String command = json.optString("Command", null);    // extract command
                                if (command == null) {
                                    out.println(GMessage.errorToJSON(3, "No command given about job " + job));
                                } else {
                                    if (command.equals("Retrieve")) {
                                        // System.out.println("Sending this to client: " + request.toJSON());
                                        out.println(request.toJSON());                  // pass on job info
                                    } else if (command.equals("Output")) {
                                        JSONObject result = request.getResult();        // request to retrieve output of completed job
                                        if (result != null) {                           // if available... probably in storage
                                            // System.out.println("Sending this to client: " + result);
                                            out.println(GMessage.server2clientReJob(job, "Result", result));
                                        } else {
                                            out.println(GMessage.errorToJSON(2));       // job not available; pass error
                                        }
                                    } else if (command.equals("Status")) {              // what's the job doing
                                        // System.out.println("Sending this to client: " + request.getStatus());
                                        out.println(GMessage.server2clientReJob(job, "Status", request.getStatus()));
                                    } else if (command.equals("Place")) {
                                        // System.out.println("Sending this to client: " + queue.getPlace(request));
                                        out.println(GMessage.server2clientReJob(job, "Place", queue.getPlace(request)));
                                    } else if (command.equals("Cancel")) {
                                        out.println(GMessage.server2clientReJob(job, "Cancel", queue.cancel(request)));
                                    } else {
                                        out.println(GMessage.errorToJSON(3, "Invalid command: " + command));
                                    }
                                }
                            } else
                                out.println(GMessage.errorToJSON(2, "No such job in queue: " + job));
                        } else { // no "Job"
                            String command = json.optString("Command", "");    // extract command
                            if (command.equals("Status")) {
                                JSONObject jreport = new JSONObject();
                                JSONArray jtable = new JSONArray();
                                int cnt = 1;
                                for (GRequest greq : queue.getJobs()) {
                                    JSONObject jjob = greq.job2JSON();
                                    jjob.put("Place", greq.isWaiting() ? cnt++ : 0);
                                    jtable.put(jjob);
                                }
                                jreport.put("Jobs", jtable);
                                jreport.put("Clients", clients.size());
                                out.println(jreport);
                            } else {
                                // Second type of commands: a new, actual compute job so needs to be managed and may be queued
                                try {
                                    GRequest greq = commandCentral.createRequest(json);
                                    // System.out.println("OK--command detected: " + greq);
                                    // decision: run it now, or queue it for later...
                                    if (greq.isQueued()) {      // the request requires to be queued
                                        queue.add(greq);        // add to job queue
                                        // System.out.println("Informing client: Job " + greq.getJob() + " has been dispatched to queue, cancel with {\"Job\":\"" + greq.getJob() + "\",\"Command\":\"Cancel\"}");
                                        // TODO: find more info about what resources are required for job so client (and compute) can be advised
                                        out.println(GMessage.server2clientReJob(greq.getJob(), "Queued"));
                                    } else { // run now ... not in separate thread
                                        greq.runnow();
                                        JSONObject result = greq.getResult();           // request to retrieve output of completed job
                                        if (result != null) {                           // if available... probably in storage
                                            out.println(GMessage.server2clientReJob(job, "Result", result));
                                        } else {
                                            out.println(GMessage.errorToJSON(2));       // job not available; pass error
                                        }
                                    }
                                } catch (CommandCentral.GRequestRuntimeException syntax) {
                                    out.println(GMessage.errorToJSON(4, syntax.getMessage()));
                                }
                            }
                        }
                        // System.out.println("Accepts new request");
                    } catch (JSONException e) {
                        out.println(GMessage.errorToJSON(1));
                        System.out.println("[" + clientSocket.toString() + "]: " + GMessage.errorToJSON(1));
                    }
                    inputLine = in.readLine();
                }
                in.close();
                clientSocket.close();
            } catch (IOException e) {
                System.err.println("[" + clientSocket.toString() + "]: socket thread failed");
            }
            System.out.println("[" + clientSocket.toString() + "]: closing socket and associated thread");
            out.close();
            clients.remove(this);
        }
    }
}
