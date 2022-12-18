package asr;

import api.GRequest;
import json.JSONException;
import json.JSONObject;

import java.net.*;
import java.io.*;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class GServer {

    static int MAX_CLIENTS = 20; // maximum number of clients that can connect to server

    private ServerSocket serverSocket;
    protected Map<ClientHandler, Date> clients = new HashMap<>();
    final GRequest.JobQueue queue;

    /**
     * Open the socket on the server using specified port
     * @param port
     * @throws IOException
     */
    public void start(int port) throws IOException {
        serverSocket = new ServerSocket(port);
        while (true) {
            Socket client = serverSocket.accept();
            ClientHandler handler = new ClientHandler(client);
            if (clients.size() < MAX_CLIENTS) {
                clients.put(handler, new Date());
                handler.start(); // spawn-off client thread
            } else {
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
    public GServer() {
        queue = new GRequest.JobQueue();
    }


    public static void main(String args[]) {
        Integer SERVERPORT = 4072; // default UQ's post-code
        if (args.length == 1) {
            SERVERPORT = Integer.parseInt(args[1]);
        }
        GServer server = new GServer();
        System.out.println("Server initialises a socket open for requests");
        try {
            server.start(SERVERPORT);
            System.out.println("Server proceeds to terminate its socket");
            server.stop();
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    private class ClientHandler extends Thread {
        private Socket clientSocket;
        private PrintWriter out;
        private BufferedReader in;

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
            System.out.println("Denied client, server's closing their socket and associated thread");
            out.println("Denying connection: too many clients currently.");
            try {
                in.close();
                out.close();
                clientSocket.close();
            } catch (IOException e) {
                System.err.println("Closing; client socket thread failed: \"" + clientSocket.toString() + "\"");
            }
        }

        public void run() { // client-specific thread running on server
            try {
                System.out.println("In thread; socket thread on server started: \"" + clientSocket.toString() + "\"");
                String inputLine = in.readLine();
                while (inputLine != null) {
                    System.out.println("In thread; server responds to client \"" + clientSocket.toString() + "\": " + inputLine);
                    try {
                        JSONObject json = new JSONObject(inputLine);
                        // Distinguish between Server/Job related commands
                        Integer job = json.optInt("Job");
                        if (job != null) { // this command has job number, so needs to be handled by server-aware code
                             String command = json.optString("Command");
                             if (command.equals("Retrieve")) {
                                 GRequest request = queue.getRequest(job);
                                 if (request != null)
                                     out.println(request.toJSON());
                                 else
                                     out.println("No such job in queue: " + job);
                             } else if (command.equals("Output")) {

                             }
                        } else { // this command is a new, actual job so needs to be managed and may be queued
                            GRequest greq = GRequest.fromJSON(json);
                            out.println("OK--command detected: " + greq);
                            // decision: run it now, or queue it for later...
                            if (greq.isQueued()) { // queue the request
                                queue.add(greq);
                                System.out.println("Job " + greq.getJob() + " has been dispatched to queue, cancel with {\"Job\":\" + greq.getJob() + \",\"Command\":\"Cancel\"}");
                            } else { // run now ... not in separate thread
                                System.out.println("Once complete use {\"Job\":" + greq.getJob() + ",\"Command\":\"Output\"}");
                                greq.run();
                            }
                        }
                        System.out.println("Accepts new request");
                    } catch (JSONException e) {
                        out.println("OK--no JSON detected. Waiting for another command");
                    }
                    inputLine = in.readLine();
                }
                System.out.println("Finished with client, server's closing their socket and associated thread");
                in.close();
                out.close();
                clientSocket.close();
                clients.remove(this);
            } catch (IOException e) {
                System.err.println("In thread; client socket thread failed: \"" + clientSocket.toString() + "\"");
            }
        }
    }
}
