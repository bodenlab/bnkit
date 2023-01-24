package asr;

import java.net.*;
import java.io.*;

public class GClient extends Thread {

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
        System.out.println("asr.GClient [-s <server-name>] [-p <port>] [-h]");
        System.out.println("\twhere <server-name> is the URL for the server (default is 127.0.0.1 or localhost)");
        System.out.println("\t<port> is the port number that the server is using (default is 4072)");
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
        String SERVERNAME = "127.0.0.1"; // localhost

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("p") && args.length > a + 1) {
                    try {
                        SERVERPORT = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(1, args[a] + "port must be an integer 0 - 65535 (default is 4072)");
                    }
                } else if (arg.equalsIgnoreCase("s") && args.length > a + 1) {
                    SERVERNAME = args[++a];
                } else if (arg.equalsIgnoreCase("h")) {
                    usage();
                }
            }
        }
        GClient client = new GClient(SERVERNAME, SERVERPORT);
        System.out.println("Client initialises a socket open for requests");
        client.start();
        System.out.println("Input stream (user typing) and output stream (server responding)\nare read and written asynchronously. ");
    }

    // initialize socket and input output streams
    private Socket socket            = null;
    private ClientServerOutputReader csor = null;
    private ClientUserInputReader cuir = null;

    private InputStreamReader stdin = null;
    private InputStream in_stream = null;
    private PrintWriter out     = null;

    // constructor to put ip address and port
    public GClient(String SERVERNAME, int SERVERPORT) {
        try {
            // establish a connection
            socket = new Socket(SERVERNAME, SERVERPORT);
            csor = new ClientServerOutputReader(socket);
            cuir = new ClientUserInputReader(socket);
        } catch (UnknownHostException e) {
            System.err.println("Invalid host " + SERVERNAME);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("Couldn't get I/O for the connection to " + SERVERNAME);
            System.exit(1);
        }
    }

    @Override
    public void run() {
        csor.start();
        cuir.start();
    }
}

class ClientServerOutputReader extends Thread {
    Socket serverSocket;
    public ClientServerOutputReader(Socket serverSocket){
        this.serverSocket = serverSocket;
    }

    public void run() {

        try {
            BufferedReader in = new BufferedReader(new InputStreamReader(serverSocket.getInputStream()));
            String outputFromServer="";
            while((outputFromServer=in.readLine())!= null){
                // This output may overlap the user input from console
                System.out.println(outputFromServer);
            }
        } catch (SocketException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

class ClientUserInputReader extends Thread {
    Socket serverSocket;
    public ClientUserInputReader(Socket serverSocket){
        this.serverSocket = serverSocket;
    }
    public void run(){
        BufferedReader stdIn = new BufferedReader(new InputStreamReader(System.in));
        PrintWriter out;
        int row = 1;
        try {
            out = new PrintWriter(serverSocket.getOutputStream(), true);
            // System.out.print("[" + row + "]>");
            String userInput;
            while ((userInput = stdIn.readLine()) != null) {
                out.println(userInput);
                row += 1;
                // System.out.print("[" + row + "]>");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
