package asr;

import java.net.*;
import java.io.*;

public class GServer {

    private ServerSocket serverSocket;

    public void start(int port) throws IOException {
        serverSocket = new ServerSocket(port);
        while (true) {
            Socket client = serverSocket.accept();
            ClientHandler handler = new ClientHandler(client);
            handler.start(); // spawn-off client thread
        }

    }

    public void stop() throws IOException {
        serverSocket.close();
    }

    // constructor with port
    public GServer() {
    }

    public static void main(String args[]) {
        GServer server = new GServer();
        try {
            System.out.println("Server initialises a socket open for requests");
            server.start(4072);
            System.out.println("Server proceeds to terminate its socket");
            server.stop();
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    private static class ClientHandler extends Thread {
        private Socket clientSocket;
        private PrintWriter out;
        private BufferedReader in;

        public ClientHandler(Socket socket) {
            this.clientSocket = socket;
            System.out.println("Client socket created: \"" + clientSocket.toString() + "\"");
        }

        public void run() { // client thread
            try {
                System.out.println("In thread; socket thread on server started: \"" + clientSocket.toString() + "\"");
                OutputStream out_stream = clientSocket.getOutputStream();
                InputStream in_stream = clientSocket.getInputStream();
                out = new PrintWriter(out_stream, true);
                in = new BufferedReader(new InputStreamReader(in_stream));
                String inputLine = in.readLine();
                while (inputLine != null) {
                    if (inputLine.startsWith("Exit")) {
                        out.println("Confirming that you are wilfully exiting the server.");
                        break;
                    }
                    System.out.println("In thread; server responds to client \"" + clientSocket.toString() + "\": " + inputLine.toUpperCase());
                    if (inputLine.contains("break"))
                        out.println("Server is breaking-up with client");
                    else
                        out.println("OK. Waiting for another command");
                    inputLine = in.readLine();
                }
                System.out.println("Finished with client, server's closing their socket and associated thread");
                in.close();
                out.close();
                clientSocket.close();
            } catch (IOException e) {
                System.err.println("In thread; client socket thread failed: \"" + clientSocket.toString() + "\"");
            }
        }
    }
}
/*
public class EchoMultiServer {
    private ServerSocket serverSocket;

    public void start(int port) {
        serverSocket = new ServerSocket(port);
        while (true)
            new EchoClientHandler(serverSocket.accept()).start();
    }

    public void stop() {
        serverSocket.close();
    }

    private static class EchoClientHandler extends Thread {
        private Socket clientSocket;
        private PrintWriter out;
        private BufferedReader in;

        public EchoClientHandler(Socket socket) {
            this.clientSocket = socket;
        }

        public void run() {
            out = new PrintWriter(clientSocket.getOutputStream(), true);
            in = new BufferedReader(
              new InputStreamReader(clientSocket.getInputStream()));

            String inputLine;
            while ((inputLine = in.readLine()) != null) {
                if (".".equals(inputLine)) {
                    out.println("bye");
                    break;
                }
                out.println(inputLine);
            }

            in.close();
            out.close();
            clientSocket.close();
    }
}
 */