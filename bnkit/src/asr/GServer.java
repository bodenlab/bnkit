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
            server.start(4072);
            //server.stop();
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
                System.out.println("Client socket thread started: \"" + clientSocket.toString() + "\"");
                out = new PrintWriter(clientSocket.getOutputStream(), true);
                in = new BufferedReader(
                        new InputStreamReader(clientSocket.getInputStream()));
                String inputLine;
                while ((inputLine = in.readLine()) != null) {
                    if (".".equals(inputLine)) {
                        out.println("bye");
                        break;
                    }
                    out.println("Server responds to client \"" + clientSocket.toString() + "\": " + inputLine.toUpperCase());
                }
                in.close();
                out.close();
                clientSocket.close();
            } catch (IOException e) {
                System.err.println("Client socket thread failed: \"" + clientSocket.toString() + "\"");
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