package asr;

import java.net.*;
import java.io.*;

public class GClient
{
    // initialize socket and input output streams
    private Socket socket            = null;
    private InputStreamReader stdin = null;
    private InputStream in_stream = null;
    private PrintWriter out     = null;

    // constructor to put ip address and port
    public GClient(String address, int port) {
        // establish a connection
        try {
            socket = new Socket(address, port);
            System.out.println("Client connected to server");
            // stream from terminal
            stdin = new InputStreamReader(System.in);
            in_stream = socket.getInputStream();
            // sends output to the socket
            out = new PrintWriter(socket.getOutputStream(), true);

        }
        catch(UnknownHostException u) {
            System.out.println("Problem host: " + u);

        } catch(IOException i) {
            System.out.println("Problem I/O:" + i);
        }

        BufferedReader reader = new BufferedReader(stdin);
        BufferedReader server_input = new BufferedReader(new InputStreamReader(in_stream));

        // string to read message from input
        String line = "";
        // keep reading until "Over" is input
        while (!line.startsWith("Exit")) {
            try {
                line = reader.readLine();
                System.out.println("Trying to send: " + line.toUpperCase());
                out.println(line);
                String server_line = server_input.readLine();
                System.out.println(server_line);
                if (server_line.equals("Breaking-up"))
                    break;
            } catch(IOException e) {
                System.err.println("Client could not read and write: " + e);
            }
        }
        // close the connection
        try {
            stdin.close();
            out.close();
            socket.close();
        } catch(IOException i) {
            System.out.println(i);
        }
    }

    public static void main(String args[]) {
        GClient client = new GClient("127.0.0.1", 4072);

    }
}