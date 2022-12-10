package asr;

import java.net.*;
import java.io.*;

public class GClient
{
    // initialize socket and input output streams
    private Socket socket            = null;
    private DataInputStream  input   = null;
    private PrintWriter out     = null;

    // constructor to put ip address and port
    public GClient(String address, int port) {
        // establish a connection
        try {
            socket = new Socket(address, port);
            System.out.println("Connected");
            // takes input from terminal
            input  = new DataInputStream(System.in);
            // sends output to the socket
            out = new PrintWriter(socket.getOutputStream(), true);
        }
        catch(UnknownHostException u) {
            System.out.println(u);
        } catch(IOException i) {
            System.out.println(i);
        }

        // string to read message from input
        String line = "";
        // keep reading until "Over" is input
        while (!line.equals("Over")) {
            try {
                line = input.readLine();
                System.out.println("Trying to send: " + line.toUpperCase());
                out.println(line + "\n");
                out.flush();
            } catch(IOException e) {
                System.err.println("Client could not read and write: " + e);
            }
        }
        // close the connection
        try {
            input.close();
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