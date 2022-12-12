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
            // stream from server
            in_stream = socket.getInputStream();
            // stream from terminal
            stdin = new InputStreamReader(System.in);
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

        String SERVERNAME = "127.0.0.1"; // localhost
        Integer SERVERPORT = 4072; // default UQ's post-code

        if (args.length == 2) {
            SERVERNAME = args[0];
            SERVERPORT = Integer.parseInt(args[1]);
        }

        try {
            Socket serverSocket = new Socket(SERVERNAME, SERVERPORT);
            ClientServerOutputReader csor = new ClientServerOutputReader(serverSocket);
            csor.start();
            ClientUserInputReader cuir = new ClientUserInputReader(serverSocket);
            cuir.start();
        } catch (UnknownHostException e) {
            System.err.println("Don't know about host " + SERVERNAME);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("Couldn't get I/O for the connection to " +
                    SERVERNAME);
            System.exit(1);
        }

        //GClient client = new GClient("127.0.0.1", 4072);

    }
}

class ClientServerOutputReader extends Thread {
    Socket serverSocket;
    public ClientServerOutputReader(Socket serverSocket){
        this.serverSocket = serverSocket;
    }

    public void run() {
        try {
            BufferedReader in = new BufferedReader(
                    new InputStreamReader(serverSocket.getInputStream()));

            String outputFromServer="";
            while((outputFromServer=in.readLine())!= null){
                //This part is printing the output to console
                //Instead it should be appending the output to some file
                //or some swing element. Because this output may overlap
                //the user input from console
                System.out.println(outputFromServer);
            }
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
        BufferedReader stdIn = new BufferedReader(
                new InputStreamReader(System.in));
        PrintWriter out;
        try {
            out = new PrintWriter(serverSocket.getOutputStream(), true);
            String userInput;

            while ((userInput = stdIn.readLine()) != null) {
                out.println(userInput);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
