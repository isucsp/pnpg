
import java.net.*;
import java.io.*;

public class GuiServer{
    public static void main(String[] args) throws IOException {
        ServerSocket serverSocket = null;

        int port;
        while(true){
            try {
                port=Math.random()*(49152-1024)+1024;
                serverSocket = new ServerSocket(port);
                System.out.println("Bind to port: " + port );
                break;
            } catch (IOException e) {
                System.err.println("Could not listen on port: 4444.");
                System.exit(-1);
            }
        }

        while (true)
	    new KKMultiServerThread(serverSocket.accept()).start();

        serverSocket.close();
    }
}
