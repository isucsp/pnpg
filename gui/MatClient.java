
import java.io.*;
import java.net.*;
import javax.swing.*;

public class MatClient
{
    Socket sock;
    PrintWriter out;
    BufferedReader in;
    DataInputStream din;
    String sendRequest, fromServer;
    JLabel status;
    int port;
    String host;
    public MatClient(String host, int port, JLabel status){
        sock = null;
        out  = null;
        in   = null;
        sendRequest = "";
        fromServer  = "";
        this.host=host;
        this.port=port;
        this.status=status;
    }

    public void connect(){
        try{
            if(sock==null){
                System.out.println("Establish socket to server at " + host + ":" + port);
                sock = new Socket(host, port );
                out = new PrintWriter(( sock.getOutputStream()),true);
                din = new DataInputStream(new BufferedInputStream( sock.getInputStream()));
                in = new BufferedReader( new InputStreamReader(din) );
                try{
                    System.out.println("Get from matlab: " + in.readLine());
                }catch(IOException e){ e.printStackTrace(); }
                System.out.println("Connected");
                System.out.println("isClosed()=" + sock.isClosed());
                System.out.println("isConnected()=" + sock.isConnected());
            }else{
                System.out.println("isClosed()=" + sock.isClosed());
                System.out.println("isConnected()=" + sock.isConnected());
                System.out.println("isInputShutdown()=" + sock.isInputShutdown());
                System.out.println("isOutputShutdown()=" + sock.isOutputShutdown());
            }
        }catch ( UnknownHostException e ){
            System.err.println("no localhost"); System.exit(1);
        }catch (IOException e) {
            System.err.println( e ); System.exit(1);
        }

    }

    public void disconnect(){
        if(sock!=null)
            try{
                exec("quit","");
                sock.close();
                in.close();
                out.close();
                sock=null;
            }catch(IOException e){ e.printStackTrace();}
    }

    public void exec(String cmd, String content){
        out.printf("%8s%s\n",cmd, content);
    }
    public void exec(String content){
        out.printf("%8s%s\n","cmd", content);
    }
    public String retString(){
        String ret=null;
        try{
            ret=in.readLine();
            /*while(false){
                if(in.ready())
                    ret=in.readLine();
                else
                    System.out.println("wait for newline");
            }*/
            //else Thread.sleep(10);
        }catch(IOException e){ e.printStackTrace(); }
        //catch(InterruptedException e){ e.printStackTrace(); }
        return ret;
    }

    public int retArray(double[] img, int cnt){
        int i=0;
        try{
            while(i<cnt){
                img[i++]=din.readDouble();
                //if(din.available()!=0)
                //    img[i++]=din.readDouble();
                //else
                //    System.out.println("wait for double");
                //System.out.println("i="+i);
            }
            System.out.println("end receiving for double");
        }catch(IOException e){ e.printStackTrace(); }
        //catch(InterruptedException e){ e.printStackTrace(); }
        return i;
    }

    public double readDouble(){
        double re=0;
        try{
            re = din.readDouble();
        }catch(IOException e){ e.printStackTrace(); }
        return re;
    }

    public void finishJob() throws IOException
    {
        exec("fclose(server); mServer; break; mServer");
        in.close();
        out.close();
        sock.close();
    }

    public double createJob( String j ) throws IOException
    {
        sendRequest = j;
        // uncomment for debugging purposes
        // System.out.println("Send Matlab request" );
        out.println( sendRequest ); 
        double rez = 0.0;
        while(( fromServer = in.readLine() ) != null ) 
        {
            // comment this out when not debugging 
            System.out.println( "server: " + fromServer );
            if( fromServer.startsWith( "bye" ) )
            {
                finishJob();
                break;
            }
            else
            {
                // process answer from server
                if( sendRequest.startsWith( "testEE" ) )
                    rez = Double.parseDouble( fromServer );
                if( sendRequest.startsWith( "testDM" ) )
                    rez = Double.parseDouble( fromServer );

                // uncomment for debugging purposes
                // System.out.println("Result Returned");
                break;
            }
        }
        return rez;
    }
}
