import java.io.*;
import java.net.*;
import java.util.StringTokenizer;

// to test server alone, start command prompt and type:
// %  telnet localhost 4444
// send it requests as if it were clientRequest
// alternatively, use MatClient.java and testClient.java 
public class MatServer 
{
  MatlabControl mc = new MatlabControl(); 
  public void body() throws IOException {
    class Caller extends Thread {
      public void run()
      {
        try 
        {
          body2();
        }
        catch( IOException e ) { System.err.println( e ); }
      }
      public void body2()  throws IOException
      {
        ServerSocket serverSocket = null;
        boolean listening = true;
        try
        {
          serverSocket = new ServerSocket( 4444 ); 
        }
        catch( IOException e ) { System.err.println( e ); System.exit(1); }
    
        while( listening )
        {
          new workerThread( serverSocket.accept() ).run();
        }
        serverSocket.close();
      } // end_run
    } // end_caller
    Caller c = new Caller();
    c.start();
  }

  public class workerThread 
  {
    private Socket socket = null;
    public workerThread( Socket socket ) { this.socket = socket; }
  
    public void run()
    {
      try
      {
        PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
        BufferedReader in = new BufferedReader(
                                      new InputStreamReader(
                                      socket.getInputStream()));
        String inputLine, outputLine;
        while(( inputLine = in.readLine()) != null) 
        {
          outputLine = processInput( inputLine );
          out.println( outputLine );
          if( outputLine.equals( "bye" ))
            break;
        }
  
        out.close();
        in.close();
        socket.close();
      }
      catch( IOException e ) { System.err.println( e ); }
    }
  
    // the server currently handles 4 requests:
    // 1. help -- checks status of connection
    // 2. bye -- disconnects with client
    // 3. testEE -- matlab script with specific syntax
    // 4. testDM -- matlab script with specific syntax
    public String processInput( String req ) 
    {
      String rez = "";
      StringTokenizer st;
      if( req.startsWith( "help" ) ) { rez = "ok"; }
      else if( req.startsWith( "bye" ) ) { rez = "bye"; }
      else if( req.startsWith( "testEE" ) )
      {
        // parse request: "testEE Q B S O iter"
        st = new StringTokenizer( req, " " );
        String scriptname = "";
        if( st.hasMoreTokens() ) { scriptname = st.nextToken(); }
        Integer[] params = new Integer[5];
        String c;
        int i = 0;
        while( st.hasMoreTokens() )
        {
          c = st.nextToken();
          int n = Integer.parseInt( c );
          params[i] = new Integer( n );
          i++;
        }
        if( i < 5 ) { rez = "not sent"; }
        else
        {
          // send it to matlab
          Object[] args = new Object[5];
          for( i=0; i<5; i++ ) args[i] = params[i];
          if( scriptname.equals( "testEE" ) ) 
          {
            try
            {
              Object[] returnVals = new Object[1];
              returnVals[0] = mc.blockingFeval( scriptname, args ); 
            }
            catch( Exception e ) { System.err.println( e ); }
          }
          rez = "1";
        }
      }// end_if(testEE)
      else if( req.startsWith( "testDM" ) )
      {
        // parse request: "testDM Q iter"
        st = new StringTokenizer( req, " " );
        String scriptname = "";
        if( st.hasMoreTokens() ) { scriptname = st.nextToken(); }
        Integer[] params = new Integer[2];
        String c;
        int i = 0;
        while( st.hasMoreTokens() )
        {
          c = st.nextToken();
          int n = Integer.parseInt( c );
          params[i] = new Integer( n );
          i++;
        }
        if( i < 2 ) { rez = "not sent"; }
        else
        {
          // send it to matlab
          Object[] args = new Object[2];
          for( i=0; i<2; i++ ) args[i] = params[i];
          Object[] returnVals = new Object[1];
          if( scriptname.equals( "testDM" ) ) 
          {
            try
            {
              returnVals[0] = mc.blockingFeval( scriptname, args ); 
              //mc.feval(new String("disp"), returnVals);
            }
            catch( Exception e ) { System.err.println( e ); }
          }
          // parse output
          String filename = new String("cp.txt");
          BufferedReader br = null;
          try
          {
            br = new BufferedReader(new FileReader( filename )); 
            String text;
            while(( text = br.readLine() ) != null )
            {
              st = new StringTokenizer( text );
              while( st.hasMoreTokens() )
              {
                c = st.nextToken();
                rez = c;
              }
            }
          }
          catch( IOException ioe ) { System.err.println( ioe ); }
        }
      }// end_if(testDM)
      else rez = "nothing";
      return rez;
    }// end_processInput
  }// end_workerThread_class
}// end_MatServer_class
