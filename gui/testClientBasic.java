import java.io.*;
import java.net.*;

public class testClientBasic
{
  public static void main( String args[] )
  {
    MatClient c = new MatClient();
    try
    {
      c.createJob( "help" );
    }
    catch( Exception e ) { System.err.println( e ); }
    // continue sending matlab commands

    try
    {
      c.createJob( "bye." );
    }
    catch( Exception e ) { System.err.println( e ); }
  }
}
