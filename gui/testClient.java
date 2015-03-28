import java.io.*;
import java.net.*;

public class testClient
{
  public static void main( String args[] )
  {
    InputStreamReader isr = new InputStreamReader ( System.in );
    StreamTokenizer   st  = new StreamTokenizer ( isr );
    MatClient c = new MatClient();
    /*
    try
    {
      c.createJob( "help" );
    }
    catch( Exception e ) { System.err.println( e ); }
    */

int iter = 2;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
    //getUserinput( st );
iter++;

    double pop = 0;
    String cmd = "";
    try 
    { 
      cmd = "testDM 1 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    //getUserinput( st );
    try
    {
      if( pop == 1 ) c.createJob( "testEE 1 2 3 11 "+iter );
      else           c.createJob( "testEE 1 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
/*
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
*/

    /* a bunch of dm
double pop = 0;
String cmd = "";
int iter = 2;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 2 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    */

    /* a bunch of ee
int iter = 2;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try
    {
      c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    */

    /*
    try
    {
      if( pop == 1 ) c.createJob( "testEE 2 2 3 11 "+iter );
      else           c.createJob( "testEE 2 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 3 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      if( pop == 1 ) c.createJob( "testEE 3 2 3 11 "+iter );
      else           c.createJob( "testEE 3 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 4 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      if( pop == 1 ) c.createJob( "testEE 4 2 3 11 "+iter );
      else           c.createJob( "testEE 4 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    try 
    { 
      cmd = "testDM 5 "+iter;
      pop = c.createJob( cmd );
      System.out.println( "TC "+iter+": dm " + pop );
    }
    catch( Exception e ) { System.err.println( e ); }
    try
    {
      if( pop == 1 ) c.createJob( "testEE 5 2 3 11 "+iter );
      else           c.createJob( "testEE 5 1 3 11 "+iter );
      System.out.println( "TC "+iter+": ee" );
    }
    catch( Exception e ) { System.err.println( e ); }
iter++;
    */

    try
    {
      c.createJob( "bye." );
    }
    catch( Exception e ) { System.err.println( e ); }
  }

  public static void getUserinput( StreamTokenizer s )
  {
    int c;
    try
    {
      c = s.nextToken();
    }
    catch ( IOException ioe ) { System.err.println("error in userinput."); }
  }
}
