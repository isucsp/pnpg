import java.net.*;
import java.io.*;
import java.util.regex.*;
import java.util.*;

public class ReadNist {
    public static void main(String[] args) throws Exception {
        ReadNist reader=new ReadNist();
        ArrayList<MaterialConst> mcTable=reader.readTable1();
        Iterator<MaterialConst> itr=mcTable.iterator();
        System.out.print("density=[\n");
        while(itr.hasNext()){
            System.out.println(itr.next().density+"; ");
        }
        System.out.println("];\n");
        Iterator<double[]> itr3;
        double[] tmp;
        for(int i=1; i<1+mcTable.size(); i++){
            System.out.printf("mac{%d} = [\n",i);
            itr3=reader.readTable3(i).iterator();
            while(itr3.hasNext()){
                tmp=itr3.next();
                System.out.printf("%9e, %8e, %8e;\n", tmp[0],tmp[1],tmp[2]);
            }
            System.out.printf("];\n\n");
        }
    }

    ArrayList<double[]> readTable3(int id) throws Exception{
        URL tab1 = new URL(
                String.format("http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z%02d.html",id));
        BufferedReader inRead = new BufferedReader(
                new InputStreamReader(tab1.openStream()));
        //BufferedReader inRead = new BufferedReader(
        //        new InputStreamReader(new FileInputStream("table1.html")));

        String inputLine;
        StringBuilder content=new StringBuilder();
        while ((inputLine = inRead.readLine()) != null)
            content.append(inputLine+"\n");
        inRead.close();
        //System.out.println(content);
        //System.out.println(parseMark("table",content.toString(),1)[0]);
        String in=content.toString();

        ParseHTML table=new ParseHTML(in, "table");
        if(table.hasNext())
            table.refine();
        if(table.hasNext())
            table.refine("tr");


        ParseHTML rowHTML;
        Iterator<String> itr;

        ArrayList<double[]> macTable=new ArrayList<double[]>();
        int c; String cell;
        double[] data;
        while(table.hasNext()){
            rowHTML=table.copy();
            rowHTML.refine("t(d|h)");
            itr=rowHTML.findAll().iterator();

            c=0;
            data=new double[3];
            while(itr.hasNext()){
                cell=itr.next().trim();
                try{
                    data[c]=Double.parseDouble(cell);
                }catch(Exception e){
                    continue;
                }
                c++;
            }
            if(c==3){
                macTable.add(data);
                //System.out.println(Arrays.toString(data));
            }
        }
        return macTable;
    }

    ArrayList<MaterialConst> readTable1() throws Exception{
        URL tab1 = new URL("http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html");
        BufferedReader inRead = new BufferedReader(
                new InputStreamReader(tab1.openStream()));
        //BufferedReader inRead = new BufferedReader(
        //        new InputStreamReader(new FileInputStream("table1.html")));

        String inputLine;
        StringBuilder content=new StringBuilder();
        while ((inputLine = inRead.readLine()) != null)
            content.append(inputLine+"\n");
        inRead.close();
        //System.out.println(content);
        //System.out.println(parseMark("table",content.toString(),1)[0]);
        String in=content.toString();

        ParseHTML table=new ParseHTML(in,"table");
        if(! table.hasNext() )
            return null;
        table.refine("tr");

        MaterialConst singlerow;
        int col;
        String[] heads={"\\d{1,2}", "[A-Z][a-z]?", "(\\w|,|\\s)+", 
            "[0-9.+\\-E]+", "[0-9.+\\-E]+", "[0-9.+\\-E]+"};
        //System.out.println(content);
        ArrayList<MaterialConst> mcTable=new ArrayList<MaterialConst>();
        int id=0;
        double za=0,i=0,density=0;
        String sym="", name="", cell;

        ParseHTML rowHTML;
        Iterator<String> itr;

        while(table.hasNext()){
            rowHTML=table.copy();
            rowHTML.refine("t(d|h)");
            itr=rowHTML.findAll().iterator();

            col=0;
            while(itr.hasNext()){
                cell=itr.next().trim();
                if(cell.matches(heads[col])){
                    switch(col){
                        case 0: 
                            id=Integer.parseInt(cell);
                            break;
                        case 1: 
                            sym=cell;
                            break;
                        case 2: 
                            name=cell;
                            //System.out.println(name);
                            break;
                        case 3: 
                            za=Double.parseDouble(cell);
                            break;
                        case 4: 
                            i=Double.parseDouble(cell);
                            break;
                        case 5: 
                            density=Double.parseDouble(cell);
                            break;
                    }
                    col++;
                }
            }
            if(col==heads.length){
                singlerow=new MaterialConst(id,sym,name,za,i,density);
                mcTable.add(singlerow);
                //System.out.println(singlerow);
            }
        }
        return mcTable;
    }

}

class MaterialConst{
    int id;
    double za,i,density;
    String sym,name;
    MaterialConst(int id,String sym, String name, double za,
            double i, double density){
        this.id=id;
        this.sym=sym;
        this.name=name;
        this.za=za;
        this.i=i;
        this.density=density;
    }
    public String toString(){
        return String.format("%2d%4s%20s%9.5f%7.1f%18e",
            id,sym, name, za, i, density);
    }
}

