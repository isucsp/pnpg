import java.net.*;
import java.io.*;
import java.util.regex.*;
import java.util.*;

public class ReadNist {
    String html;

    public static void main(String[] args) throws Exception {
        ReadNist reader=new ReadNist();
        String url;
        html = ParseHTML.readHTML("http://www.nist.gov/pml/data/xraycoef/");

        url = ParseHTML.getLink(html,"\\s*Table 2.\\s*");
        ArrayList<Material> mcTable=reader.readTable2(url);
        Iterator<Material> itr=mcTable.iterator();
        /*
        System.out.print("density=[\n");
        while(itr.hasNext()){
            System.out.println(itr.next().density+"; ");
        }
        System.out.println("];\n");

        url = ParseHTML.getLink(html,"\\s*Table 1.\\s*");
        ArrayList<Material> mcTable=reader.readTable1(url);
        Iterator<Material> itr=mcTable.iterator();
        System.out.print("density=[\n");
        while(itr.hasNext()){
            System.out.println(itr.group().density+"; ");
        }
        System.out.println("];\n");

        Iterator<double[]> itr3;
        double[] tmp;
        for(int i=1; i<1+mcTable.size(); i++){
            System.out.printf("mac{%d} = [\n",i);
            itr3=reader.readTable3(i).iterator();
            while(itr3.hasNext()){
                tmp=itr3.group();
                System.out.printf("%9e, %8e, %8e;\n", tmp[0],tmp[1],tmp[2]);
            }
            System.out.printf("];\n\n");
        }
        */
    }

    ArrayList<Material> readTable1(String url) throws Exception{
        String in=ParseHTML.readHTML(url);

        ParseHTML table=new ParseHTML(in,"table");
        if(! table.hasNext() )
            return null;
        table.refine("tr");

        Material singlerow;
        int col;
        //System.out.println(content);
        ArrayList<Material> mcTable=new ArrayList<Material>();
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
                if(cell.matches(Material.heads[col])){
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
                            singlerow=new Material(id,sym,name,za,i,density);
                            mcTable.add(singlerow);
                            //System.out.println(singlerow);
                            break;
                    }
                    col++;
                }
            }
        }
        return mcTable;
    }

    ArrayList<Material> readTable24(String url, String codeGen){
        String table4;

        if(codeGen!=null && codeGen.equals("")) return;
        else
            table4=ParseHTML.readHTML(ParseHTML.getLink(html,"\\s*Table 4.\\s*"));

        String in=ParseHTML.readHTML(url);

        ParseHTML table=new ParseHTML(in,"table");
        if(! table.hasNext() )
            return null;
        table.refine("tr");

        Composition singlerow;
        int col,row=0;
        //System.out.println(content);
        ArrayList<Material> mcTable=new ArrayList<Material>();
        int id=0;
        double za=0,i=0,density=0;
        String sym="", name="", cell;

        ParseHTML rowHTML;
        Iterator<String> itr;

        while(table.hasNext()){
            row++;
            //System.out.println(table.group());
            if(row<=2) continue;
            rowHTML=table.copy();
            rowHTML.refine("t(d|h)");
            itr=rowHTML.findAll().iterator();

            col=0;
            while(itr.hasNext()){
                cell=itr.next().trim();
                if(cell.matches(Composition.heads[col])){
                    //System.out.println(cell);
                    switch(col){
                        case 0: 
                            name=cell;
                            //System.out.println(name);
                            break;
                        case 1: 
                            za=Double.parseDouble(cell);
                            break;
                        case 2: 
                            i=Double.parseDouble(cell);
                            break;
                        case 3: 
                            density=Double.parseDouble(cell);
                            break;
                        case 4: 
                            sym=findComUrl(table4,name,codeGen);
                            sym=name.replaceFirst("^[^A-Z]*","").split("[\\(\\),/\\s]+")[0].replaceAll("-","");
                            singlerow=new Composition(id,sym,name,za,i,density);
                            id++;
                            Matcher components = Pattern.compile("(\\d{1,2})\\s*:\\s*([0-9.]+)").matcher(cell);
                            while(components.find())
                                singlerow.add(Integer.parseInt(components.group(1)),
                                        Float.parseFloat(components.group(2)));
                            System.out.println(singlerow);
                            mcTable.add(singlerow);
                            break;
                    }
                    col++;
                }
            }
        }
        return mcTable;
    }

    String findComUrl(String name, String codeGen){
        if(codeGen!=null && codeGen.equals("")) return;
        name = name.toLowerCase();

    }

    ArrayList<Material> readTable4(String url){
        String in=ParseHTML.readHTML(url);


    }


    ArrayList<double[]> readAttenTable(url) throws Exception{
        String in=ParseHTML.readHTML(url);

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


}

class Material{
    static final String[] heads={"\\d{1,2}", "[A-Z][a-z]?", "(\\w|,|\\s)+", 
        "[0-9.+\\-Ee]+", "[0-9.+\\-Ee]+", "[0-9.+\\-Ee]+"};
    int id;
    double za,i,density;
    String sym,name;
    Material(int id,String sym, String name, double za,
            double i, double density){
        this.id=id;
        this.sym=sym;
        this.name=name;
        this.za=za;
        this.i=i;
        this.density=density;
    }
    public String toString(){
        return String.format("%2d %4s %20s %9.5f %7.1f %18e",
            id,sym, name, za, i, density);
    }
}

class Composition extends Material{
    static final String[] heads={"[a-zA-Z\\-0-9\\(\\),/\\s]+", 
        "[0-9.+\\-Ee]+", "[0-9.+\\-Ee]+", "[0-9.+\\-Ee]+", "((\\d{1,2})\\s*:\\s*([0-9.]+)[^0-9]*)+"};
    ArrayList<Integer> comp;
    ArrayList<Float> prop;
    Composition(int id, String sym, String name, double za,
            double i, double density){
        super(id,sym,name,za,i,density);
        comp = new ArrayList<Integer>();
        prop = new ArrayList<Float>();
    }
    public void add(int c, float p){
        comp.add(c); prop.add(p);
    }
    public String toString(){
        String str= String.format("%2d %10s %30s %9.5f %7.1f %18e",
            id,sym, name, za, i, density) + "\n";
        for(int i=0; i<comp.size(); i++)
            str += "\t\t\t" + comp.get(i) +" : " + prop.get(i) + "\n";
        return str;
    }
}

