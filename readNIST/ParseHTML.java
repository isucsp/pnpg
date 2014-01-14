
import java.util.regex.*;
import java.util.*;
public class ParseHTML{
    private static String content;
    private String tag;
    private Matcher matcher;
    private boolean recursive;
    int[] innerRegion = new int[2], outerRegion= new int[2];
    int count=0;

    ParseHTML(String content){
        this(content,"html");
    }

    ParseHTML(String content, String tag){
        this.content=content;
        this.tag=tag;
        String regex="(?ims)</?" + tag + "((\\s[^>]*)|([^>]*))>";
        matcher=Pattern.compile(regex).matcher(content);
        //System.out.println(pattern.flags()+" "
        //        + (Pattern.CASE_INSENSITIVE | Pattern.DOTALL));
        recursive=setRec(tag);
        count=0;
    }

    void setTag(String tag){
        this.tag=tag;
        String regex="(?ims)</?" + tag + "((\\s[^>]*)|([^>]*))>";
        matcher.usePattern(Pattern.compile(regex));
        recursive=setRec(tag);
        count=0;
    }

    void setRegion(int s, int e){ 
        matcher.region(s,e);
        count=0;
    }

    void refine(){
        matcher.region(innerRegion[0],innerRegion[1]);
        count=0;
    }

    void refine(String tag){
        setTag(tag);
        refine();
    }

    ParseHTML copy(){
        ParseHTML res=new ParseHTML(content,tag);
        res.innerRegion=innerRegion; res.outerRegion=outerRegion;
        return res;
    }

    String next(){
        return content.substring(innerRegion[0],innerRegion[1]);
    }

    void reset(){
        matcher.region(matcher.regionStart(), matcher.regionEnd());
        count=0;
    }

    int find(int ith){
        reset();
        while(matcher.find()){
            count++;
            if(count==ith) break;
        }
        return count;
    }

    boolean hasNext(){
        if(matcher.find()){
            int[] openTag={matcher.start(), matcher.end()};
            int[] closeTag=new int[2];
            int depth=0, innerCount=0;
            innerRegion[0]=matcher.end();
            outerRegion[0]=matcher.start();
            if(matcher.group().charAt(1)=='/'){
                System.err.printf("L:%d, no opening tag for %s!\n",
                        getLineNumber(openTag[0]), content.substring(openTag[0],openTag[1]));
                System.err.print("depth: " + depth);
                System.err.println("; innerCount: " + innerCount);
                System.err.println("last found: " + content.substring(openTag[0],openTag[1]));
                //System.err.println("Search from: " +
                //        content.substring(matcher.regionStart(),matcher.regionEnd()));
                //reset
                return false;
            }
            depth++; innerCount++; count++;
            while(depth>0){
                //System.out.println("depth=" + depth);
                if(matcher.find()){
                    count++; innerCount++;
                    closeTag[0]=matcher.start(); closeTag[1]=matcher.end();
                    if(matcher.group().charAt(1)=='/')
                        depth--;
                    else{
                        if(recursive)
                            depth++;
                        else{
                            innerRegion[1]=matcher.start();
                            outerRegion[1]=matcher.start();
                            find(count-1);
                            System.err.println("alone tag: "+matcher.group() +
                                    " at L:" + getLineNumber(matcher.start()));
                            //System.err.println("Search from: " +
                            //        content.substring(matcher.regionStart(),matcher.regionEnd()));
                            return true;
                        }
                    }
                }else{
                    System.err.printf("L:%d, no closing tag for %s!\n",
                            getLineNumber(openTag[0]), content.substring(openTag[0],openTag[1]));
                    System.err.print("depth: " + depth);
                    System.err.println("; count: " + count);
                    //System.err.printf("last found: " + content.substring(lastFind[0],lastFind[1]));
                    //System.err.println("Search from: " +
                    //        content.substring(matcher.regionStart(),matcher.regionEnd()));
                    reset();
                    return false;
                }
            }
            innerRegion[1]=matcher.start();
            outerRegion[1]=matcher.end();
            return true;
        }else
            reset();
        return false;
    }

    ArrayList<String> findAll(){
        return findAll(0);
    }

    ArrayList<String> findAll(int limit){
        ArrayList<String> output = new ArrayList<String>();
        int depth=0, count=0, lastFind=0;
        while(hasNext()){
            output.add(next()); count++;
            if(limit>0 && count>=limit) break;
        }
        return output;
    }

    int getLineNumber(int pos){
        Matcher m=Pattern.compile("(?m)$").matcher(content);
        m.region(0,pos);
        int count=0;
        while(m.find()){
            //System.out.println("Find \"" + m.group() + "\" at " + m.start());
            count++;
        }
        return count;
    }

    boolean setRec(String tag){
        tag=tag.toLowerCase();
        if(tag.equals("html") || tag.equals("table"))
            return true;
        return false;
    }
}
