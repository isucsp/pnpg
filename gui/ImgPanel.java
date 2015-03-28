import javax.swing.*;
import javax.swing.border.*;
import java.awt.event.*;
import java.awt.color.*;
import javax.imageio.*;
import java.awt.*;
import java.io.*;
import java.awt.image.*;

public class ImgPanel extends JPanel
    implements MouseListener, MouseMotionListener, MouseWheelListener, Runnable{

    double[] raw;
    int[] iraw;
    double pMin, pMax;
    boolean rangeSet=false;
    private static final double zoomRatio = 0.2;
    int width, height, bw;
    MyImage mi;
    int tool = 3;
    Rectangle source, sink;
    DataBuffer db;
    SampleModel sm;
    ColorModel cm;
    Image cImage;
    WritableRaster wr;

    int itr=0;
    double cost=-1, rse=-1;

    public ImgPanel(MyImage mi){
        this.mi= mi;
        //mi.addParent(this);
        //setBorder(BorderFactory.createLineBorder(Color.black));
        //setBorder(new EmptyBorder(10, 10, 10, 10));
        //setBackground(Color.red);
        //raisedetched = 
        //BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
        //loweredetched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
        //raisedbevel = BorderFactory.createRaisedBevelBorder();
        //loweredbevel = BorderFactory.createLoweredBevelBorder();
        //TitledBorder border =BorderFactory.createTitledBorder("Image 
        //Reconstructed");
        //bw=(int)(border.getMinimumSize(this).getHeight()/2);
        //setBorder(border);
        addMouseListener(this);
        addMouseMotionListener(this);
        addMouseWheelListener(this);
        //setPreferredSize(new Dimension(mi.getWidth(),mi.getHeight()));
        sink=new Rectangle(0,0,getWidth(),getHeight());
        new Thread(this).start();
    }

    public void run(){
        while(true){
            if(mi.hasNewImg()){
                repaint();
            }else
                try{Thread.sleep(50);}catch(InterruptedException ee){
                    System.out.println(ee);
                }
        }
    }

    private void setCustomCursor(int i){
        Toolkit toolkit = Toolkit.getDefaultToolkit();
        String[] icons = {"res/zoom_in_16.png", "res/zoom_out_16.png",
            "res/cursor_hand.png", "res/cursor.png","res/cursor_drag_hand.png"};
        Image image = toolkit.getImage(icons[i]);
        Cursor c = toolkit.createCustomCursor(image , new Point(0,0), "img");
        //System.out.format("setCustomCursor(%d)\n",i);
        setCursor (c);
    }

    public void addStatistics(int itr, double cost, double rse){
        this.itr = itr;
        this.cost = cost;
        this.rse = rse;
    }

    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        int w= getWidth(), h= getHeight();
        System.out.format("size of imgPanel: (%d, %d)\n", w, h);
        if(w>h)
            sink.setRect(((double)w-h)/2,0, h,h);
        else
            sink.setRect(0, ((double)h-w)/2, w,w);
        Rectangle source = mi.getWindow();

        //System.out.println("from paintComponent: " + source);
        
        g.drawImage(mi.getImage(), (int)sink.getX(), (int)sink.getY(),
                (int)(sink.getX()+sink.getWidth()),
                (int)(sink.getY()+sink.getHeight()),
                (int)(source.getX()), (int)source.getY(),
                (int)(source.getX()+source.getWidth()),
                (int)(source.getY()+source.getHeight()),
                this);
        String str = "Itr: "+itr +"   cost: " + cost + "   RSE: " + rse;
        g.setColor(Color.white);
        g.drawString(str,15,15);
        /*g.drawImage(mi.getImage(), (int)sink.getX(), (int)sink.getY(),
                (int)(sink.getX()+sink.getWidth()),
                (int)(sink.getY()+sink.getHeight()),
                this);
                */
    }

    public double[] getSourcePoint(Point sinkPoint){
        return getSourcePoint(sinkPoint,mi.getWindow());
    }
    public double[] getSourcePoint(Point sinkPoint, Rectangle source){
        double ratio = source.getWidth()/sink.getWidth();
        double x = ratio*(sinkPoint.getX()-sink.getX())+source.getX();
        double y = ratio*(sinkPoint.getY()-sink.getY())+source.getY();
        if(x<0) x=0; else if(x>mi.getWidth()) x=mi.getWidth();
        if(y<0) y=0; else if(y>mi.getHeight()) y=mi.getHeight();
        return (new double[]{x, y});
    }

    public void zoomAt(double[] c, double ratio){
        double cx = c[0];
        double cy = c[1];
        Rectangle source = mi.getWindow();
        System.out.println("from zoomAt: " + source);
        System.out.println("ratio = " + ratio);
        double w = ratio*source.getWidth(),
               h = ratio*source.getHeight();
        w = Math.min(w, mi.getWidth());
        h = Math.min(h, mi.getHeight());
        double lb,rb,tb,bb;
        lb=cx-w/2; tb=cy-h/2;
        rb=lb+w; bb=tb+h;
        if(lb<0) lb=0; else if(rb>mi.getWidth()) lb=lb-(rb-mi.getWidth());
        if(tb<0) tb=0; else if(bb>mi.getHeight()) tb=tb-(bb-mi.getHeight());
        mi.setRect(lb,tb,w,h);
    }

    private void moveImage(Point orig, Point curr){
        Rectangle source = mi.getWindow();
        System.out.println(source);
        System.out.println(sink);
        double ratio = snapshot.getWidth()/sink.getWidth();
        double dx = curr.getX()-orig.getX();
        double dy = curr.getY()-orig.getY();
        System.out.format("(%f,%f)\n",dx,dy);
        System.out.format("ratio=%e\n",ratio);
        dx=-dx*ratio;
        dy=-dy*ratio;
        double lb,rb,tb,bb;
        lb = dx + snapshot.getX();
        tb = dy + snapshot.getY();
        rb = lb + snapshot.getWidth();
        bb = tb + snapshot.getHeight();
        if(lb<0) lb=0; else if(rb>mi.getWidth()) lb=lb-(rb-mi.getWidth());
        if(tb<0) tb=0; else if(bb>mi.getHeight()) tb=tb-(bb-mi.getHeight());
        mi.mvRectTo((int)lb,(int)tb);
    }

    public void mouseClicked(MouseEvent e){
        Point sinkPoint = e.getPoint();
        double[] sourcePoint = getSourcePoint(sinkPoint);
        System.out.format("mouse clicked, tool=%d\n",tool);
        switch (tool) {
            case 0: 
                if(sink.contains(sinkPoint)){
                    zoomAt(sourcePoint, 1-zoomRatio);
                    repaint();
                }
                break;
            case 1:
                if(sink.contains(sinkPoint)){
                    zoomAt(sourcePoint, 1+zoomRatio);
                    repaint();
                }
                break;
            default: break;
        }
    }
    public void mouseEntered(MouseEvent e){}
    public void mouseExited(MouseEvent e){}
    Point anchor;
    Rectangle snapshot;
    public void mousePressed(MouseEvent e){
        anchor = e.getPoint();
        snapshot = new Rectangle(mi.getWindow());
        switch (tool) {
            case 2: 
                if(sink.contains(anchor)){
                    setCustomCursor(4);
                }else{
                    anchor = null; snapshot = null;
                }
                break;
            default:
                    break;
        }
    }
    public void mouseReleased(MouseEvent e){
        anchor=null; snapshot = null;
        switch (tool) {
            case 2: 
                setCustomCursor(2);
                break;
            default:
                    break;
        }
    }
    public void mouseDragged(MouseEvent e){
        //System.out.println("dragged");
        Point sinkPoint = e.getPoint();
        switch (tool) {
            case 2: 
                if(sink.contains(anchor)){
                    moveImage(anchor, sinkPoint);
                    repaint();
                    //setCustomCursor(2);
                }else
                    ;
                break;
            default:
                    break;
        }
    }
    public void mouseMoved(MouseEvent e){
        processMouseMoving(e);
    }
    public void processMouseMoving(MouseEvent e){
        Point sinkPoint = e.getPoint();
        //System.out.format("mouse moving, tool=%d\n",tool);
        switch (tool) {
            case 2: 
                if(sink.contains(sinkPoint))
                    //setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                    setCustomCursor(2);
                else{
                    //System.out.format("set cursor default\n");
                    setCursor(Cursor.getDefaultCursor());
                }
                break;
            default:
                setCursor(Cursor.getDefaultCursor());
                break;
        }
    }
    public void mouseWheelMoved(MouseWheelEvent e){
        Point sinkPoint = e.getPoint();
        //System.out.format("mouse rotated %d\n",e.getWheelRotation());
        if(sink.contains(sinkPoint)){
            double[] sourcePoint = getSourcePoint(sinkPoint);
            zoomAt(sourcePoint, 1+zoomRatio*e.getWheelRotation());
            repaint();
        }
    }

    void update(int[] raw){
        iraw=raw;
        //mapColor();
        //mis.newPixels(0,0,width,height);
    }

    void update(){
        //mapColor();
        //mis.newPixels(0,0,width,height);
    }

    void setTool(int tool){
        this.tool=tool;
    }
       

    public static void main(String[] str){
        int w=500, h=500;
        int[] img = new int[w*h];
        JFrame frame = new JFrame("test ImgPanel");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        MyImage mi = new MyImage(img, w, h);
        ImgPanel test = new ImgPanel(mi);
        frame.getContentPane().add(test);
        //frame.add(new JLabel("frame button"));
        frame.pack();
        frame.setVisible(true);
        int base=0;
        while(true){
            for(int i=0; i< h; i++)
                for(int j=0; j< w; j++){
                    if(i>base && i<base+10){
                        img[i*w+j]=(255);
                    }else if(i<=base && i>base-10)
                        img[i*w+j]=(int)(255f*(i-base+10)/10);
                    else if(i<base+20 && i>=base+10)
                        img[i*w+j]=(int)(255f*(base+20-i)/10);
                    else
                        img[i*w+j]=0;
                }
            base+=1;
            base%=h;
            mi.updateRawData(img);
        }
    }
}
