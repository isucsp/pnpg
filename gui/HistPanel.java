
import javax.swing.*;
import javax.swing.border.*;
import java.awt.color.*;
import javax.imageio.*;
import java.awt.*;
import java.io.*;
import java.awt.image.*;

public class HistPanel extends JPanel implements Runnable{
    private static final double EPS=0.0000001;
    private static final Dimension HIST_MIN_SIZE = new Dimension(512, 100);
    private static final Dimension HIST_PRE_SIZE = new Dimension(512, 100);
    private int MARGIN=6;

    MyImage mi;
    int[] histArray;
    int max;
    double barWidth;
    int height, width;
    double[] range= new double[2];
    public HistPanel(MyImage mi, boolean vertical){
        this.mi = mi;
        setMinimumSize(HIST_MIN_SIZE);
        setPreferredSize(HIST_PRE_SIZE);
        setBorder(BorderFactory.createLineBorder(Color.black));
        new Thread(this).start();
    }

    private int getMax(int[] h){
        int max=h[0];
        for (int i = 1; i < h.length; i++)
            if(h[i]>max) max=h[i];
        return max;
    }

    public void run(){
        while(true){
            if(mi.hasNewHist()){
                repaint();
            }else
                try{Thread.sleep(50);}catch(InterruptedException ee){
                    System.out.println(ee);
                }
        }
    }

    protected void paintComponent(Graphics g){
        super.paintComponent(g);

        int width=getWidth()-2*MARGIN, height=getHeight();
        int[] histArray = mi.getHistArray();
        barWidth = 0.8*width/histArray.length;
        double h;
        System.out.println("Histogram: ("+width+", "+height+") "
                + barWidth + ", "+histArray[200]);
        Graphics2D g2d = (Graphics2D)g; 
        g2d.setColor(Color.RED);
        //histogram = Histogram.smooth(histogram);

        int nbin = histArray.length;
        int max = getMax(histArray);
        for (int i = 0; i < nbin; i++) {
            h=((double)histArray[i])*height/max;
            g2d.fill(new Rectangle.Double(-barWidth/2+(i+0.5)*width/nbin+MARGIN, height - h,
                    barWidth, h));
        }
    }
}

