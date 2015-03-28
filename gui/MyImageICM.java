import java.awt.*;
import java.awt.image.*;
import java.awt.color.ColorSpace;
import javax.swing.*;
import java.util.*;

public class MyImage{
    double[] rawData;
    byte[] imgData;
    //MemoryImageSource mis;
    IndexColorModel cm;
    BufferedImage bi;
    Image image;
    ArrayList<JComponent> parent = new ArrayList<JComponent>();

    private int nbin = 512;
    private static final double EPS=0.0000001;
    private static final int colorBits = 8;
    private static final int cmapSize = (1 << colorBits);

    private byte[] rgb = new byte[cmapSize];
    private int width, height;
    private int[] histArray;
    private boolean localHist = true;
    private boolean autoRange = true;
    private boolean showHist = true;
    private double[] range = new double[2];
    private double[] ir = new double[2]; //image range
    private Rectangle source;

    public MyImage(double[] rawData, int w, int h, boolean lh){
        pseudoCons(rawData,w,h,lh);
    }
    public MyImage(double[] rawData, int w, int h){
        pseudoCons(rawData, w, h, false);
    }
    public MyImage(int[] rawData, int w, int h){
        double[] temp = new double[rawData.length];
        for(int i=0; i < rawData.length; i++)
            temp[i] = rawData[i];
        pseudoCons(temp, w, h, false);
    }

    private void pseudoCons(double[] rawData, int w, int h, boolean lh){
        width = w; height = h;
        localHist = lh;
        imgData = new byte[rawData.length];
        this.rawData = rawData;
        getImageRange(rawData); // range of raw data
        setColorMap(new double[]{0, 255});     // range how image is shown 
        cm = new IndexColorModel(colorBits,cmapSize, rgb, rgb, rgb);
        bi = new BufferedImage(width,height,
                BufferedImage.TYPE_BYTE_INDEXED, cm);
        image=bi;
        updateImgData(imgData);
        source = new Rectangle(0,0,width,height);
        histArray = new int[nbin];
        if(showHist) updateHist(rawData);
    }

    public void addParent(JComponent comp){ parent.add(comp); }
    public void rmParent(JComponent comp){ parent.remove(comp); }

    public void getImageRange(final double[] rawData){
        ir[0]=rawData[0]; ir[1]=rawData[0];
        for(int i=1; i< rawData.length; i++)
            if(ir[0]> rawData[i])
                ir[0]=rawData[i];
            else if(ir[1] < rawData[i])
                ir[1]=rawData[i];
    }

    public void setColorMap(double[] ir){
        range[0] = ir[0]; range[1] = ir[1];
        double r = range[1]-range[0];
        for(int i=0; i < rgb.length; i++){
            if( i > range[1])
                rgb[i]=rgb.length-1;
            else if( i<range[0] )
                rgb[i]=0;
            else
                rgb[i] = (byte)((i-range[0])*(rgb.length-1)/r);
                //(byte)((i-0)*r/rgb.length+range[0]);
        }
        //IndexColorModel cm = new IndexColorModel(colorBits,
        //        cmapSize, rgb, rgb, rgb);
    }

    public void updateRawData(double[] rawData){
        getImageRange(rawData);
        this.rawData = rawData;
        if(autoRange) setColorMap(ir);
        else updateImgData(rawData);
        if(showHist){
            updateHist(rawData);
        }
    }

    public void updateRawData(int[] rawData){
        if(this.rawData==null)
            this.rawData = new double[rawData.length];
        for(int i=0; i < rawData.length; i++)
            this.rawData[i] = rawData[i];
        updateRawData(this.rawData);
    }

    private void updateImgData(double[] rawData){
        for(int i=0; i < rawData.length; i++){
            if(rawData[i]>=range[1])
                imgData[i]=(byte) 255;
            else if(rawData[i] < range[0])
                imgData[i]=0;
            else
                imgData[i] = (byte)(255f*(rawData[i]-range[0])
                        /(range[1]-range[0]));
            if(bi!=null) bi.setRGB(i/width,i%width,imgData[i]);
        }
        //if(mis!=null) mis.newPixels();
    }

    private void updateHist(final double[] rawData){
        if(!localHist){
            double span = (ir[1]-ir[0])*nbin/(nbin-EPS);
            for(int i=0; i< rawData.length; i++)
                histArray[(int)((rawData[i]-ir[0])*nbin/span)]++;
        }else{
        }
    }

    private Image setupImage(byte[] imgData){

        //ColorModel cm =
        //    new ComponentColorModel(
        //            ColorSpace.getInstance(ColorSpace.CS_GRAY),
        //            (new int[]{8}), false, true,
        //            Transparency.OPAQUE,
        //            DataBuffer.TYPE_BYTE
        //            );
        //mis = new MemoryImageSource(width, height, cm, imgData, 0, width);
        //mis.setAnimated(true);
        //System.out.println(mis);
        //System.out.println("sizeof(int)=" +
        //        DataBuffer.getDataTypeSize(DataBuffer.TYPE_INT));
        //if(parent==null)
        //    System.out.println("null");

        //image = (new JPanel()).createImage(mis);
        System.out.println(image);
        return image;
    }

    public void setRect(double lb, double tb, double w, double h){
        source.setRect(lb,tb,w,h);
        System.out.println("from MyImage setRect: " + source);
    }
    public void mvRectTo(int lb, int tb){
        source.setLocation(lb,tb);
    }

    public int[] getHistArray(){ return histArray; }
    public int getWidth(){ return width; }
    public int getHeight(){ return height; }
    public double[] getRange() {return range; }
    public double[] getRawRange() {return Arrays.copyOf(ir,2); }
    public Rectangle getWindow(){ return source; }
    public Image getImage(){ return image; }
}
