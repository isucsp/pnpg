import javax.swing.*;
import java.awt.color.*;
import javax.imageio.*;
import java.awt.*;
import java.io.*;
import java.awt.image.*;

public class ImgPanel extends JPanel{
    double[] raw;
    int[] iraw;
    double pMin, pMax;
    int width, height;
    byte[] img;
    DataBuffer db;
    SampleModel sm;
    ColorModel cm;
    Image cImage;
    MemoryImageSource mis;
    WritableRaster wr;
    public ImgPanel(double[] raw, int width, int height){
        this.raw=raw;
        img= new byte[raw.length];
        pMin=raw[0]; pMax=raw[0];
        for(int i=1; i< raw.length; i++){
            if(raw[i] < pMin) pMin=raw[i];
            if(raw[i] > pMax) pMax=raw[i];
        }
        this.width=width; this.height=height;
    }

    public ImgPanel(int[] raw, int width, int height){
        iraw=raw;
        img= new byte[raw.length];
        int pMin=raw[0], pMax=raw[0];
        for(int i=1; i< raw.length; i++){
            if(raw[i] < pMin) pMin=raw[i];
            if(raw[i] > pMax) pMax=raw[i];
        }
        this.pMax = (double)pMax; this.pMin=(double)pMin;
        this.width=width; this.height=height;
    }

    public void setupRange(double min, double max){
        pMin=min; pMax=max;
    }
    public void setupRange(int min, int max){
        pMin=(double)min; pMax=(double)max;
    }

    private void mapColor(){
        if(raw!=null)
            for(int i=0; i< raw.length; i++){
                if(raw[i]>pMax)
                    img[i]=(byte) 255;
                else if(raw[i] < pMin)
                    img[i]=0;
                else
                    img[i] = (byte)(255f*(raw[i]-pMin)/(pMax-pMin));
            }
        else
            for(int i=0; i< iraw.length; i++){
                if(iraw[i]>pMax)
                    img[i]=(byte) 255;
                else if(iraw[i] < pMin)
                    img[i]=0;
                else
                    img[i] = (byte)(255f*(iraw[i]-pMin)/(pMax-pMin));
            }
    }

    public void setupImage(){
        mapColor();
        cm = new ComponentColorModel(ColorSpace.getInstance(ColorSpace.CS_GRAY),
                (new int[]{8}), false, true,
                Transparency.OPAQUE,
                DataBuffer.TYPE_BYTE);
        mis = new MemoryImageSource(width, height, cm, img, 0, width);
        mis.setAnimated(true);
        System.out.println(mis);
        System.out.println("sizeof(int)=" + DataBuffer.getDataTypeSize(DataBuffer.TYPE_INT));
        setPreferredSize(new Dimension(width,height));
        cImage=createImage(mis);
        System.out.println(cImage);
        //getGraphics().drawImage(cImage,0,0,this);
        //for(int i=0; i< img2.length; i++) img2[i]=(0xff<<24 | 127);
        //add(new JLabel(new ImageIcon(cImage)));
    }

    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        g.drawImage(cImage,0,0,this);
    }

    void update(int[] raw){
        iraw=raw;
        mapColor();
        mis.newPixels(0,0,width,height);
    }

    public static void main(String[] str){
        int w=500, h=500;
        int[] img = new int[w*h];
        JFrame frame = new JFrame("test ImgPanel");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        ImgPanel test = new ImgPanel(img, w, h);
        test.setupImage();
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
            test.update(img);
        }
    }
}
