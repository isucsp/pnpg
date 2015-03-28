import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.awt.*;
import java.io.*;
import java.awt.image.*;
import javax.imageio.*;
import com.feradz.*;
import slider.*;
//import sun.awt.image.*;

public class ImgDialog extends JDialog{
    private static final int nbin = 512;
    private static final Dimension DIALOG_SIZE = new Dimension(1024/2, 1024/2);
    static final String[] icons = {"res/zoom_in_16.png", "res/zoom_out_16.png",
        "res/cursor_hand.png", "res/cursor.png","res/cursor_drag_hand.png"};
    // to host the histogram of the image.
    private ImgPanel imgPanel;
    private JToggleButton[] button = new JToggleButton[4];
    public ImgDialog(JFrame parent, String title, MyImage img){
        super(parent, title);
        //getContentPane().setLayout(new BorderLayout());
        int width = img.getWidth(), height = img.getHeight();

        Border raisedetched = BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
        Border loweredetched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
        Border raisedbevel = BorderFactory.createRaisedBevelBorder();
        Border loweredbevel = BorderFactory.createLoweredBevelBorder();
        Border blackline = BorderFactory.createLineBorder(Color.black);

        JPanel toolSection = new JPanel(new FlowLayout(FlowLayout.LEFT));
        //toolSection.setBorder(raisedbevel);
        ButtonGroup bg = new ButtonGroup();
        int gap = 1;
        Insets insert = new Insets(gap,gap,gap,gap);
        ActionListener al = new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                for(int i=0; i< button.length; i++)
                    if(button[i].equals(e.getSource())){
                        imgPanel.setTool(i);
                        System.out.format("%d-th button pressed\n",i);
                        break;
                    }}};
        for(int i=0 ; i < button.length; i++){
            button[i]=new JToggleButton( new ImageIcon(icons[i]));
            System.out.println(button[i].getMargin());
            button[i].setMargin(insert);
            button[i].addActionListener(al);
            bg.add(button[i]);
            toolSection.add(button[i]);
        }
        add(toolSection, BorderLayout.PAGE_START);

        JPanel imgSection = new JPanel(new BorderLayout());
        imgSection.setBorder(BorderFactory.createTitledBorder("Image Reconstructed"));

        imgPanel = new ImgPanel(img);
        add(imgSection, BorderLayout.CENTER);
        imgSection.add(imgPanel, BorderLayout.CENTER);

        HistSection histSection = new HistSection(img);
        add(histSection, BorderLayout.PAGE_END);
        //histSection.add(new Button("test"),BorderLayout.PAGE_END);

        pack();
        button[3].doClick();
    }

    public void updateImg(){
        imgPanel.update();
    }

    public void updateImgR(){
        //imgPanel.setupRange();
        imgPanel.update();
    }
    public static void main(String[] str){
        int w=512, h=w;
        double[] img=new double[w*h];

        File file = new File("test.png");
        BufferedImage bi;
        try{
            bi = ImageIO.read(file);
        w=bi.getWidth();
        w=512;
        h=bi.getHeight();
        if(w<h) h=w; else w=h;
        img=new double[w*h];
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
                img[i*w+j]=(double) (bi.getRGB(j,i) & 0xff);
            }
        }

        img=new double[w*h];
        for(int i=0; i<h; i++){
            for(int j=0; j<w; j++){
                img[i*w+j]=(i+j)  % 712;
            }
        }
        MyImage mi = new MyImage(img,w,h);
        ImgDialog test = new ImgDialog(null, "test", mi);
        test.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        test.setVisible(true);
        }catch(IOException e){ e.printStackTrace();}
    }
}

