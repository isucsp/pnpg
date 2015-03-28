
import javax.swing.*;
import javax.swing.border.*;
import java.awt.event.*;
import java.awt.*;
import java.util.*;

public class MACPanel extends JPanel{
    MatClient mClient;
    ImgDialog imgDialog=null;
    JLabel statusPanel;
    JFrame parent;
    private ImgPanel imgPanel;
    JTextField urlAddress, port;
    JPanel mainPanel;
    JPanel imgSection, rSection, bSection, brSection;
    double imgPercent=0.8;
    MyImage mi;
    public static final int IMG_ARRAY_SIZE=1024*1024;
    double[] img= new double[IMG_ARRAY_SIZE];
    public MACPanel(JFrame parent){
        this.parent=parent;
        setLayout(new BorderLayout());
        JButton button;

        JPanel panel = new JPanel();
        panel.setPreferredSize(new Dimension(250, parent.getHeight()-50));
        add(panel,BorderLayout.LINE_START);
        panel.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.gridwidth = GridBagConstraints.REMAINDER;
        c.weightx = 1;
        c.fill = GridBagConstraints.HORIZONTAL;

        JPanel paraPanel = new ParameterPanel();
        panel.add(paraPanel,c);

        JPanel viewPanel = new JPanel();
        viewPanel.setBorder(BorderFactory.createTitledBorder("View Options"));
        //viewPanel.setPreferredSize(CONNECTION_PANEL_SIZE);
        //viewPanel.setSize(CONNECTION_PANEL_SIZE);
        JCheckBox keepUpdating = new JCheckBox("Updating reconstruction");
        viewPanel.add(keepUpdating);
        //keepUpdating.addItemListener(new ItemListener(){
        //    public void itemStateChanged(ItemEvent e) {
        //        imgDialog.setVisible(e.getStateChange()==ItemEvent.SELECTED);
        //    }
        //});
        panel.add(viewPanel,c);

        c.weighty=1; c.anchor=GridBagConstraints.PAGE_START;
        JPanel ctrlPanel = new ControlPanel();
        panel.add(ctrlPanel,c);

        panel = new JPanel();
        panel.setBorder(BorderFactory.createEmptyBorder());
        //panel.setBorder(BorderFactory.createEtchedBorder());
        //Border raisedbevel = BorderFactory.createRaisedBevelBorder();
        //Border loweredbevel = BorderFactory.createLoweredBevelBorder();
        ((FlowLayout)panel.getLayout()).setVgap(0);
        ((FlowLayout)panel.getLayout()).setAlignment(FlowLayout.LEFT);
        //panel.setBackground(Color.cyan);
        statusPanel = new JLabel("Ready...");
        panel.add(statusPanel);
        //panel.setPreferredSize(new Dimension(200, 30));
        add(panel,BorderLayout.PAGE_END);

        panel = new JPanel();
        panel.setBorder(BorderFactory.createEmptyBorder());
        //((FlowLayout)panel.getLayout()).setVgap(0);
        //((FlowLayout)panel.getLayout()).setAlignment(FlowLayout.LEFT);

        JLabel label = new JLabel("Server:");
        urlAddress = new JTextField("localhost",20);
        panel.add(label); panel.add(urlAddress);
        label = new JLabel("Port:");
        port = new JTextField("30000",5);
        panel.add(label); panel.add(port);
        add(panel,BorderLayout.PAGE_START);

        int w = (int)Math.sqrt(IMG_ARRAY_SIZE), h;
        if(w*w==IMG_ARRAY_SIZE) h=w;
        else{
            h=IMG_ARRAY_SIZE/w;
            System.err.println("Img size incorrect");
            return;
        }
        mi = new MyImage(img,w,h);
        imgDialog=new ImgDialog(parent, "Reconstruected Image", mi);

        button=new JButton("Connect");
        panel.add(button);
        button.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                //new Thread(new Runnable() {
                //SwingUtilities.invokeLater(new Runnable() {
                //    public void run() {
                        String url = urlAddress.getText();
                        int port = Integer.parseInt(MACPanel.this.port.getText());
                        mClient=new MatClient(url,port,statusPanel);
                        statusPanel.setText("Connecting to MATLAB ...");
                        mClient.connect();
                        //try{ Thread.sleep(20000); }catch(Exception eee){}
                        statusPanel.setText("Connected to MATLAB.");
                //    }
                //});
                //}).start();
            }
        });
        button=new JButton("Disconnect");
        panel.add(button);
        button.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                statusPanel.setText("Disconnecting to MATLAB ...");
                mClient.disconnect();
                mClient = null;
                statusPanel.setText("Ready.");
            }
        });

        mainPanel = new JPanel();
        mainPanel.setBackground(Color.cyan);
        mainPanel.setLayout(null);
        add(mainPanel, BorderLayout.CENTER);

        imgSection = new JPanel(new GridBagLayout());
        imgSection.setBorder(BorderFactory.createTitledBorder("Image Reconstructed"));
        imgPanel = new ImgPanel(mi);
        c = new GridBagConstraints();
        c.weightx=1; c.weighty=1; c.fill=GridBagConstraints.BOTH;
        imgSection.add(imgPanel,c);
        mainPanel.add(imgSection);

        rSection = new JPanel();
        rSection.setBorder(BorderFactory.createTitledBorder("Vertical profile"));
        mainPanel.add(rSection);

        bSection = new JPanel();
        bSection.setBorder(BorderFactory.createTitledBorder("Horizontal profile"));
        mainPanel.add(bSection);

        brSection = new HistSection(mi);
        mainPanel.add(brSection);
    }
    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        int aa = (int)(imgPercent*Math.min(mainPanel.getWidth(),mainPanel.getHeight()));
        System.out.format("size of image: %d\n", aa);
        imgSection.setBounds(0,0,aa,aa);
        System.out.format("size of imgSection: (%d, %d)\n", imgSection.getWidth(), imgSection.getHeight());
        rSection.setBounds(aa,0,mainPanel.getWidth()-aa,aa);
        bSection.setBounds(0,aa,aa,mainPanel.getHeight()-aa);
        brSection.setBounds(aa,aa,mainPanel.getWidth()-aa,mainPanel.getHeight()-aa);
    }
    void start(){
    }
    public static void main(String[] args){

    }

    class ParameterPanel extends JPanel{
        final int gaps=3;
        ArrayList<JComponent> comp = new ArrayList<JComponent>();
        public ParameterPanel(){
            GridBagConstraints cl = new GridBagConstraints();
            GridBagConstraints cr = new GridBagConstraints();
            cl.weightx=0.5; cl.anchor=GridBagConstraints.LINE_END;
            cl.weighty=0.0; cl.fill=GridBagConstraints.NONE;
            cl.insets = new Insets(gaps,0,gaps,gaps);
            cr.weightx=0.5; cr.anchor=GridBagConstraints.LINE_START;
            cr.weighty=0.0; cr.fill=GridBagConstraints.HORIZONTAL;
            cr.gridwidth = GridBagConstraints.REMAINDER;
            cr.insets = new Insets(gaps,gaps,gaps,gaps);

            setLayout(new GridBagLayout());
            setBorder(BorderFactory.createTitledBorder("Parameters"));
            //System.out.println(parent.getSize());

            String[] maskTypeStr = { "Circular", "Convex Hull", "Tight"};
            JComboBox<String> comboBox = new JComboBox<String>(maskTypeStr);
            comboBox.setSelectedIndex(0);
            //petList.addActionListener(this);
            add(new JLabel("Mask Type:"),cl);
            add(comboBox,cr);
            comp.add(comboBox);

            String[] imageNameStr= {"twoMaterials", "realct", "pellet", "castSim", "phantom"};
            comboBox = new JComboBox<String>(imageNameStr);
            add(new JLabel("Image Name:"),cl);
            add(comboBox,cr);
            comp.add(comboBox);

            JPanel panel = new JPanel();
            //panel.setBorder(BorderFactory.createTitledBorder("Parameters"));
            //panel.setBackground(Color.cyan);
            ((FlowLayout)panel.getLayout()).setVgap(0);
            ButtonGroup bg = new ButtonGroup();
            JRadioButton spark = new JRadioButton("with",true);
            JRadioButton nonSpark = new JRadioButton("without",true);
            bg.add(spark); bg.add(nonSpark);
            panel.add(spark); panel.add(nonSpark);
            add(new JLabel("Characteristic Lines:"),cr);
            cl.gridwidth = GridBagConstraints.REMAINDER;
            add(panel,cl);
            comp.add(spark); comp.add(nonSpark);
            cl.gridwidth = 1;

            String[] initSignalStr={"FBP recon.", "zeros"};
            add(new JLabel("Initialization:"),cl);
            add(comboBox=new JComboBox<String>(initSignalStr),cr);

            String[] jtfLables = {"Max # of Itr",
                "<html><font face=serif><i>v</i></font></html>",
                "<html><font face=serif><i>u</i></font></html>",
                "<html><font face=serif><i>J</i></font></html>",
                "log span", "Refinement"};
            String[] jtfValue = { "2000", "1", "1e-4", "17", "1e3", "1"};
            JTextField textField;
            for(int i=0; i< jtfLables.length; i++){
                add(new JLabel(jtfLables[i]+":"),cl);
                textField = new JTextField(jtfValue[i],100);
                add(textField,cr);
                comp.add(textField);
            }
        }
    }
    class ControlPanel extends JPanel{
        final String[] icons = {"res/zoom_in_16.png", "res/zoom_out_16.png",
            "res/cursor_hand.png", "res/cursor.png","res/cursor_drag_hand.png"};
        // to host the histogram of the image.
        private JToggleButton[] button = new JToggleButton[4];
        public ControlPanel(){
            //Border raisedetched = BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
            //Border loweredetched = BorderFactory.createEtchedBorder(EtchedBorder.LOWERED);
            //Border raisedbevel = BorderFactory.createRaisedBevelBorder();
            //Border loweredbevel = BorderFactory.createLoweredBevelBorder();
            //Border blackline = BorderFactory.createLineBorder(Color.black);
            setBorder(BorderFactory.createTitledBorder("Control Panel"));
            setLayout(new GridBagLayout());
            GridBagConstraints c = new GridBagConstraints();
            c.gridwidth = GridBagConstraints.REMAINDER;

            JButton startButton=new JButton("Start");
            System.out.println("MACPanel current thread ID: " + Thread.currentThread().getId());
            startButton.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e) {
                    statusPanel.setText("Calling matlab function ...");
                    mClient.exec("setupPath");
                    mClient.exec("close all");
                    mClient.exec("runRealCT");
                    //EventQueue.invokeLater(
                    (new Thread(){
                        // start to update the progress information in dialog windows
                        public void run(){
                            double rse, cost;
                            int itr;
                            System.out.println("run: " + Thread.currentThread().getId());
                            while(true){
                                mClient.exec("progress","");
                                String ret="image"; //mClient.retString().trim().toLowerCase();
                                System.out.println("get: " + ret);
                                if(ret.equals("image")){
                                    //mClient.retArray(img,IMG_ARRAY_SIZE);
                                    System.out.println(mClient.retArray(img,IMG_ARRAY_SIZE) +
                                        "doubles: " + img[0]+", "+img[1] + " ...");
                                }else{
                                    if(ret.equals("idle")){
                                        break;
                                    }
                                }
                                itr=(int) mClient.readDouble();
                                cost = mClient.readDouble();
                                rse = mClient.readDouble();
                                mi.updateRawData();
                                imgPanel.addStatistics(itr, cost, rse);
                            }
                        }
                    }).start();
                    System.out.println("ActionListener: " + Thread.currentThread().getId());
                    System.out.println("System action performed");
                }
            });

            // Stop Button implementation
            JButton stopButton=new JButton("Stop");
            stopButton.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e) {
                    String mesg="stop matlab function";
                    System.out.println(mesg);
                }
            });

            JPanel toolSection = new JPanel();
            toolSection.add(startButton);
            toolSection.add(stopButton);

            add(toolSection,c);

            toolSection = new JPanel(new FlowLayout(FlowLayout.LEFT));
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
            add(toolSection, c);
            button[3].setSelected(true);

            /*
               JPanel imgSection = new JPanel(new BorderLayout());
               imgSection.setBorder(BorderFactory.createTitledBorder("Image Reconstructed"));

               imgPanel = new ImgPanel(img);
               add(imgSection, BorderLayout.CENTER);
               imgSection.add(imgPanel, BorderLayout.CENTER);

               HistSection histSection = new HistSection(img);
               add(histSection, BorderLayout.PAGE_END);
            //histSection.add(new Button("test"),BorderLayout.PAGE_END);
            */
        }
    }
}
