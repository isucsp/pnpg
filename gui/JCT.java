
import javax.swing.*;
import java.awt.event.*;

public class JCT extends JFrame implements ActionListener{
    private JMenu[] menus = {
        new JMenu("New Session"),
    };
    private JMenuItem[] algs;
    private JTabbedPane tabbedPane = new JTabbedPane();
    public JCT(String title){
        super(title);
        int i=0;
        algs=new JMenuItem[MenuDescriptor.size];
        for(MenuDescriptor menuItem: MenuDescriptor.values()){
            algs[i]=new JMenuItem(menuItem.menu);
            //algs[i].setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, InputEvent.ALT_MASK));
            menus[0].add(algs[i]);
            algs[i].addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    JMenuItem m=(JMenuItem)e.getSource();
                    //System.out.println(e.getSource());
                    System.out.println("menu \"" + m.getText() + "\" was pressed.");
                    for(MenuDescriptor menuItem: MenuDescriptor.values())
                        if(menuItem.menu.equals(m.getText())){
                            tabbedPane.addTab(menuItem.abr, new MACPanel(JCT.this));
                            tabbedPane.setTabComponentAt(tabbedPane.getTabCount()-1,
                                new ButtonTabComponent(tabbedPane));
                            break;
                        }
                }
            });
            i++;
        }
        JMenuBar mb=new JMenuBar();
        for(i=0; i< menus.length; i++){
            mb.add(menus[i]);
        }
        setJMenuBar(mb);
        add(tabbedPane);
        algs[0].doClick();
    }
    public JCT(){
        this("CT Image Reconstruction Tool");
    }
    public void actionPerformed(ActionEvent e){

    }
    public static void main(String[] args) {
        JCT jctWin = new JCT("CT Image Reconstruction Tool");
        int windowWidth = 800;
        // Window width in pixels
        int windowHeight = 600;
        // Window height in pixels
        jctWin.setBounds(50, 100,
                // Set position
                windowWidth, windowHeight);
        // and size
        jctWin.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        jctWin.setLocationRelativeTo(null);
        jctWin.setVisible(true);
        System.out.println("main ends");
    }
}

