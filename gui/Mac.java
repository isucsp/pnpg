
import javax.swing.*;
import java.awt.event.*;

public class Mac extends JPanel{
    public Mac(){
        JButton button;
       
        // Start button implementation
        button=new JButton("Start");
        add(button);
        button.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                String mesg="call matlab function";
                System.out.println(mesg);
            }
        });

        // Stop Button implementation
        button=new JButton("Stop");
        add(button);
        button.addActionListener(new ActionListener(){
            public void actionPerformed(ActionEvent e) {
                String mesg="stop matlab function";
                System.out.println(mesg);
            }
        });


    }
    public static void main(String[] args){

    }
}
