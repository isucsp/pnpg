import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import java.awt.color.*;
import javax.imageio.*;
import java.awt.*;
import java.io.*;
import java.awt.image.*;
import slider.*;

public class HistSection extends JPanel{
    private ImgPanel imgPanel;
    private JLabel histLabel = new JLabel();
    private double[] bound;
    private RangeSlider rangeSlider = new RangeSlider();
    private static final Dimension HIST_MIN_SIZE = new Dimension(512, 100);
    private static final Dimension HIST_PRE_SIZE = new Dimension(512, 100);
    private MyImage mi;

    public HistSection(MyImage img){
        super(new BorderLayout());
        //histSection.setBorder(BorderFactory.createLineBorder(Color.black));
        mi = img;

        //rangeSlider.setPreferredSize(
        //        new Dimension(getWidth(),
        //            rangeSlider.getPreferredSize().height)
        //        );
        //double[] range = img.getRange();
        //rangeSlider.setMinimum(0);
        //rangeSlider.setMaximum(255);
        //rangeSlider.setValue(0);
        //rangeSlider.setUpperValue(255);
        //
        //// Add listener to update display.
        //rangeSlider.addChangeListener(new ChangeListener() {
        //    public void stateChanged(ChangeEvent e) {
        //        RangeSlider slider = (RangeSlider) e.getSource();
        //        double[] ir = mi.getRawRange();
        //        System.out.format("get from slider: (%d, %d)\n",
        //            slider.getValue(), slider.getUpperValue());
        //        double range = slider.getMaximum()-slider.getMinimum();
        //        double miRange = ir[1]-ir[0];

        //        double[] newIr = new double[2];
        //        newIr[0] = (slider.getValue()-slider.getMinimum())/range*miRange+ir[0];
        //        newIr[1] = (slider.getUpperValue()-slider.getMinimum())/range*miRange+ir[0];
        //        System.out.format("set range from (%f, %f) to (%f, %f)\n",
        //            ir[0], ir[1], newIr[0], newIr[1]);
        //        mi.setRange(newIr);
        //    }
        //});

        HistPanel histPanel = new HistPanel(img, false);
        setBorder(BorderFactory.createTitledBorder("Histogram"));
        //setMinimumSize(HIST_MIN_SIZE);
        add(histPanel,BorderLayout.CENTER);
        //add(rangeSlider,BorderLayout.PAGE_END);
    }
}
