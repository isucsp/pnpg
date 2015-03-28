/**
 * File: HistogramWindow.java
 * Date: January 5, 2013
 * Author: Ferad Zyulkyarov (feradz@gmail.com)
 * 
 * Implements a GUI for the displaying histograms.
 */
package com.feradz;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.LayoutManager;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.File;
import java.awt.dnd.DropTarget;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

public class HistogramWindow extends JFrame {

	private static final long serialVersionUID = 1L;
	
	/**
	 * This is the main content pane. The main content pane is organized in two
	 * columns. The left displaying the image inside {@link #_imagePane} and the 
	 * right displaying the the histogram inside {@link #_histogramPane}.
	 */
	private JPanel _contentPane;
	
	/**
	 * This is the content pain where the image is displayed. 
	 */
	private JPanel _imagePane;
	
	/**
	 * The content pane where the histograms are plotted.
	 */
	private JPanel _histogramPane;
	
	/**
	 * We display the image inside this label by setting it as the label's image icon.
	 */
	private JLabel _imageContainer;
	
	/**
	 * A label which shows the name of the file which histogram is drawn.
	 */
	private JLabel _imageFileNameLabel;
	
	private JLabel _dragAndDropLabel;
	
	/*
	 * These label objects are used as a container for the histogram plots
	 * which are drawn as BufferedImage objects. 
	 */
	private JLabel _luminanceHistContainer;
	private JLabel _grayHist1Container;
	private JLabel _grayHist2Container;
	private JLabel _redHistContainer;
	private JLabel _greenHistContainer;
	private JLabel _blueHistContainer;
	
	/*
	 * These labels are attached to the plotted histograms.
	 */
	private JLabel _luminanceHistLabel;
	private JLabel _grayHist1Label;
	private JLabel _grayHist2Label;
	private JLabel _redHistLabel;
	private JLabel _greenHistLabel;
	private JLabel _blueHistLabel;
	
	
	private static final int WINDOW_WIDTH = 1024;
	private static final int WINDOW_HEIGHT = 800;
	private static final int SCALE_LONG_EDGE = 600;
	private static final int PLOT_WIDTH = 256;
	private static final int PLOT_HEIGHT = 100;

	/**
	 * Create the frame.
	 */
	public HistogramWindow() {
		setTitle("Image Histogram");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, WINDOW_WIDTH, WINDOW_HEIGHT);
		
		//
		// initialize contentPane
		//
		_contentPane = new JPanel();
		_contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(_contentPane);
		GridBagLayout contentPaneLayout = new GridBagLayout();
		_contentPane.setLayout(contentPaneLayout);
		_contentPane.setBackground(Color.WHITE);
		
		//
		// initialize imagePane
		//
		
		_imagePane = new JPanel();
		_imagePane.setBackground(Color.WHITE);
		_imagePane.setBorder(new EmptyBorder(10, 10, 10, 10));
		//LayoutManager imagePaneLayout = new BoxLayout(imagePane, BoxLayout.Y_AXIS);
		LayoutManager imagePaneLayout = new BorderLayout(10, 10);
		_imagePane.setLayout(imagePaneLayout);
		
		_dragAndDropLabel = new JLabel();
		_imagePane.add(_dragAndDropLabel, BorderLayout.PAGE_START);
		JButton loadImageButton = new JButton("Load Image");
		_imagePane.add(loadImageButton, BorderLayout.LINE_START);
		_imageFileNameLabel = new JLabel();
		_imagePane.add(_imageFileNameLabel, BorderLayout.CENTER);
		_imageContainer = new JLabel();
		_imageContainer.setText("Drag and drop an image file to render its histogram.");
		_imagePane.add(_imageContainer, BorderLayout.PAGE_END);
		
		//
		// initialize histogramPane
		//
		_histogramPane = new JPanel();
		_histogramPane.setBorder(new EmptyBorder(10, 10, 10, 10));
		BoxLayout hitoramPaneLayout = new BoxLayout(_histogramPane, BoxLayout.Y_AXIS); 
		_histogramPane.setLayout(hitoramPaneLayout);
		_histogramPane.setBackground(Color.WHITE);
		_luminanceHistContainer = new JLabel();
		_grayHist1Container = new JLabel();
		_grayHist2Container = new JLabel();
		_redHistContainer = new JLabel();
		_greenHistContainer = new JLabel();
		_blueHistContainer = new JLabel();
		
		_luminanceHistLabel = new JLabel();
		_grayHist1Label = new JLabel();
		_grayHist2Label = new JLabel();
		_redHistLabel = new JLabel();
		_greenHistLabel = new JLabel();
		_blueHistLabel = new JLabel();
		
		_histogramPane.add(_luminanceHistLabel);
		_histogramPane.add(_luminanceHistContainer);
		
		_histogramPane.add(_grayHist1Label);
		_histogramPane.add(_grayHist1Container);
		
		_histogramPane.add(_grayHist2Label);
		_histogramPane.add(_grayHist2Container);
		
		_histogramPane.add(_redHistLabel);
		_histogramPane.add(_redHistContainer);
		
		_histogramPane.add(_greenHistLabel);
		_histogramPane.add(_greenHistContainer);
		
		_histogramPane.add(_blueHistLabel);
		_histogramPane.add(_blueHistContainer);
		
		//
		// add the imagePane and histogram pane to the contentPane
		//
		_contentPane.add(_imagePane);
		_contentPane.add(_histogramPane);
		
		DropTarget dt = new ImageFileDropTarget(this);
		_contentPane.setDropTarget(dt);
		
		//
		// Set event listener for the load image button.
		//
		loadImageButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JButton btn = (JButton)e.getSource();
				final JFileChooser fc = new JFileChooser();
				fc.setFileFilter(new ImageFileFilter());
				int res = fc.showOpenDialog(btn);
				if (res == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					try {
						loadImage(selectedFile.getPath());
						plotHistograms(selectedFile.getPath());						
					} catch (IOException ex) {
						System.err.println(ex.getMessage());
						ex.printStackTrace();
						JOptionPane.showMessageDialog(btn, 
								"There was a problem trying to open the image file " + selectedFile.getName() + ".\nProbably the file is corrupt or you do not have permission to read it.", 
								"Error", 
								JOptionPane.ERROR_MESSAGE);
					}
					
				}
			}
		});
	}
	
	/**
	 * Load an image.
	 * @param fileName the file name of the image to be displayed.
	 * @throws IOException an exception which may be thrown if the 
	 *    file does not exists or if there is a problem reading it.
	 */
	public void loadImage(String fileName) throws IOException {	
		File f = new File(fileName);
		BufferedImage img = ImageIO.read(f);
		int imgWidth = img.getWidth();
		int imgHeight = img.getHeight();
		int longEdge = Math.max(imgWidth, imgHeight);
		int scaleRatio = longEdge / SCALE_LONG_EDGE;
		int imgNewWidth = imgWidth / scaleRatio;
		int imgNewHeight = imgHeight / scaleRatio;
		ImageIcon imageIcon = new ImageIcon(img.getScaledInstance(imgNewWidth, imgNewHeight, Image.SCALE_SMOOTH));
		_imageContainer.setIcon(imageIcon);
		_imageFileNameLabel.setText(fileName);
	}
	
	/**
	 * Plots the histograms for the specified image file.
	 * @param fileName the file to plot the histograms.
	 * @throws IOException an exception which may be thrown if the 
	 *    file does not exists or if there is a problem reading it.
	 */
	public void plotHistograms(String fileName) throws IOException {
		Histogram h = Histogram.getHisrogram(fileName);
		
		int[] luminanceH = h.getLuminanceHistogram();
		int[] computedGrayH1 = h.getComputedGrayHistogram1();
		int[] computedGrayH2 = h.getComputedGrayHistogram2();
		int[] redH = h.getRedHistogram();
		int[] greenH = h.getGreenHistogram();
		int[] blueH = h.getBlueHistogram();
		
		if (_luminanceHistLabel.getText().equals(null) || _luminanceHistLabel.getText().length() == 0) {
			setHistogramLabels();
			_imageContainer.setText("");
		}
		
		BufferedImage plotLuminanceH = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.GRAY, luminanceH);
		BufferedImage plotGrayH1 = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.GRAY, computedGrayH1);
		BufferedImage plotGrayH2 = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.GRAY, computedGrayH2);
		BufferedImage plotRedH = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.RED, redH);
		BufferedImage plotGreenH = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.GREEN, greenH);
		BufferedImage plotBlueH = plotHistogram(PLOT_WIDTH, PLOT_HEIGHT, Color.BLUE, blueH);
		
		ImageIcon iconLuminanceH = new ImageIcon(plotLuminanceH);
		ImageIcon iconGrayH1 = new ImageIcon(plotGrayH1);
		ImageIcon iconGrayH2 = new ImageIcon(plotGrayH2);
		ImageIcon iconRedH = new ImageIcon(plotRedH);
		ImageIcon iconGreenH = new ImageIcon(plotGreenH);
		ImageIcon iconBlueH = new ImageIcon(plotBlueH);
		
		_luminanceHistContainer.setIcon(iconLuminanceH);
		_grayHist1Container.setIcon(iconGrayH1);
		_grayHist2Container.setIcon(iconGrayH2);
		_redHistContainer.setIcon(iconRedH);
		_greenHistContainer.setIcon(iconGreenH);
		_blueHistContainer.setIcon(iconBlueH);
	}
	
	/**
	 * Smoothens and normalizes a histogram and then plots the histogram
	 * into a {@link BufferedImage} object.
	 * @param width the width of the plot.
	 * @param height the height of the plot.
	 * @param c the color.
	 * @param histogram the histogram to be plotted.
	 * @return a {@link BufferedImage} object which contains the plot of the
	 * histogram.
	 */
	private BufferedImage plotHistogram(int width, int height, Color c, int[] histogram) {
		BufferedImage plot = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = plot.createGraphics();
		g.setColor(c);
		int barWidth = width/histogram.length;
		
		histogram = Histogram.smooth(histogram);
		histogram = Histogram.normalize(height, histogram);
		
		for (int i = 0; i < histogram.length; i++) {
			g.fillRect(i*barWidth, height - histogram[i], barWidth, histogram[i]);
		}
		
		return plot;
	}
	
	private void setHistogramLabels() {
		_luminanceHistLabel.setText("Luminance");
		_grayHist1Label.setText("Grayscale 1");
		_grayHist2Label.setText("Grayscale 2");
		_redHistLabel.setText("Red");
		_greenHistLabel.setText("Green");
		_blueHistLabel.setText("Blue");
		_dragAndDropLabel.setText("Drag and drop an image file to render its histogram.");
	}

}
