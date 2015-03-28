/**
 * File: Histogram.java
 * Date: December 30, 2012
 * Author: Ferad Zyulkyarov (feradz@gmail.com)
 * 
 * A program which plots a gray scale histogram and histograms
 * for the red, green and blue channels of an image.
 * 
 * It accepts one command line argument - an image file name.
 */
package com.feradz;

import java.awt.EventQueue;
import java.io.IOException;

public class HistogramMain {
	
	/**
	 * The main method (entry point) of the program.
	 * @param args command line arguments.
	 */
	public static void main(final String[] args) {		
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					HistogramWindow frame = new HistogramWindow();
					frame.setVisible(true);
					if (args.length == 1) {
						String fileName = args[0];
						frame.loadImage(fileName);
						frame.plotHistograms(fileName);
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		});
	}
}
