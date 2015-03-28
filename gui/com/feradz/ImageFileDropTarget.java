/**
 * File: ImageFileDroptTarget.java
 * Date: January 6, 2013
 * Author: Ferad Zyulkyarov (feradz@gmail.com)
 * 
 * Implements a drop target in order to open files through drag and drop.
 */

package com.feradz;

import java.awt.datatransfer.DataFlavor;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDropEvent;
import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileFilter;

public class ImageFileDropTarget extends DropTarget{
	private static final long serialVersionUID = 1L;
	private HistogramWindow _frame;
	
	public ImageFileDropTarget(JFrame frame) {
		super();
		_frame = (HistogramWindow)frame;
	}
	
	public synchronized void drop(DropTargetDropEvent evt) {
		try {
			evt.acceptDrop(DnDConstants.ACTION_COPY);
			
			@SuppressWarnings("unchecked")
			List<File> droppedFiles = (List<File>)evt.getTransferable().getTransferData(DataFlavor.javaFileListFlavor);
			
			if (droppedFiles.size() > 1) {
				JOptionPane.showMessageDialog(
						_frame, 
						"You must select only one file.", 
						"Error", 
						JOptionPane.ERROR_MESSAGE);
			}
			
			File droppedFile = droppedFiles.get(0);
			FileFilter fileFilter = new ImageFileFilter();
			if (fileFilter.accept(droppedFile)) {
				String droppedFileName = droppedFile.getPath();
				try {
					_frame.loadImage(droppedFileName);
					_frame.plotHistograms(droppedFileName);
				} catch (IOException ex) {
					JOptionPane.showMessageDialog(_frame, 
							"There was a problem trying to open the image file " + droppedFileName + ".\nProbably the file is corrupt or you do not have permission to read it.", 
							"Error", 
							JOptionPane.ERROR_MESSAGE);
				}
			}
			else {
				JOptionPane.showMessageDialog(
						_frame, 
						"This file type is not supported. Supported file types are: " + fileFilter.getDescription(), 
						"Error", 
						JOptionPane.ERROR_MESSAGE);
			}
		} catch (Exception e) {
			JOptionPane.showMessageDialog(_frame, 
					"Bad nasty error happened. " + e.getMessage(), 
					"Error", 
					JOptionPane.ERROR_MESSAGE);
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}
}
