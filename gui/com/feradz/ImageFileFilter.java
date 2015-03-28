/**
 * File: ImageFileFilter.java
 * Date: January 5, 2013
 * Author: Ferad Zyulkyarov (feradz@gmail.com)
 * 
 * Implements a file filter for image files.
 */
package com.feradz;

import java.io.File;

import javax.swing.filechooser.FileFilter;

public class ImageFileFilter extends FileFilter {

	public static final String[] SUPPORTED_EXTENSIONS = {"jpg"};
	public static final String SUPPORTED_EXTENSIONS_DESCR = "jpg";
	
	@Override
	public boolean accept(File f) {
		boolean isSupportedFileType = false;
		
		if (f.isDirectory()) {
			return true;
		}
		
		String fileExt = getFileExtensions(f);
		for (String s: SUPPORTED_EXTENSIONS) {
			if (fileExt.equals(s)) {
				isSupportedFileType = true;
				break;
			}
		}
		
		return isSupportedFileType;
	}

	@Override
	public String getDescription() {
		return SUPPORTED_EXTENSIONS_DESCR;
	}
	
	private static String getFileExtensions(File f) {
		String ext = "";
		String fileName = f.getName();
		int extDotIndex = fileName.lastIndexOf('.');
		if (extDotIndex > 0 && extDotIndex < fileName.length() - 1) {
			ext = fileName.substring(extDotIndex + 1);
			ext = ext.toLowerCase();
		}
		
		return ext;
	}

}
