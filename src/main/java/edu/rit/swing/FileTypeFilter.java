//******************************************************************************
//
// File:    FileTypeFilter.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.FileTypeFilter
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.swing;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * Class FileTypeFilter implements a Swing file filter that accepts only files
 * whose file name matches a specified suffix.
 *
 * @author  Alan Kaminsky
 * @version 20-Jul-2006
 */
public class FileTypeFilter
	extends FileFilter
	{

// Hidden data members.

	private String myExtension;
	private int myExtensionLength;
	private String myDescription;

// Exported constructors.

	/**
	 * Construct a new file type filter. Initially, the file type filter has no
	 * file type set.
	 */
	public FileTypeFilter()
		{
		}

// Exported operations.

	/**
	 * Set the file extension that this file type filter will accept. Any file
	 * whose file name ends with <TT>theExtension</TT> (case insensitive) will
	 * be accepted.
	 *
	 * @param  theExtension    File extension, e.g. <TT>".png"</TT>.
	 * @param  theDescription  Description of files accepted, e.g. <TT>"PNG
	 *                         image files"</TT>.
	 */
	public void setType
		(String theExtension,
		 String theDescription)
		{
		myExtension = theExtension;
		myExtensionLength = theExtension.length();
		myDescription = theDescription;
		}

	/**
	 * Determine whether this file type filter accepts the given file.
	 *
	 * @param  f  File.
	 *
	 * @return  True if <TT>f</TT> is accepted, false otherwise.
	 */
	public boolean accept
		(File f)
		{
		String name = f.getName();
		int n = name.length();
		return
			f.isDirectory() ||
			(n >= myExtensionLength &&
				name.substring (n - myExtensionLength)
					.equalsIgnoreCase (myExtension));
		}

	/**
	 * Get a description of the files this file type filter accepts.
	 *
	 * @return  Description.
	 */
	public String getDescription()
		{
		return myDescription;
		}

	}
