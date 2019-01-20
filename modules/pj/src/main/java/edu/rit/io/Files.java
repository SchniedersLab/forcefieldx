//******************************************************************************
//
// File:    Files.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.Files
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.io;

import java.io.File;

/**
 * Class Files provides static methods for various file related operations.
 *
 * @author Alan Kaminsky
 * @version 26-Nov-2007
 */
public class Files {

// Prevent construction.
    private Files() {
    }

// Exported operations.
    /**
     * Append the given rank to the given file. The rank goes before the file
     * extension if any. For example,
     * <code>Files.fileForRank(new&nbsp;File("image.pjg"),2)</code> returns a File
     * whose name is <code>"image_2.pjg"</code>;
     * <code>Files.fileForRank(new&nbsp;File("image"),2)</code> returns a File whose
     * name is <code>"image_2"</code>.
     *
     * @param file File.
     * @param rank Rank.
     * @return File with rank appended.
     */
    public static File fileForRank(File file,
            int rank) {
        return fileAppend(file, "_" + rank);
    }

    /**
     * Append the given rank to the given file name. The rank goes before the
     * file extension if any. For example,
     * <code>Files.fileNameForRank("image.pjg",2)</code> returns
     * <code>"image_2.pjg"</code>; <code>Files.fileNameForRank("image",2)</code> returns
     * <code>"image_2"</code>.
     *
     * @param filename File name.
     * @param rank Rank.
     * @return File name with rank appended.
     */
    public static String fileNameForRank(String filename,
            int rank) {
        return fileNameAppend(filename, "_" + rank);
    }

    /**
     * Append the given suffix to the given file. The suffix goes before the
     * file extension if any. For example,
     * <code>Files.fileAppend(new&nbsp;File("image.pjg"),"_new")</code> returns a
     * File whose name is <code>"image_new.pjg"</code>;
     * <code>Files.fileAppend(new&nbsp;File("image"),"_new")</code> returns a File
     * whose name is <code>"image_new"</code>.
     *
     * @param file File.
     * @param suffix Suffix.
     * @return File with suffix appended.
     */
    public static File fileAppend(File file,
            String suffix) {
        return new File(fileNameAppend(file.getPath(), suffix));
    }

    /**
     * Append the given suffix to the given file name. The suffix goes before
     * the file extension if any. For example,
     * <code>Files.fileNameAppend("image.pjg","_new")</code> returns
     * <code>"image_new.pjg"</code>; <code>Files.fileNameAppend("image","_new")</code>
     * returns <code>"image_new"</code>.
     *
     * @param filename File name.
     * @param suffix Suffix.
     * @return File name with suffix appended.
     */
    public static String fileNameAppend(String filename,
            String suffix) {
        int i = filename.lastIndexOf('.');
        return i == -1
                ? filename + suffix
                : filename.substring(0, i) + suffix + filename.substring(i);
    }

    /**
     * Prepend the given prefix to the given file. The prefix goes after the
     * directory if any. For example,
     * <code>Files.filePrepend(new&nbsp;File("/home/ark/image.pjg"),"new_")</code>
     * returns a File whose name is <code>"/home/ark/new_image.pjg"</code>;
     * <code>Files.filePrepend(new&nbsp;File("image.pjg"),"new_")</code> returns a
     * File whose name is <code>"new_image.pjg"</code>. The system-dependent file
     * name separator character is used to detect the end of the directory if
     * any.
     *
     * @param file File.
     * @param prefix Prefix.
     * @return File with prefix prepended.
     */
    public static File filePrepend(File file,
            String prefix) {
        return new File(fileNamePrepend(file.getPath(), prefix));
    }

    /**
     * Prepend the given prefix to the given file name. The prefix goes after
     * the directory if any. For example,
     * <code>Files.fileNamePrepend("/home/ark/image.pjg","new_")</code> returns
     * <code>"/home/ark/new_image.pjg"</code>;
     * <code>Files.fileNamePrepend("image.pjg","new_")</code> returns
     * <code>"new_image.pjg"</code>. The system-dependent file name separator
     * character is used to detect the end of the directory if any.
     *
     * @param filename File name.
     * @param prefix Prefix.
     * @return File name with prefix prepended.
     */
    public static String fileNamePrepend(String filename,
            String prefix) {
        int i = filename.lastIndexOf(File.separatorChar);
        return i == -1
                ? prefix + filename
                : filename.substring(0, i + 1) + prefix + filename.substring(i + 1);
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.io.Files <I>filename</I> <I>rank</I>
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		String filename = args[0];
//		int rank = Integer.parseInt (args[1]);
//		File file = new File (filename);
//		System.out.println
//			("Files.fileForRank (file, rank) = " +
//			  Files.fileForRank (file, rank));
//		System.out.println
//			("Files.fileNameForRank (filename, rank) = " +
//			  Files.fileNameForRank (filename, rank));
//		}
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.io.Files <I>filename</I> <I>str</I>
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		String filename = args[0];
//		String str = args[1];
//		File file = new File (filename);
//		System.out.println
//			("Files.fileAppend (file, str) = " +
//			  Files.fileAppend (file, str));
//		System.out.println
//			("Files.fileNameAppend (filename, str) = " +
//			  Files.fileNameAppend (filename, str));
//		System.out.println
//			("Files.filePrepend (file, str) = " +
//			  Files.filePrepend (file, str));
//		System.out.println
//			("Files.fileNamePrepend (filename, str) = " +
//			  Files.fileNamePrepend (filename, str));
//		}
}
