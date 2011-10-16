//******************************************************************************
//
// File:    Job.java
// Package: edu.rit.pj.job
// Unit:    Class edu.rit.pj.job.Job
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

package edu.rit.pj.job;

import edu.rit.io.Stdio;

import java.io.Externalizable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInput;
import java.io.ObjectOutput;
import java.io.PrintStream;

import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

import java.util.ArrayList;

/**
 * Class Job encapsulates a job and its attributes. A job is typically created
 * by a {@linkplain JobGenerator}. The job's attributes are:
 * <UL>
 * <P><LI>
 * Job number. A JobGenerator assigns a unique job number to each job it
 * creates.
 * <P><LI>
 * Description. An arbitrary string. The {@linkplain Runner} program prints the
 * job number and description whenever a job is started.
 * <P><LI>
 * Fully-qualified name of the main program class. To run a job, the main
 * program class's <TT>public static void main(String[])</TT> method is called.
 * <P><LI>
 * List of argument strings for the <TT>main()</TT> method.
 * <P><LI>
 * Standard input redirection. The standard input may be redirected from a file.
 * The job redirects the per-thread standard input stream in class {@linkplain
 * edu.rit.io.Stdio Stdio}. For redirection to work, the <TT>main()</TT> method
 * must use <TT>Stdio.in()</TT> for the standard input.
 * <P><LI>
 * Standard output redirection. The standard output may be redirected to a file,
 * either overwriting or appending to the file. The job redirects the per-thread
 * standard output stream in class {@linkplain edu.rit.io.Stdio Stdio}. For
 * redirection to work, the <TT>main()</TT> method must use <TT>Stdio.out()</TT>
 * for the standard output.
 * <P><LI>
 * Standard error redirection. The standard error may be redirected to a file,
 * either overwriting or appending to the file. The standard error may also be
 * redirected to the same place as the standard output. The job redirects the
 * per-thread standard error stream in class {@linkplain edu.rit.io.Stdio
 * Stdio}. For redirection to work, the <TT>main()</TT> method must use
 * <TT>Stdio.err()</TT> for the standard error.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 09-Oct-2010
 */
public class Job
	implements Runnable, Externalizable
	{

// Hidden data members.

	//private static final long serialVersionUID = 717404081766242374L;

	private int myJobNumber;
	private String myDescription;
	private String myMainClassName;
	private ArrayList<String> myArguments = new ArrayList<String>();
	private Stdin myStdinRedirection = Stdin.NONE;
	private Stdout myStdoutRedirection = Stdout.NONE;
	private Stderr myStderrRedirection = Stderr.NONE;
	private File myStdinFile = null;
	private File myStdoutFile = null;
	private File myStderrFile = null;

// Exported constructors.

	/**
	 * Construct a new uninitialized job. This constructor is for use only by
	 * object deserialization.
	 */
	public Job()
		{
		}

	/**
	 * Construct a new job. The job has the given job number and description.
	 * The job will run the <TT>public static void main(String[])</TT> method in
	 * the given class.
	 *
	 * @param  theJobNumber      Job number.
	 * @param  theDescription    Job description.
	 * @param  theMainClassName  Fully qualified name of main class.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMainClassName</TT> is null.
	 */
	public Job
		(int theJobNumber,
		 String theDescription,
		 String theMainClassName)
		{
		if (theMainClassName == null)
			{
			throw new NullPointerException
				("Job(): Main class name is null");
			}
		myJobNumber = theJobNumber;
		myDescription = theDescription;
		myMainClassName = theMainClassName;
		}

// Exported operations.

	/**
	 * Returns this job's number.
	 *
	 * @return  Job number.
	 */
	public int getJobNumber()
		{
		return myJobNumber;
		}

	/**
	 * Returns this job's description.
	 *
	 * @return  Job description.
	 */
	public String getDescription()
		{
		return myDescription;
		}

	/**
	 * Add the given argument string to this job.
	 *
	 * @param  arg  Argument string.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>arg</TT> is null.
	 */
	public void addArgument
		(String arg)
		{
		if (arg == null)
			{
			throw new NullPointerException
				("Job.addArgument(): Argument string is null");
			}
		myArguments.add (arg);
		}

	/**
	 * Read this job's standard input from the given file. It is an error if the
	 * file does not exist.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 */
	public void stdinFromFile
		(File file)
		{
		if (file == null)
			{
			throw new NullPointerException
				("Job.stdinFromFile(): File is null");
			}
		myStdinRedirection = Stdin.FILE;
		myStdinFile = file;
		}

	/**
	 * Store this job's standard output in the given file. The file is created
	 * if it does not exist. The file is overwritten if it does exist.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 */
	public void stdoutToFile
		(File file)
		{
		if (file == null)
			{
			throw new NullPointerException
				("Job.stdoutToFile(): File is null");
			}
		myStdoutRedirection = Stdout.FILE;
		myStdoutFile = file;
		}

	/**
	 * Append this job's standard output to the end of the given file. The file
	 * is created if it does not exist.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 */
	public void stdoutAppendToFile
		(File file)
		{
		if (file == null)
			{
			throw new NullPointerException
				("Job.stdoutAppendToFile(): File is null");
			}
		myStdoutRedirection = Stdout.FILE_APPEND;
		myStdoutFile = file;
		}

	/**
	 * Store this job's standard error in the given file. The file is created
	 * if it does not exist. The file is overwritten if it does exist.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 */
	public void stderrToFile
		(File file)
		{
		if (file == null)
			{
			throw new NullPointerException
				("Job.stderrToFile(): File is null");
			}
		myStderrRedirection = Stderr.FILE;
		myStderrFile = file;
		}

	/**
	 * Append this job's standard error to the end of the given file. The file
	 * is created if it does not exist.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 */
	public void stderrAppendToFile
		(File file)
		{
		if (file == null)
			{
			throw new NullPointerException
				("Job.stderrAppendToFile(): File is null");
			}
		myStderrRedirection = Stderr.FILE_APPEND;
		myStderrFile = file;
		}

	/**
	 * Redirect this job's standard error to the same place as this job's
	 * standard output.
	 */
	public void stderrToStdout()
		{
		myStderrRedirection = Stderr.STDOUT;
		myStderrFile = null;
		}

	/**
	 * Run this job.
	 * <OL TYPE=1>
	 * <P><LI>
	 * The per-thread standard input, standard output, and standard error in
	 * class {@linkplain edu.rit.io.Stdio Stdio} are redirected as necessary.
	 * <P><LI>
	 * This job's main class is found, using the calling thread's context class
	 * loader.
	 * <P><LI>
	 * The main class's <TT>public static void main(String[])</TT> method is
	 * found.
	 * <P><LI>
	 * The <TT>main()</TT> method is called, passing in this job's arguments.
	 * </OL>
	 * <P>
	 * If an I/O error occurs during Step 1, a RuntimeException wrapping an
	 * IOException is thrown. If any exception is thrown during Steps 2 through
	 * 4, an exception stack trace is printed on the per-thread standard error
	 * (which may have been redirected), but the exception is not propagated.
	 */
	public void run()
		{
		// Redirect standard input, standard output, and standard error.
		try
			{
			switch (myStdinRedirection)
				{
				case FILE:
					Stdio.in (new FileInputStream (myStdinFile));
					break;
				}
			switch (myStdoutRedirection)
				{
				case FILE:
					Stdio.out
						(new PrintStream
							(new FileOutputStream (myStdoutFile, false),
							 true));
					break;
				case FILE_APPEND:
					Stdio.out
						(new PrintStream
							(new FileOutputStream (myStdoutFile, true),
							 true));
					break;
				}
			switch (myStderrRedirection)
				{
				case FILE:
					Stdio.err
						(new PrintStream
							(new FileOutputStream (myStderrFile, false),
							 true));
					break;
				case FILE_APPEND:
					Stdio.err
						(new PrintStream
							(new FileOutputStream (myStderrFile, true),
							 true));
					break;
				case STDOUT:
					Stdio.err (Stdio.out());
					break;
				}
			}
		catch (IOException exc)
			{
			throw new RuntimeException (exc);
			}

		try
			{
			// Find main class.
			Class<?> mainClass =
				Class.forName
					(myMainClassName,
					 true,
					 Thread.currentThread().getContextClassLoader());

			// Find public static void main(String[]) method.
			Method mainMethod = mainClass.getMethod ("main", String[].class);
			int modifiers = mainMethod.getModifiers();
			if (! Modifier.isStatic (mainMethod.getModifiers()))
				{
				throw new NoSuchMethodException
					("Job.run(): Class "+mainClass.getName()+
					 " main() method is not static");
				}
			if (mainMethod.getReturnType() != Void.TYPE)
				{
				throw new NoSuchMethodException
					("Job.run(): Class "+mainClass.getName()+
					 " main() method does not return void");
				}

			// Call main() method with arguments.
			String[] args =
				myArguments.toArray (new String [myArguments.size()]);
			mainMethod.invoke (null, (Object) args);
			}
		catch (Throwable exc)
			{
			exc.printStackTrace (Stdio.err());
			}
		}

	/**
	 * Write this job to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeInt (myJobNumber);
		out.writeUTF (myDescription);
		out.writeUTF (myMainClassName);
		out.writeInt (myArguments.size());
		for (String arg : myArguments)
			{
			out.writeUTF (arg);
			}
		out.writeObject (myStdinRedirection);
		out.writeObject (myStdoutRedirection);
		out.writeObject (myStderrRedirection);
		out.writeObject (myStdinFile);
		out.writeObject (myStdoutFile);
		out.writeObject (myStderrFile);
		}

	/**
	 * Read this job from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if an object's class could not be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		int n;
		myJobNumber = in.readInt();
		myDescription = in.readUTF();
		myMainClassName = in.readUTF();
		myArguments.clear();
		n = in.readInt();
		for (int i = 0; i < n; ++ i)
			{
			myArguments.add (in.readUTF());
			}
		myStdinRedirection  = (Stdin)  in.readObject();
		myStdoutRedirection = (Stdout) in.readObject();
		myStderrRedirection = (Stderr) in.readObject();
		myStdinFile  = (File) in.readObject();
		myStdoutFile = (File) in.readObject();
		myStderrFile = (File) in.readObject();
		}

	}
