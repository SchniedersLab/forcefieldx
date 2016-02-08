/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.ui;

import java.io.File;
import java.util.logging.Logger;

import org.apache.commons.lang3.SystemUtils;
import org.apache.commons.lang3.builder.ToStringBuilder;

/**
 * FFXExec encapsulates a native replacement for the JDK System.exec() method.
 * TINKER programs are executed in their own thread through a call to a Native
 * method called "FFXExec" which in turn calls the function "system()". The
 * reason we are not using the System.exec() methods is that some TINKER
 * routines execute indefinitely. Users may want to exit Force Field X and shut
 * down the JVM after launching a dynamics run, for example. In this case the
 * thread should not be dependent on a JVM instance.
 *
 * @author Michael J. Schnieders
 *
 */
public class FFXExec implements Runnable {

    private static final Logger logger = Logger.getLogger(FFXExec.class.getName());
    private static String path;
    private static String ld_library_path;
    private static String classpath;

    // Set up the PATH (to TINKER/bin), CLASSPATH and LD_LIBRARY_PATH
    private void setEnv() {
        path = MainPanel.ffxDir.getAbsolutePath();
        classpath = MainPanel.classpath;
        // java.home should be the jre directory.
        ld_library_path = System.getProperty("java.home", ".");
        if (!SystemUtils.IS_OS_WINDOWS) {
            ld_library_path = ld_library_path + "/lib/i386/client:"
                    + ld_library_path + "/lib/i386:" + ld_library_path
                    + "/lib/i386/native_threads";
        } else {
            ld_library_path = ld_library_path + "\\bin\\client";
            path = path + File.pathSeparator + ld_library_path;
        }
    }
    private FFXSystem system;
    private String name;
    private String args;
    private String dir;
    private MainPanel mainPanel;
    private File newFile;
    private boolean alive = true;
    private boolean openOnto;
    private int returnValue = 0;

    /**
     * Constructor
     *
     * @param s FFXSystem the Native command will execute on
     * @param n Name of the log file
     * @param a Command to execute
     * @param d Directory to launch the command in
     * @param m MainPanel
     * @param file File to open
     * @param o Load the resulting version file onto the passed FFXSystem
     */
    public FFXExec(FFXSystem s, String n, String a, String d, MainPanel m,
            File file, boolean o) {
        system = s;
        name = n;
        args = a;
        dir = d;
        mainPanel = m;
        newFile = file;
        openOnto = o;
        logger.info(toString());
    }

    /**
     * <p>
     * Getter for the field <code>returnValue</code>.</p>
     *
     * @return a int.
     */
    public int getReturnValue() {
        return returnValue;
    }

    /**
     * <p>
     * isAlive</p>
     *
     * @return a boolean.
     */
    public boolean isAlive() {
        return alive;
    }

    /**
     * nativeExec method for launching native executables
     *
     * @param argv String
     * @param dir String
     * @param path String
     * @param classpath String
     * @param jre String
     * @return int
     */
    private native int nativeExec(String argv, String dir, String path,
            String classpath, String jre);

    /**
     * Executes the native call to "System()" and notifies the ResultPanel upon
     * completion. This should only be called indirectly by Thread.Start()
     */
    public void run() {
        setEnv();
        if (args == null || dir == null || path == null || classpath == null
                || ld_library_path == null) {
            Logger.getLogger("ffx").severe(
                    "Native executable could not be executed." + "\nCommand: "
                    + args + "\nDIR: " + dir + "\nPATH: " + path
                    + "\nCLASSPATH: " + classpath
                    + "\nLD_LIBRARY_PATH: " + ld_library_path);
            return;
        }
        Logger.getLogger("ffx").info(
                "Native command invoked." + "\nCommand: " + args + "\nDIR: "
                + dir + "\nFFXExec - PATH: " + path + "\nCLASSPATH: "
                + classpath + "\nLD_LIBRARY_PATH: " + ld_library_path);
        returnValue = nativeExec(args, dir, path, classpath, ld_library_path);
        // Check for a bad return value
        if (returnValue < 0) {
            Logger.getLogger("ffx").warning(
                    "The following job exited with a failure status: "
                    + returnValue + "\n" + args);
        }
        // Open any created file and display the log.
        if (mainPanel != null) {
            if (newFile != null) {
                String[] labels = args.split(" +");
                String command = labels[0].toUpperCase() + " on "
                        + system.getFile().getName();
                if (openOnto) {
                    mainPanel.openOn(newFile, system, command);
                } else {
                    mainPanel.open(newFile, command);
                }
            }
            // mainPanel.getLogPanel().setDone(name);
        }
        alive = false;
    }

    /**
     * {@inheritDoc}
     *
     * Commons.Lang Style toString.
     */
    @Override
    public String toString() {
        ToStringBuilder toStringBuilder = new ToStringBuilder(this).append(path).append(classpath).append(ld_library_path);
        return toStringBuilder.toString();
    }
}
