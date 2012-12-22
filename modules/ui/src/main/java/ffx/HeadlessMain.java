/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */
package ffx;

import java.io.File;
import java.net.URL;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.lang3.builder.ToStringBuilder;
import org.apache.commons.lang3.time.StopWatch;

import ffx.ui.LogHandler;
import ffx.ui.MainPanel;

/**
 * The HeadlessMain class is the entry point to the command line version of
 * Force Field X.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class HeadlessMain {

    private static final Logger logger = Logger.getLogger(HeadlessMain.class.getName());

    /**
     * Complete initializations.
     *
     * @param commandLineFile a {@link java.io.File} object.
     * @param argList a {@link java.util.List} object.
     * @param logHandler a {@link ffx.ui.LogHandler} object.
     */
    public HeadlessMain(File commandLineFile, List<String> argList, LogHandler logHandler) {
        // Start a timer.
        stopWatch.start();

        // Construct the MainPanel, set it's LogHandler, and initialize then it.
        mainPanel = new MainPanel();
        logHandler.setMainPanel(mainPanel);
        mainPanel.initialize();

        // Open the supplied script file.
        if (commandLineFile != null) {
            if (!commandLineFile.exists()) {
                /**
                 * See if the commandLineFile is an embedded script.
                 */
                String name = commandLineFile.getName();
                name = name.replace('.', File.separatorChar);
                ClassLoader loader = getClass().getClassLoader();
                URL embeddedScript = loader.getResource("ffx/scripts/" + name + ".ffx");
                if (embeddedScript == null) {
                    embeddedScript = loader.getResource("ffx/scripts/" + name + ".groovy");
                }
                if (embeddedScript != null) {
                    try {
                        commandLineFile = new File(
                                FFXClassLoader.copyInputStreamToTmpFile(
                                embeddedScript.openStream(), commandLineFile.getName(), ".ffx"));
                    } catch (Exception e) {
                        logger.warning("Exception extracting embedded script "
                                + embeddedScript.toString() + "\n" + e.toString());
                    }
                }
            }
            if (commandLineFile.exists()) {
                mainPanel.getModelingShell().setArgList(argList);
                mainPanel.open(commandLineFile, null);
            } else {
                logger.warning(format("%s was not found.", commandLineFile.toString()));
            }
        }

        /**
         * Print start-up information.
         */
        if (logger.isLoggable(Level.FINE)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Start-up Time (msec): %s.", stopWatch.getTime()));
            Runtime runtime = Runtime.getRuntime();
            runtime.runFinalization();
            runtime.gc();
            long occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
            long KB = 1024;
            sb.append(format("\n In-Use Memory   (Kb): %d", occupiedMemory / KB));
            sb.append(format("\n Free Memory     (Kb): %d", runtime.freeMemory() / KB));
            sb.append(format("\n Total Memory    (Kb): %d", runtime.totalMemory() / KB));
            logger.fine(sb.toString());
        }
    }

    /**
     * {@inheritDoc}
     *
     * Commons.Lang Style toString.
     */
    @Override
    public String toString() {
        ToStringBuilder toStringBuilder = new ToStringBuilder(this).append(
                "Up Time: " + stopWatch).append("Logger: " + logger.getName());
        return toStringBuilder.toString();
    }
    /**
     * This is the main application container for both the GUI and CLI.
     */
    public MainPanel mainPanel;

    /**
     * Constant
     * <code>stopWatch</code>
     */
    public static StopWatch stopWatch = new StopWatch();
}
