/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx;

import java.io.File;
import java.net.URL;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.lang.builder.ToStringBuilder;
import org.apache.commons.lang.time.StopWatch;

import ffx.ui.LogHandler;
import ffx.ui.MainPanel;

/**
 * The HeadlessMain class is the entry point to the command line version of
 * Force Field X.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class HeadlessMain {

    private static final Logger logger = Logger.getLogger(HeadlessMain.class.getName());

    /**
     * Main does some window initializations.
     */
    public HeadlessMain(File commandLineFile, List<String> argList, LogHandler logHandler) {
        stopWatch.start();
        // Create the MainPanel and MainMenu, then add them to the JFrame
        mainPanel = new MainPanel();
        logHandler.setMainPanel(mainPanel);
        mainPanel.initialize();
        // Finally, open the supplied file if necessary.
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
                                embeddedScript.openStream(), ".ffx"));
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
        if (System.getProperty("ffx.timer") != null) {
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
            logger.info(sb.toString());
        }
    }

    /**
     * Commons.Lang Style toString.
     */
    @Override
    public String toString() {
        ToStringBuilder toStringBuilder = new ToStringBuilder(this).append(
                "Up Time: " + stopWatch).append("Logger: " + logger.getName());
        return toStringBuilder.toString();
    }
    /**
     * This is the main application wrapper.
     */
    public MainPanel mainPanel;
    public static StopWatch stopWatch = new StopWatch();
}
