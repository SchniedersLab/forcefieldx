/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
import java.util.logging.Logger;

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
    public HeadlessMain(File commandLineFile, LogHandler logHandler) {
        stopWatch.start();
        // Create the MainPanel and MainMenu, then add them to the JFrame
        mainPanel = new MainPanel();
        logHandler.setMainPanel(mainPanel);
        mainPanel.initialize();
        // Finally, open the supplied file if necessary.
        if (commandLineFile != null) {
            if (commandLineFile.exists()) {
                mainPanel.open(commandLineFile, null);
            } else {
                logger.warning(commandLineFile.toString() + " was not found.");
            }
        }
        if (System.getProperty("ffx.timer") != null) {
            logger.info("\nStart-up Time (msec): " + stopWatch.getTime());
            Runtime runtime = Runtime.getRuntime();
            runtime.runFinalization();
            runtime.gc();
            long occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
            long KB = 1024;
            logger.info("\nIn-Use Memory   (Kb): " + occupiedMemory / KB + "\nFree Memory     (Kb): " + runtime.freeMemory() / KB + "\nTotal Memory    (Kb): " + runtime.totalMemory() / KB);
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
