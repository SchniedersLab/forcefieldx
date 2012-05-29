/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
package ffx.ui;

import java.awt.GraphicsEnvironment;
import java.util.logging.ErrorManager;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;

/**
 * The default ConsoleHanlder publishes logging to System.err. This class
 * publishes to System.out, which is normally intercepted by the Force Field X
 * Shell. The formatter used reduces verbosity relative to the default
 * SimpleFormatter.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class LogHandler extends Handler {

    private static final boolean headless = GraphicsEnvironment.isHeadless();

    private MainPanel mainPanel = null;
    private boolean fatal = false;

    /**
     * A reference to the Force Field X MainPanel container to shut down
     * if we encounter a fatal (SEVERE) exception. If we are not in
     * Headless mode, then LogRecords can be published to the ModelingShell.
     *
     * @param mainPanel the Force Field X MainPanel.
     * @since 1.0
     */
    public void setMainPanel(MainPanel mainPanel) {
        this.mainPanel = mainPanel;
    }

    /**
     * Construct the Force Field X logging handler.
     *
     * @since 1.0
     */
    public LogHandler() {
        setFormatter(new LogFormatter(false));
        setLevel(Level.ALL);
    }

    /**
     * {@inheritDoc}
     *
     * Publish a LogRecord.
     * @since 1.0.
     */
    @Override
    public synchronized void publish(LogRecord record) {
        /**
         * Check if the record is loggable and that we have not already
         * encountered a fatal error.
         */
        if (!isLoggable(record) || fatal) {
            return;
        }
        String msg;
        try {
            msg = getFormatter().format(record);
        } catch (Exception e) {
            /**
             * We don't want to throw an exception here, but we
             * report the exception to any registered ErrorManager.
             */
            reportError(null, e, ErrorManager.FORMAT_FAILURE);
            return;
        }
        try {
            if (record.getLevel() == Level.SEVERE) {
                fatal = true;
                System.err.println(msg);
                System.err.println(" Force Field X will not continue.");
                System.err.println(" Shutting down...");
                mainPanel.exit();
            }
            ModelingShell shell = null;
            Thread thread = null;
            if (mainPanel != null) {
                shell = mainPanel.getModelingShell();
                if (shell != null) {
                    thread = shell.getRunThread();
                }
            }
            if (!headless && thread != null && thread.isAlive() && !thread.isInterrupted()) {
                shell.appendOutputNl(msg, shell.getResultStyle());
            } else {
                System.out.println(msg);
            }
        } catch (Exception e) {
            /**
             * We don't want to throw an exception here, but we
             * report the exception to any registered ErrorManager.
             */
            reportError(null, e, ErrorManager.WRITE_FAILURE);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void flush() {

        if (mainPanel.getModelingShell() == null) {
            System.out.flush();
        } else {
            // Scroll to visible!
        }
    }

    /**
     * {@inheritDoc}
     *
     * Flush, but do not close System.out or the Shell.
     * @since 1.0
     */
    @Override
    public void close() {
        flush();
    }
}
