/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

import ffx.utilities.LoggerSevereError;
import java.awt.GraphicsEnvironment;
import java.util.logging.ErrorManager;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import org.apache.commons.lang3.exception.ExceptionUtils;

/**
 * The default ConsoleHanlder publishes logging to System.err. This class
 * publishes to System.out, which is normally intercepted by the Force Field X
 * Shell. The formatter used reduces verbosity relative to the default
 * SimpleFormatter.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class LogHandler extends Handler {

    private static final boolean headless = GraphicsEnvironment.isHeadless();
    private static final boolean tryCatchSevere;
    private MainPanel mainPanel = null;
    private boolean fatal = false;

    static {
        String tryCatchSevereStr = System.getProperty("tryCatchSevere");
        if (tryCatchSevereStr != null) {
            tryCatchSevere = Boolean.parseBoolean(tryCatchSevereStr);
        } else {
            tryCatchSevere = false;
        }
    }
    
    /**
     * A reference to the Force Field X MainPanel container to shut down if we
     * encounter a fatal (SEVERE) exception. If we are not in Headless mode,
     * then LogRecords can be published to the ModelingShell.
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
     *
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
             * We don't want to throw an exception here, but we report the
             * exception to any registered ErrorManager.
             */
            reportError(null, e, ErrorManager.FORMAT_FAILURE);
            return;
        }
        try {
            if (record.getLevel() == Level.SEVERE) {
                fatal = true;
                System.err.println(msg);
                
                Throwable throwable = record.getThrown();
                if (throwable != null) {
                    System.err.println(String.format(" Exception %s logged.", throwable));
                }
                
                // If tryCatchSevere, and the throwable (if it exists) is not an Error, then...
                if (tryCatchSevere && (throwable == null || !(throwable instanceof Error))) {
                    System.err.println(" Force Field X may not continue.");
                    System.err.println(" Throwing new error...");
                    fatal = false;
                    if (throwable != null) {
                        throw new LoggerSevereError(throwable);
                    } else {
                        throw new LoggerSevereError(" Unknown exception");
                    }
                }
                
                System.err.println(" Force Field X will not continue.");
                System.err.println(" Shutting down...");
                flush();
                mainPanel.exit();
            }
            ModelingShell shell = null;
            if (mainPanel != null) {
                shell = mainPanel.getModelingShell();
            }

            if (!headless && shell != null) {
                shell.appendOutputNl(msg, shell.getResultStyle());
            } else {
                System.out.println(msg);
            }
        } catch (Exception e) {
            /**
             * We don't want to throw an exception here, but we report the
             * exception to any registered ErrorManager.
             */
            reportError(null, e, ErrorManager.WRITE_FAILURE);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void flush() {
        System.out.flush();
        System.err.flush();
        if (mainPanel.getModelingShell() != null) {
            // Scroll to visible!
        }
    }

    /**
     * {@inheritDoc}
     *
     * Flush, but do not close System.out or the Shell.
     *
     * @since 1.0
     */
    @Override
    public void close() {
        flush();
    }
}
