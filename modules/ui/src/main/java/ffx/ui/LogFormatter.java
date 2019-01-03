/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import java.text.MessageFormat;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;

import edu.rit.pj.Comm;

/**
 * A minor extension to the SimpleFormatter to reduce verbosity if debugging is
 * not turned on.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class LogFormatter extends SimpleFormatter {

    private final boolean debug;
    private final boolean mpiLogging;
    private static final int warningLevel = Level.WARNING.intValue();

    /**
     * Constructor for the LogFormatter.
     *
     * @param debug If debug is true, then LogFormatter is equivalent to
     *              {@link SimpleFormatter}.
     * @param mpiLogging Configure for MPI logging.
     * @since 1.0
     */
    public LogFormatter(boolean debug, boolean mpiLogging) {

        this.debug = debug;
        this.mpiLogging = mpiLogging;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Unless debugging is turned on or the LogRecord is of level WARNING or
     * greater, just return the message.
     * <p>
     * If more than one process is active, prepend the rank of the process to
     * each line of the message.
     *
     * @since 1.0
     */
    @Override
    public String format(LogRecord record) {

        String message = null;

        if (debug || record.getLevel().intValue() >= warningLevel) {
            message = super.format(record);
        } else {
            message = record.getMessage();
            Object objects[] = record.getParameters();
            message = MessageFormat.format(message, objects);
        }

        try {
            Comm comm = Comm.world();
            int size = comm.size();
            int rank = comm.rank();
            if (size > 1) {
                if (mpiLogging) {
                    String lines[] = message.split("\n");
                    message = String.format(" [%d]%s", rank, lines[0]);
                    int numLines = lines.length;
                    for (int i = 1; i < numLines; i++) {
                        message = message.concat(String.format("\n [%d]%s", rank, lines[i]));
                    }
                } else if (rank != 0) {
                    message = null;
                }
            }
        } catch (Exception e) {
            // If Comm.world does not exist, do not append the rank.
        }

        return message;
    }
}
