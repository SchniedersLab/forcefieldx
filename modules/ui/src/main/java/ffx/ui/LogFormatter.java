/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.ui;

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
 *
 */
public class LogFormatter extends SimpleFormatter {

    private final boolean debug;
    private static final int warningLevel = Level.WARNING.intValue();

    /**
     * Constructor for the LogFormatter.
     *
     * @param debug If debug is true, then LogFormatter is equivalent to
     * {@link SimpleFormatter}.
     * @since 1.0
     */
    public LogFormatter(boolean debug) {
        this.debug = debug;
    }

    /**
     * {@inheritDoc}
     *
     * Unless debugging is turned on or the LogRecord is of level WARNING or
     * greater, just return the message.
     *
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
        }

        try {
            Comm comm = Comm.world();
            int size = comm.size();
            if (size > 1) {
                int rank = comm.rank();
                String lines[] = message.split("\n");
                message = String.format(" [%d]%s", rank, lines[0]);
                int numLines = lines.length;
                for (int i = 1; i < numLines; i++) {
                    message = message.concat(String.format("\n [%d]%s", rank, lines[i]));
                }
            }
        } catch (Exception e) {
            // If Comm.world does not exist, do not append the rank.
        }

        return message;
    }
}
