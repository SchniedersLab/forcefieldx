//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
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

    private static final int warningLevel = Level.WARNING.intValue();
    private final boolean debug;
    private final boolean mpiLogging;

    /**
     * Constructor for the LogFormatter.
     *
     * @param debug      If debug is true, then LogFormatter is equivalent to {@link SimpleFormatter}.
     * @param mpiLogging Configure for MPI logging.
     * @since 1.0
     */
    LogFormatter(boolean debug, boolean mpiLogging) {
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
        String message;
        if (debug || record.getLevel().intValue() >= warningLevel) {
            message = super.format(record);
        } else {
            message = record.getMessage();
            Object[] objects = record.getParameters();
            message = MessageFormat.format(message, objects);
        }

        try {
            Comm comm = Comm.world();
            int size = comm.size();
            int rank = comm.rank();
            if (size > 1) {
                if (mpiLogging) {
                    String[] lines = message.split("\n");
                    message = mpiFormat(size, rank, lines);
                } else if (rank != 0) {
                    message = null;
                }
            }
        } catch (Exception e) {
            // If Comm.world does not exist, do not append the rank.
        }

        return message;
    }

    /**
     * Prepend the MPI rank to a line of text to give: " [Rank]line"
     *
     * @param size  Number of MPI processes.
     * @param rank  Rank of this MPI processs.
     * @param lines The String to format split by new line characters.
     * @return " [Rank]line"
     */
    private String mpiFormat(int size, int rank, String[] lines) {
        int rankLen = Integer.toString(size - 1).length();
        String formatString = " [%0" + rankLen + "d]%s";
        StringBuilder sb = new StringBuilder(String.format(formatString, rank, lines[0]));
        for (int i = 1; i < lines.length; i++) {
            sb.append("\n");
            sb.append(String.format(formatString, rank, lines[i]));
        }
        return sb.toString();
    }
}
