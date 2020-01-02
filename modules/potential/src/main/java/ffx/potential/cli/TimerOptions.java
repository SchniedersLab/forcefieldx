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
package ffx.potential.cli;

import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that perform timings for energy and optionally gradients.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class TimerOptions {
    /**
     * -n or --iterations to set the number of iterations
     */
    @Option(names = {"-n", "--iterations"}, paramLabel = "5",
            description = "Number of iterations.")
    private int iterations = 5;

    /**
     * --nt or --threads to set the number of SMP threads (the default of 0 specifies use of all CPU cores)
     */
    @Option(names = {"--nt", "--threads"}, paramLabel = "0",
            description = "Number of SMP threads (0 specifies use of all CPU cores).")
    private int threads = 0;

    /**
     * -g or --noGradient to ignore computation of the atomic coordinates noGradient
     */
    @Option(names = {"-g", "--noGradient"},
            description = "Ignore computation of the atomic coordinates noGradient.")
    private boolean noGradient = false;

    /**
     * -v or --verbose to suppress printing of the energy for each iteration
     */
    @Option(names = {"-v", "--verbose"},
            description = "Print the energy for each iteration.")
    private boolean verbose = false;

    /**
     * <p>Getter for the field <code>iterations</code>.</p>
     *
     * @return a int.
     */
    public int getIterations() {
        return iterations;
    }

    /**
     * <p>Getter for the field <code>threads</code>.</p>
     *
     * @return a int.
     */
    public int getThreads() {
        return threads;
    }

    /**
     * <p>Getter for the field <code>noGradient</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getNoGradient() {
        return noGradient;
    }

    /**
     * <p>Getter for the field <code>verbose</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getVerbose() {
        return verbose;
    }
}
