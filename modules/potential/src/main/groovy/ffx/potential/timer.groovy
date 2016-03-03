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

// TIMER

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.potential.ForceFieldEnergy;

// The number of iterations.
int nEvals = 5;

// Compute the atomic coordinate gradient.
boolean gradient = true;

// Print the energy for each iteraction.
boolean print = true;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc timer [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.n(longOpt:'iterations', args:1, argName:'5', 'Number of iterations');
cli.c(longOpt:'threads', args:1, argName:'all', 'Number of SMP threads (ie. default uses all CPU cores)');
cli.g(longOpt:'gradient', args:1, argName:'true', 'Compute the atomic coordinats gradeint');
cli.v(longOpt:'verbose', args:1, argName:'true', 'Print out the energy for each step');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Load the number iterations.
if (options.n) {
    nEvals = Integer.parseInt(options.n);
}

// Load the number of threads.
if (options.c) {
    System.setProperty("pj.nt", options.c);
}

// Compute the gradient for each step.
if (options.g) {
    gradient = Boolean.parseBoolean(options.g);
}

// Print the energy for each step.
if (options.v) {
    print = Boolean.parseBoolean(options.v);
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Timing energy and gradient for " + filename);

open(filename);

ForceFieldEnergy energy = active.getPotentialEnergy();

long minTime = Long.MAX_VALUE;
double sumTime2 = 0.0;
int halfnEvals = (nEvals % 2 == 1) ? (nEvals/2) : (nEvals/2) - 1; // Halfway point
for (int i=0; i<nEvals; i++) {
    long time = -System.nanoTime();
    energy.energy(gradient, print);
    time += System.nanoTime();
    minTime = time < minTime ? time : minTime;
    if (i >= (int) (nEvals/2)) {
        double time2 = time * 1.0E-9;
        sumTime2 += (time2*time2);
    }
}
++halfnEvals;
double rmsTime = Math.sqrt(sumTime2/halfnEvals);
logger.info(String.format(" Minimum time: %14.5f (sec)", minTime * 1.0E-9));
logger.info(String.format(" RMS time (latter half): %14.5f (sec)", rmsTime));
