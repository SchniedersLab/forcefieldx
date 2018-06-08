/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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

package ffx.utilities

import groovy.cli.picocli.CliBuilder

boolean testMode = false;
Optional<Double> psPerHisto = Optional.empty();
Optional<Integer> blockSizeStep = Optional.empty();
Optional<Integer> maxBlockSize = Optional.empty();
Optional<String> grepCmd = Optional.empty();

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc utilities.blockAverage [options] <logFile>');
cli.h(longOpt:'help', 'Print this help message.');
cli.dt(longOpt:'psPerHisto', args:1, argName:'1.0', 'Number of picoseconds between each input histogram.');
cli.g(longOpt:'grepCmd', args:1, argName:'grep', 'Location of grep executable on non-UNIX systems.');
cli.m(longOpt:'maxBlockSize', args:1, argName:'-1', 'Maximum block size to attempt, else uses numHistograms.');
cli.s(longOpt:'blockSizeStep', args:1, argName:'100', 'Step increment for block size.');
cli.gt(longOpt:'generateTestData', args:1, argName:'1000', 'Create correlated and uncorrelated validation sets.');
cli.t(longOpt:'operateTestData', 'Operate on headerless, two-column data.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);

if (options.gt) {
    int size = Integer.parseInt(options.gt);
    BlockAverager.generateTestData(filename, size);
    return;
}

if (options.dt) {
    psPerHisto = Optional.of(Double.parseDouble(options.t));
}
if (options.s) {
    blockSizeStep = Optional.of(Integer.parseInt(options.s));
}
if (options.m) {
    maxBlockSize = Optional.of(Integer.parseInt(options.m));
}
if (options.g) {
    grepCmd = Optional.of(options.g);
}
if (options.t) {
    testMode = true;
}

BlockAverager ba = new BlockAverager(filename, testMode, grepCmd, psPerHisto, blockSizeStep, maxBlockSize);
double[] binStdErrors = ba.computeBinUncertainties();
double totalStdError = ba.computeTotalUncertainty();
