/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

// SAVE AS PDB

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;
import ffx.potential.utils.PotentialsUtils;

String siftFileName = null;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc saveAsSIFTPDB [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.f(longOpt: 'fileName', args:1, argName:'GENENAME.out', 'File of sift scores to enter.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}

if (options.f) {
    siftFileName = options.f;
}

String[] data;

File siftFile = new File(siftFileName);

if (siftFile.exists()) {
    data = siftFile.text.split('\n')
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Writing out PDB for " + filename);

systems = open(filename);

PotentialsUtils potentialUtils = new PotentialsUtils();

filename = FilenameUtils.removeExtension(filename) + ".pdb";
if (siftFileName == null) {
    saveAsPDB(systems, new File(filename));
} else {
    potentialUtils.saveAsSIFTPDB(systems, new File(filename), data)
}
