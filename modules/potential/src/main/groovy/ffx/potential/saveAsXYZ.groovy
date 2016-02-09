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

// SAVE AS XYZ

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// FFX Imports
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;

// Groovy Imports
import groovy.util.CliBuilder;

// Things below this line normally do not need to be changed.
// ===============================================================================================
int offset = 0;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc saveAsXYZ [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.p(longOpt:'pos-offset', args:1, argName:'0', 'Positive offset of atom types in the new file.');
cli.n(longOpt:'neg-offset', args:1, argName:'0', 'Negative offset of atom types in the new file.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Positive offset atom types.
if (options.p) {
    offset = Integer.parseInt(options.p);
}

// Offset atom types.
if (options.n) {
    offset = Integer.parseInt(options.n);
    offset = -offset;
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Writing out XYZ for " + filename);

open(filename);

// Offset atom type numbers.
if (offset != 0) {
    logger.info("\n Offset atom types by " + offset);
    MolecularAssembly molecularAssembly = active;
    ForceField forceField = molecularAssembly.getForceField();
    forceField.renumberForceField(0,offset,0);
}

filename = FilenameUtils.removeExtension(filename) + ".xyz";

save(new File(filename));
