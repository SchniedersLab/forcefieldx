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

// Convert from PRM to Property

// Apache Imports
import org.apache.commons.configuration.CompositeConfiguration;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.utilities.Keyword;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc prmToProperty <prm> [prm] ...');
cli.h(longOpt:'help', 'Print this help message.');
cli.i(longOpt:'improper', args:1, argName:'1.0', 'Scale Improper Torsions');
cli.r(longOpt:'radii', args:1, argName:'false', 'Convert vdW Radii values to Diameters');
cli.s(longOpt:'sigma', args:1, argName:'false', 'Convert vdW Sigma values to R-Min');
cli.t(longOpt:'torsion', args:1, argName:'1.0', 'Scale Torsions');

def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Read in the command line file.
String xyzname = arguments.get(0);
CompositeConfiguration properties = Keyword.loadProperties(null);
properties.setProperty("parameters", xyzname);
ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);

if (options.i) {
    improperTorsionScale = Double.parseDouble(options.i);
    forceFieldFilter.setImproperTorsionScale(improperTorsionScale);
}

if (options.r) {
    convertRadiusToDiameter = Boolean.parseBoolean(options.r);
    forceFieldFilter.setConvertRadiusToDiameter(convertRadiusToDiameter);
}

if (options.s) {
    convertSigmaToRMin = Boolean.parseBoolean(options.s);
    forceFieldFilter.setConvertSigmaToRMin(convertSigmaToRMin);
}

if (options.t) {
    torsionScale = Double.parseDouble(options.t);
    forceFieldFilter.setTorsionScale(torsionScale);
}

ForceField forceField = forceFieldFilter.parse();

int prms = arguments.size();
for (int i=1; i<prms; i++) {
    xyzname = arguments.get(i);
    properties = Keyword.loadProperties(null);
    properties.setProperty("parameters", xyzname);
    forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField2 = forceFieldFilter.parse();
    forceField.append(forceField2);
}

if (forceField != null) {
    forceField.print();
}
