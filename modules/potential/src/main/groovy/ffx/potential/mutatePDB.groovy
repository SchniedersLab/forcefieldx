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

package ffx.potentials;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.RotamerOptimization;
import ffx.algorithms.RotamerOptimization.Algorithm;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.RotamerLibrary.ProteinLibrary;

import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.PotentialsUtils;
import ffx.potential.ForceFieldEnergy;
import ffx.utilities.Keyword;

// Groovy Imports
import groovy.util.CliBuilder;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc mutatePDB [options] <PDB>');
cli.h(longOpt:'help', 'Print this help message.');
cli.r(longOpt:'resid', args:1, argName:'1', 'Residue number.');
cli.n(longOpt:'resname', args:1, argName:'ALA', 'New residue name.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name (default is \' \').');
cli.p(longOpt:'repack', args:1, argName:'7.0', 'After mutation, repack all residues within # Angstroms.');
cli.pt(longOpt:'threeBodyRepack', args:1, argName:'true', 'Include three-body energies in repacking.');

boolean repack = false;
double repackDistance = 7.0;
boolean threeBodyRepack = true;

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

int resID = 1;
String resName = "ALA";
Character chain = ' ';

// Residue number.
if (options.r) {
    resID = Integer.parseInt(options.r);
}

// Residue Name.
if (options.n) {
    resName = options.n;
}

// Chain Name.
if (options.c) {
    chain = options.c.toCharacter();
}

if (options.p) {
    repack = true;
    repackDistance = Double.parseDouble(options.p);
}

if (options.pt) {
    threeBodyRepack = Boolean.parseBoolean(options.pt);
}

// Read in command line.
String filename = arguments.get(0);

logger.info("\n Mutating residue number " + resID + " of chain " + chain + " to " + resName);

File structure = new File(filename);
int index = filename.lastIndexOf(".");
String name = filename.substring(0, index);
MolecularAssembly molecularAssembly = new MolecularAssembly(name);
molecularAssembly.setFile(structure);

CompositeConfiguration properties = Keyword.loadProperties(structure);
ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
ForceField forceField = forceFieldFilter.parse();
molecularAssembly.setForceField(forceField);

PDBFilter pdbFilter = new PDBFilter(structure, molecularAssembly, forceField, properties);
pdbFilter.mutate(chain,resID,resName);
pdbFilter.readFile();
molecularAssembly.finalize(true, forceField);

if (repack) {
    logger.info("\n Repacking... \n");
    ForceFieldEnergy forceFieldEnergy = new ForceFieldEnergy(molecularAssembly);
    molecularAssembly.setPotential(forceFieldEnergy);
    
    // Do a sliding-window rotamer optimization on a single one-residue window with a radius-inclusion criterion.
    RotamerLibrary.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
    RotamerLibrary.setUseOrigCoordsRotamer(true);
    
    RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
    rotamerOptimization.setThreeBodyEnergy(threeBodyRepack);
    rotamerOptimization.setWindowSize(1);
    rotamerOptimization.setDistanceCutoff(repackDistance);
    
    startResID = resID;
    finalResID = resID;
    if (options.c) {
        rotamerOptimization.setResidues(options.c, startResID, finalResID);
    } else {
        rotamerOptimization.setResidues(startResID, finalResID);
    }    
    
    residueList = rotamerOptimization.getResidues();
    energy();
    RotamerLibrary.measureRotamers(residueList, false);
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.SLIDING_WINDOW);
    logger.info("\n Repacking successful.\n");
}

pdbFilter.writeFile(structure, false);
