
package ffx.potentials;

// Groovy Imports
import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.util.CliBuilder

import org.apache.commons.configuration.CompositeConfiguration;

// Force Field X Imports
import ffx.algorithms.RotamerOptimization;
import ffx.algorithms.RotamerOptimization.Algorithm;
import ffx.potential.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.RotamerLibrary.ProteinLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;

/**
 * The MutatePDB script mutates a residue of a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc MutatePDB [options] &lt;pdb&gt;
 */
class MutatePDB extends Script {

    /**
     * Options for the MutatePDB script.
     * <br>
     * Usage:
     * <br>
     * ffxc MutatePDB [options] &lt;filename&gt;
     */
    public class Options {

        /**
         * -h or --help to print a help message
         */
        @Option(longName='help', shortName='h', defaultValue='false', description='Print this help message.') boolean help
        /**
         * -r or --resid Residue number.
         */
        @Option(longName='resid', shortName='r', defaultValue='1', description='Residue number.') int resid
        /**
         * -n or --resname New residue name.
         */
        @Option(longName='resname', shortName='n', defaultValue='ALA', description='New residue name.') String resname
        /**
         * -c or --chain Single character chain name (default is ' ').
         */
        @Option(longName='chain', shortName='c', defaultValue=' ', description='Single character chain name (default is \' \').') Character chain
        /**
         * -p or --repack After mutation, repack all residues within the specified cutoff radius (Angstroms).
         */
        @Option(longName='repack', shortName='p', defaultValue='-1.0',
            description='After mutation, repack all residues within a cutoff radius (Angstroms).') double repack
        /**
         * -pt or --twoBodyRepack Do not include three-body energies in repacking.
         */
        @Option(longName='twoBodyRepack', shortName='pt', defaultValue='false',
            description='Include three-body energies in repacking.') boolean twoBodyRepack
        /**
         * -eR or --energyRestart Load energy restart file from a previous run (ensure that all parameters are the same).
         */
        @Option(longName='energyRestart', shortName='eR', defaultValue='filename',
            description='Load energy restart file from a previous run (ensure that all parameters are the same).') String energyRestart
        /**
         * -R or --rotamer Rotamer number to apply.
         */
        @Option(longName='rotamer', shortName='R', defaultValue='-1', description='Rotamer number to apply.') int rotamer
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }

    /**
     * Execute the script.
     */
    def run() {

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc MutatePDB [options] <PDB>')
        def options = new Options()
        cli.parseFromInstance(options, args)
        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filenames
        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
            //open(modelFilename)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }


        int resID = options.resid;
        String resName = options.resname
        Character chain = options.chain

        boolean repack = false;
        double repackDistance = options.repack
        if (repackDistance > -1) {
            repack = true
        }
        boolean threeBodyRepack = !options.twoBodyRepack
        boolean useEnergyRestart = false
        energyRestartFile = null
        if (!options.energyRestart.equalsIgnoreCase('filename')) {
            useEnergyRestart = true
            energyRestartFile = new File(options.energyRestart)
        }

        int destRotamer = 0;
        if (options.rotamer > -1) {
            if (repack) {
                logger.severe(" Can't combine repack with explicit rotamer specification.");
            }
            destRotamer = options.rotamer
        }
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary()

        logger.info("\n Mutating residue number " + resID + " of chain " + chain + " to " + resName);

        // Read in command line.
        String filename = arguments.get(0);
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
        pdbFilter.applyAtomProperties();
        molecularAssembly.finalize(true, forceField);

        if (repack) {
            logger.info("\n Repacking... \n");
            ForceFieldEnergy forceFieldEnergy = new ForceFieldEnergy(molecularAssembly);
            forcefieldEnergy.setPrintOnFailure(false, false);
            molecularAssembly.setPotential(forceFieldEnergy);

            // Do a sliding-window rotamer optimization on a single one-residue window with a radius-inclusion criterion.
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
            rLib.setUseOrigCoordsRotamer(true);

            RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, null);
            rotamerOptimization.setThreeBodyEnergy(threeBodyRepack);
            rotamerOptimization.setForcedResidues(resID, resID);
            rotamerOptimization.setWindowSize(1);
            rotamerOptimization.setDistanceCutoff(repackDistance);
            if (useEnergyRestart) {
                rotamerOptimization.setEnergyRestartFile(energyRestartFile);
            }

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

        if (destRotamer > -1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);
            Polymer polymer = molecularAssembly.getChain(chain.toString());
            Residue residue = polymer.getResidue(resID);
            Rotamer[] rotamers = residue.getRotamers(rLib);
            RotamerLibrary.applyRotamer(residue, rotamers[destRotamer]);
        }
        pdbFilter.writeFile(structure, false);
    }
}

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