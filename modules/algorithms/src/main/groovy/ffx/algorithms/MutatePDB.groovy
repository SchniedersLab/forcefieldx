package ffx.algorithms

import org.apache.commons.configuration.CompositeConfiguration

import ffx.algorithms.cli.AlgorithmsScript
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.ForceFieldFilter
import ffx.potential.parsers.PDBFilter
import ffx.utilities.Keyword

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The MutatePDB script mutates a residue of a PDB file.
 * <br>
 * Usage:
 * <br>
 * ffxc MutatePDB [options] &lt;pdb&gt;
 */
@Command(description = " Mutate a PDB residue.", name = "ffxc MutatePDB")
class MutatePDB extends AlgorithmsScript {

    /**
     * -r or --resid Residue number.
     */
    @Option(names = ['--resid', '-r'], paramLabel = '1',
            description = 'Residue number.')
    int resID = 1
    /**
     * -n or --resname New residue name.
     */
    @Option(names = ['--resname', '-n'], paramLabel = 'ALA',
            description = 'New residue name.')
    String resName = "ALA"
    /**
     * -c or --chain Single character chain name (default is ' ').
     */
    @Option(names = ['--chain', '-c'], paramLabel = ' ',
            description = 'Single character chain name (default is \' \').')
    Character chain = ' '
    /**
     * -p or --repack After mutation, repack all residues within the specified cutoff radius (Angstroms).
     */
    @Option(names = ['--repack', '-p'], paramLabel = '-1.0',
            description = 'After mutation, repack all residues within a cutoff radius (Angstroms).')
    double repackDistance = -1.0
    /**
     * -pt or --twoBodyRepack Do not include three-body energies in repacking.
     */
    @Option(names = ['--threeBody', '--tB'], paramLabel = 'false',
            description = 'Include three-body energies in repacking.')
    boolean threeBody = false

    /**
     * -eR or --energyRestart Load energy restart file from a previous run (ensure that all parameters are the same).
     */
    @Option(names = ['--energyRestart', '--eR'], paramLabel = 'filename',
            description = 'Load energy restart file from a previous run (ensure that all parameters are the same).')
    String energyRestart = null
    /**
     * -R or --rotamer Rotamer number to apply.
     */
    @Option(names = ['--rotamer', '-R'], paramLabel = '-1', description = 'Rotamer number to apply.')
    int rotamer = -1

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files", description = "A PDB input files.")
    private List<String> filenames

    /**
     * Execute the script.
     */
    def run() {

        if (!init()) {
            return
        }

        boolean repack = false;
        if (repackDistance > -1) {
            repack = true
        }

        boolean threeBodyRepack = !twoBodyRepack
        boolean useEnergyRestart = false
        File energyRestartFile = null
        if (energyRestart != null) {
            useEnergyRestart = true
            energyRestartFile = new File(energyRestart)
        }

        int destRotamer = 0
        if (rotamer > -1) {
            if (repack) {
                logger.severe(" Can't combine repack with explicit rotamer specification.")
            }
            destRotamer = rotamer
        }
        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary()

        logger.info("\n Mutating residue number " + resID + " of chain " + chain + " to " + resName)

        // Read in command line.
        String filename = filenames.get(0)
        File structure = new File(filename)
        int index = filename.lastIndexOf(".")
        String name = filename.substring(0, index)
        MolecularAssembly molecularAssembly = new MolecularAssembly(name)
        molecularAssembly.setFile(structure)

        CompositeConfiguration properties = Keyword.loadProperties(structure)
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties)
        ForceField forceField = forceFieldFilter.parse()
        molecularAssembly.setForceField(forceField)

        PDBFilter pdbFilter = new PDBFilter(structure, molecularAssembly, forceField, properties)
        pdbFilter.mutate(chain, resID, resName)
        pdbFilter.readFile()
        pdbFilter.applyAtomProperties()
        molecularAssembly.finalize(true, forceField)

        if (repack) {
            logger.info("\n Repacking... \n")
            ForceFieldEnergy forceFieldEnergy = ForceFieldEnergy.energyFactory(molecularAssembly)
            molecularAssembly.setPotential(forceFieldEnergy)
            forceFieldEnergy.setPrintOnFailure(false, false)

            // Do a sliding-window rotamer optimization on a single one-residue window with a radius-inclusion criterion.
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson)
            rLib.setUseOrigCoordsRotamer(true)

            // This does break encapsulation of our modules.
            RotamerOptimization rotamerOptimization = new RotamerOptimization(molecularAssembly, forceFieldEnergy, algorithmListener)
            rotamerOptimization.setThreeBodyEnergy(threeBody)
            rotamerOptimization.setForcedResidues(resID, resID)
            rotamerOptimization.setWindowSize(1)
            rotamerOptimization.setDistanceCutoff(repackDistance)
            if (useEnergyRestart) {
                rotamerOptimization.setEnergyRestartFile(energyRestartFile)
            }

            rotamerOptimization.setResidues(chain, resID, resID)
            ArrayList<Residue> residueList = rotamerOptimization.getResidues()

            algorithmFunctions.energy(molecularAssembly)

            RotamerLibrary.measureRotamers(residueList, false)
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW)
            logger.info("\n Repacking successful.\n")
        }

        if (destRotamer > -1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson)
            Polymer polymer = molecularAssembly.getChain(chain.toString())
            Residue residue = polymer.getResidue(resID)
            Rotamer[] rotamers = residue.getRotamers(rLib)
            if (rotamers != null && rotamers.length > 0) {
                RotamerLibrary.applyRotamer(residue, rotamers[destRotamer])
            } else {
                logger.info(" No rotamer to apply.")
            }
        }
        pdbFilter.writeFile(structure, false)
    }
}

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
