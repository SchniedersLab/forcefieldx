package ffx.potential

import org.apache.commons.io.FilenameUtils

import ffx.crystal.Crystal
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Molecule
import ffx.potential.bonded.Polymer
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The MoveIntoUnitCell script moves the center of mass of each molecule into the unit cell.
 * <br>
 * Usage:
 * <br>
 * ffxc MoveIntoUnitCell &lt;filename&gt;
 */
@Command(description = " Move all molecules into the unit cell.", name = "ffxc MoveIntoUnitCell")
class MoveIntoUnitCell extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    /**
     * Execute the script.
     */
    def run() {

        if (!init()) {
            return
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            assemblies = [activeAssembly]
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Moving molecular centers of mass into the unit cell for " + modelFilename + "\n");

        // Loop over each system.
        for (int i = 0; i < assemblies.length; i++) {
            MolecularAssembly molecularAssembly = assemblies[i]
            Crystal crystal = molecularAssembly.getCrystal().getUnitCell()

            int nAtoms = 0
            double[] com = new double[3]
            double[] translate = new double[3]

            // Move the polymers together
            List<Polymer> polymers = molecularAssembly.getChains()
            if (polymers != null && polymers.size() > 0) {

                // Find the center of mass
                for (polymer in polymers) {
                    List<Atom> atoms = polymer.getAtomList()
                    nAtoms = nAtoms + atoms.size()
                    for (atom in atoms) {
                        com[0] = com[0] + atom.getX()
                        com[1] = com[1] + atom.getY()
                        com[2] = com[2] + atom.getZ()
                    }
                }
                com[0] = com[0] / nAtoms
                com[1] = com[1] / nAtoms
                com[2] = com[2] / nAtoms

                // Calculate the translation vector for the center of mass
                crystal.toPrimaryCell(com, translate)
                translate[0] = translate[0] - com[0]
                translate[1] = translate[1] - com[1]
                translate[2] = translate[2] - com[2]

                // Move each atom
                for (polymer in polymers) {
                    List<Atom> atoms = polymer.getAtomList()
                    for (atom in atoms) {
                        atom.move(translate)
                    }
                }
            }

            // Loop over each molecule
            List<Molecule> molecules = molecularAssembly.getMolecules()
            for (molecule in molecules) {
                List<Atom> atoms = molecule.getAtomList()
                // Find the center of mass
                com[0] = 0.0
                com[1] = 0.0
                com[2] = 0.0
                for (atom in atoms) {
                    com[0] = com[0] + atom.getX()
                    com[1] = com[1] + atom.getY()
                    com[2] = com[2] + atom.getZ()
                }

                nAtoms = atoms.size()
                com[0] = com[0] / nAtoms
                com[1] = com[1] / nAtoms
                com[2] = com[2] / nAtoms

                // Calculate the translation vector for the center of mass
                crystal.toPrimaryCell(com, translate)
                translate[0] = translate[0] - com[0]
                translate[1] = translate[1] - com[1]
                translate[2] = translate[2] - com[2]

                // Move each atom
                for (atom in atoms) {
                    atom.move(translate)
                }
            }

            // Loop over each water
            List<MSNode> waters = molecularAssembly.getWaters()
            for (water in waters) {
                List<Atom> atoms = water.getAtomList()
                // Find the center of mass
                com[0] = 0.0
                com[1] = 0.0
                com[2] = 0.0
                for (atom in atoms) {
                    com[0] = com[0] + atom.getX()
                    com[1] = com[1] + atom.getY()
                    com[2] = com[2] + atom.getZ()
                }

                nAtoms = atoms.size()
                com[0] = com[0] / nAtoms
                com[1] = com[1] / nAtoms
                com[2] = com[2] / nAtoms

                // Calculate the translation vector for the center of mass
                crystal.toPrimaryCell(com, translate)
                translate[0] = translate[0] - com[0]
                translate[1] = translate[1] - com[1]
                translate[2] = translate[2] - com[2]

                // Move each atom
                for (atom in atoms) {
                    atom.move(translate)
                }
            }

            // Loop over each ion
            List<MSNode> ions = molecularAssembly.getIons()
            for (ion in ions) {
                List<Atom> atoms = ion.getAtomList()
                // Find the center of mass
                com[0] = 0.0
                com[1] = 0.0
                com[2] = 0.0
                for (atom in atoms) {
                    com[0] = com[0] + atom.getX()
                    com[1] = com[1] + atom.getY()
                    com[2] = com[2] + atom.getZ()
                }

                nAtoms = atoms.size()
                com[0] = com[0] / nAtoms
                com[1] = com[1] / nAtoms
                com[2] = com[2] / nAtoms

                // Calculate the translation vector for the center of mass
                crystal.toPrimaryCell(com, translate)
                translate[0] = translate[0] - com[0]
                translate[1] = translate[1] - com[1]
                translate[2] = translate[2] - com[2]

                // Move each atom
                for (atom in atoms) {
                    atom.move(translate)
                }
            }
        }

        String ext = FilenameUtils.getExtension(modelFilename)
        modelFilename = FilenameUtils.removeExtension(modelFilename)
        if (ext.toUpperCase().contains("XYZ")) {
            potentialFunctions.saveAsXYZ(assemblies[0], new File(modelFilename + ".xyz"))
        } else {
            potentialFunctions.saveAsPDB(assemblies, new File(modelFilename + ".pdb"))
        }
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
