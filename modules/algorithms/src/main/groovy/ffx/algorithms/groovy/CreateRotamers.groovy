//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.algorithms.optimize.Minimize
import ffx.algorithms.optimize.Minimize.MinimizationEngine
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.*
import ffx.potential.bonded.RotamerLibrary.ProteinLibrary
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The CreateRotamers script creates a set of conformation dependent rotamers.
 * <br>
 * Usage:
 * <br>
 * ffxc CreateRotamers [options] &lt;filename&gt;
 */
@Command(description = " Creates a set of conformation dependent rotamers.", name = "ffxc CreateRotamers")
class CreateRotamers extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    // TODO: instead @Mixin a subset of current ManyBodyOptions.

    /**
     * -L or --library Choose either Ponder and Richards (1) or Richardson (2)
     * rotamer library.
     */
    @Option(names = ["-L", "--library"], paramLabel = "2",
            description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    int library = 2

    /**
     * The final argument should be a filename.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    public ForceFieldEnergy forceFieldEnergy = null
    public ForceFieldEnergy secondForceFieldEnergy = null

    @Override
    CreateRotamers run() {

        if (!init()) {
            return this
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info(" Running CreateRotamers on " + filename)

        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        // Set all atoms to be "inactive".
        for (int i = 0; i < nAtoms; i++) {
            atoms[i].setActive(false)
        }

        // For now, always use the original coordinates as a (fixed) rotamer.
        boolean useOriginalRotamers = true

        // AA Library
        RotamerLibrary rotamerLibrary = new RotamerLibrary(ProteinLibrary.intToProteinLibrary(library), useOriginalRotamers)

        // Initialize Default NA Coordinates
        Polymer[] polymers = activeAssembly.getChains()
        RotamerLibrary.initializeDefaultAtomicCoordinates(polymers)

        // Get the residue list.
        List<Residue> residues = activeAssembly.getResidueList().stream().
                filter({
                    Residue r -> Rotamer[] rots = r.getRotamers(rotamerLibrary)
                    return rots != null && rots.length > 1 }).collect(Collectors.toList())

        logger.info(String.format(" Number of residues: %d\n", residues.size()))

        // Loop over Residues and set side chain atoms to not be used.
        for (Residue residue : residues) {
            for (Atom atom : residue.getVariableAtoms()) {
                atom.setUse(false)
            }
        }

        // Create .rot file name: should match input file name and end in ".rot"
        String rotFileName = String.format("%s.rot", FilenameUtils.removeExtension(filename))

        BufferedWriter bw = null
        try {
            bw = new BufferedWriter(new FileWriter(new File(rotFileName)))

            // TODO: Make this ALGORITHM:[ALGORITHM]:[box/window number] instead of assuming global:1.
            bw.write("ALGORITHM:GLOBAL:1")
            bw.newLine();

            // Loop over Residues
            for (Residue residue : residues) {

                StringBuilder resLine = new StringBuilder(" RES:")
                resLine.append(residue.getChainID()).append(":")
                resLine.append(residue.getSegID()).append(":")
                resLine.append(residue.getName()).append(":")
                resLine.append(residue.getResidueNumber()).append("\n")
                bw.write(resLine.toString())

                // Get this residue's rotamers.
                Rotamer[] rotamers = residue.getRotamers(rotamerLibrary)

                assert rotamers != null && rotamers.length > 1

                // Configure "active" and "use" flags.
                // .getVariableAtoms returns backbone atoms for
                // nucleic acids, which is what we want
                List<Atom> sideChainAtoms = residue.getVariableAtoms()
                for (Atom atom : sideChainAtoms) {
                    atom.setActive(true)
                    atom.setUse(true)
                }

                // Define "all previously saved rotamers" arrayList to be used for RMSD comparison
                ArrayList<MolecularAssembly> keptRotamers = new ArrayList<>()
                // Loop over rotamers for this Residue.
                for (int i = 0; i < rotamers.length; i++) {
                    Rotamer rotamer = rotamers[i]

                    // Apply the rotamer (i.e. amino acid side-chain or nucleic acid suite).
                    RotamerLibrary.applyRotamer(residue, rotamer)

                    bw.write(String.format("  ROT:%d\n", i))

                    if (i > 0 || !useOriginalRotamers) {
                        // -Dplatform=omm
                        MinimizationEngine engine = Minimize.defaultEngine(activeAssembly.getPotentialEnergy())
                        Minimize minimize = Minimize.minimizeFactory(activeAssembly,
                                activeAssembly.getPotentialEnergy(), algorithmListener, engine)
                        // Locally minimize.
                        minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())
                    } else {
                        logger.info(" Skipping minimization of original-coordinates rotamer.")
                    }

                    // Only save if the "new" rotamer isn't within 0.1 kcal/mol of any other rotamer
                    // Update "all previously saved rotamers" arrayList
                    if(i == 0){
                        // Add 0th rotamer to the assemblies arrayList
                        keptRotamers.add(activeAssembly)

                        // Save out coordinates to a rotamer file.
                        for (Atom atom : sideChainAtoms) {
                            double x = atom.getX()
                            double y = atom.getY()
                            double z = atom.getZ()
                            logger.info(String.format(" %s %16.8f %16.8f %16.8f", atom.toString(), x, y, z))
                            StringBuilder atomLine = new StringBuilder("   ATOM:")
                            atomLine.append(atom.getName()).append(":")
                            atomLine.append(x).append(":")
                            atomLine.append(y).append(":")
                            atomLine.append(z).append("\n")
                            bw.write(atomLine.toString())
                        }
                        bw.write("  ENDROT\n")
                    } else{
                        // For all but the 0th rotamer, do RMSD calculations to determine if the "new" rotamer
                        // is within 0.1 kcal/mol of any other rotamer. Only save it if it's not.
                        //MolecularAssembly newRotamer = activeAssembly
                        println("Rotamer number: "+i)

                        println("Number of rotamers in keptRotamers: "+keptRotamers.size())

                        forceFieldEnergy = activeAssembly.getPotentialEnergy()
                        Atom[] rotamerAtoms = activeAssembly.getAtomArray()
                        println("rotamerAtoms: "+rotamerAtoms.length)
                        //for (int atmcount = 0; atmcount < rotamerAtoms.length; atmcount++){
                        //    print(rotamerAtoms[atmcount].name+":")
                        //}
                        println()
                        int nVars = forceFieldEnergy.getNumberOfVariables()
                        println("nVars is: "+nVars)
                        double[] x1 = new double[nVars]
                        forceFieldEnergy.getCoordinates(x1)

                        double[] x2 = new double[nVars]
                        double[] mass = new double[rotamerAtoms.length]
                        //println("Mass array: "+mass.length)

                        int nRotamerAtoms = rotamerAtoms.length
                        for (int j = 0; j < nRotamerAtoms; j++) {
                            mass[j] = rotamerAtoms[j].getMass()
                        }

                        // Indices of atoms used in alignment and RMSD calculations.
//                        int[] usedIndices = atomIndexStream.toArray()
//                        int nUsed = usedIndices.length
//                        int nUsedVars = nUsed * 3
//                        double[] massUsed = Arrays.stream(usedIndices).
//                                mapToDouble({ int j -> rotamerAtoms[j].getAtomType().atomicWeight }).toArray()
//                        double[] xUsed = new double[nUsedVars]
//                        double[] x2Used = new double[nUsedVars]

                        // Superpose.rmsd(systemFilter, nUsed, usedIndices, x, x2, xUsed, x2Used, massUsed, 1)
                        // Define empty array to store RMSDs and RMSD threshold value boolean
                        boolean withinRange = false
                        for (int k = 0; k < keptRotamers.size(); k++){
                            secondForceFieldEnergy = keptRotamers[k].getPotentialEnergy()
                            secondForceFieldEnergy.getCoordinates(x2)
                            for(int xc1 = 0; xc1< x1.length;xc1++){
                                print(x1[xc1]+":")
                            }
                            println()
                            println()
                            for(int xc = 0; xc < x2.length; xc++){
                                print(x2[xc]+":")
                            }
                            println()
                            double origRMSD = ffx.potential.utils.Superpose.rmsd(x1, x2, mass)
                            println("RMSD: "+origRMSD)
                            if(origRMSD <= 0.01){withinRange = true}
                        }

                        // If the new rotamer is within 0.1 kcal/mol of any previously kept rotamer
                        // (i.e.: if any of the values in "rmsdValues" are 0.1 or less),
                        // don't keep that rotamer.  Keep it otherwise
                        if(withinRange){
                            // Rotamer not written because it's too energetically similar to another rotamer in the set
                            logger.info("Rotamer not kept")
                        } else{
                            // Save out coordinates to a rotamer file.
                            for (Atom atom : sideChainAtoms) {
                                double x = atom.getX()
                                double y = atom.getY()
                                double z = atom.getZ()
                                logger.info(String.format(" %s %16.8f %16.8f %16.8f", atom.toString(), x, y, z))
                                StringBuilder atomLine = new StringBuilder("   ATOM:")
                                atomLine.append(atom.getName()).append(":")
                                atomLine.append(x).append(":")
                                atomLine.append(y).append(":")
                                atomLine.append(z).append("\n")
                                bw.write(atomLine.toString())
                            }
                            bw.write("  ENDROT\n")

                            // Add it to keptRotamers list
                            keptRotamers.add(activeAssembly)
                        }

                    }
                    // Original location of "Save out coordinates to a rotamer file" logic
                }

                // Set the Residue conformation back to rotamer 0.
                RotamerLibrary.applyRotamer(residue, rotamers[0])

                // Revert the active and use flags.
                for (Atom atom : sideChainAtoms) {
                    atom.setActive(false)
                    atom.setUse(false)
                }
            }

        } finally {
            bw?.flush()
            bw?.close()
        }

        return this
    }
}