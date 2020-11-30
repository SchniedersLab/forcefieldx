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
package ffx.potential.groovy

import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.cli.PotentialScript
import ffx.potential.utils.PACCOMFunctions
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.Superpose
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.*

/**
 * The PACCOM script calculates intracrystollographic distances.
 *
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler, Kaleb Bierstedt, and Micheal Schnieders 2020
 * <br>
 * Usage:
 * <br>
 * ffxc test.PACCOM &lt;filename&gt &lt;filename&gt;
 */
@Command(description = " Compare crystal packings based on intermolecular distances.", name = "ffxc test.PACCOM")
class PACCOM extends PotentialScript {

    /**
     * --nm or --numberMolecules Number of molecules to include from each crystal in RMSD comparison.
     */
    @Option(names = ['--nm', '--numberMolecules'], paramLabel = '50',
            description = 'Determines crystal shell size for comparison.')
    int n_mol = 50

    /**
     * --na or --numberAlign Number of molecules to align between each crystal in RMSD comparison.
     */
    @Option(names = ['--na', '--numberAlign'], paramLabel = '3',
            description = 'Determines number of molecules to overlap during comparison.')
    int numMolAlign = 3

    /**
     * --ls or --largeShell Number of molecules in the large shell to compare with first crystal.
     */
    @Option(names = ['--ls', '--largeShell'], paramLabel = '500',
            description = 'Determines number of molecules in initial shell.')
    int largeShell = 500

    /**
     * --all or --allVsAll Compute RMSD of between each file.
     */
    @Option(names = ['--all', '--allVsAll'], paramLabel = "false", defaultValue = "false",
            description = 'Calculate the RMSD between each shell (instead of first shell)')
    private boolean allVsAll = false

    //TODO: Fix original PACCOM
    /**
     * --op or --originalPACCOM Use Algorithm supplied by original PACCOM.
     */
    @Option(names = ['--op', '--originalPACCOM'], paramLabel = "false", defaultValue = "false",
            description = 'Utilize a version of PACCOM that mirrors original.')
    private static boolean original = false

    /**
     * -s or --save Save out individual XYZ/PDB files for each created shell.
     */
    @Option(names = ['-s', '--save'], paramLabel = "false", defaultValue = "false",
            description = 'Save each shell as an aperiodic system to a PDB file.')
    private boolean saveFiles = false

    /**
     * --si or --saveIntermediates Save out intermediate crystals as PDB files.
     */
    @Option(names = ['--si', '--saveIntermediates'], paramLabel = "false", defaultValue = "false",
            description = 'Save each shell as an aperiodic system to a PDB file.')
    private static boolean saveIntermediates = false

    //TODO: implement closestDistance
    /**
     * --cd or --closestDistance Neighbor molecules calculated via closest atom (not currently implemented).
     */
    @Option(names = ['--cd', '--closestDistance'], paramLabel = "false", defaultValue = "false",
            description = 'Neighbors calculated via closest atom.')
    private static boolean closestDistance = false

    /**
     * Select individual atoms to match and calculate RMSD rather than using all atoms.
     */
    @Option(names = ['--ca', '--centerAtoms'], arity = "1..*", paramLabel = "atom integers",
            description = 'Specify atoms to match and calculate RMSD. Otherwise, all atom.')
    private int[] centerAtoms = null

    /**
     * The final argument(s) should be two or more filenames.
     */
    @Parameters(arity = "2..*", paramLabel = "files",
            description = 'Atomic coordinate files to compare in XYZ format.')
    List<String> filenames = null

    private MolecularAssembly[] assemblies

    public double[][] coordinates = null

    private File baseDir = null

    /**
     * Execute the script.
     */
    @Override
    PACCOM run() {
        System.setProperty("vdwterm", "false")
        // READ IN STRUCTURES
        if (!init()) {
            return
        }
        // Ensure file exists/create active molecular assembly
        assemblies = new MolecularAssembly[filenames.size()]
        if (filenames != null && filenames.size() > 0) {
            for (int i = 0; i < filenames.size(); i++) {
                assemblies[i] = potentialFunctions.open(filenames.get(i))
            }
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            assemblies = [activeAssembly]
        }

        Binding binding = getBinding()
        baseDir = assemblies[0].getFile()
        potentialFunctions = (PotentialsFunctions) binding.getVariable("functions")

        // Now that the molecular assembly is obtained, we will center each system to the origin and
        //  record distances from each molecule to the origin for each crystal

        MolecularAssembly[] expandedAssemblies = new MolecularAssembly[assemblies.size()]
        int[] minIndicies = new int[assemblies.size()]
        //Determine list of which atoms to compare/calculate RMSD.
        int[] comparisonAtoms
        if (centerAtoms != null) {
            comparisonAtoms = centerAtoms.clone()
            for (int i = 0; i < comparisonAtoms.size(); i++) {
                comparisonAtoms[i] = comparisonAtoms[i] - 1
                if (comparisonAtoms[i] < 0 ||
                        comparisonAtoms[i] > assemblies[0].getMolecules().get(0).getAtomList(false).size()) {
                    logger.severe(" Selected atoms are outside of molecular size.")
                }
            }
        } else {
            comparisonAtoms = new int[assemblies[0].getMolecules().get(0).getAtomList().size()]
            for (int i = 0; i < comparisonAtoms.size(); i++) {
                comparisonAtoms[i] = i
            }
        }

        if (original) { //TODO: Ensure this default is accurate
            numMolAlign = 3
        }
        int compareAtomsSize = comparisonAtoms.size()

        logger.info(format(" # Compare Atoms: %d", compareAtomsSize))
        for (int integerVal : comparisonAtoms) {
            logger.fine(format(" Compare Atom #: %d", integerVal + 1))
        }

        MolecularAssembly baseAssembly = null
        MolecularAssembly currentAssembly = null
        int moleculeSize = assemblies[0].getMolecules().get(0).getAtomList().size()
        int numAssemblies = assemblies.size()
        double[] baseCoords = new double[compareAtomsSize * numMolAlign * 3]
        double[] centerOfFirst = new double[moleculeSize * 3]
        int massIndex = 0
        int coordIndex = 0
        double[] assemblyMasses = null
        double[] assemblyCoords = null
        for (int i = 0; i < numAssemblies; i++) {
            String fileName = FilenameUtils.getBaseName(assemblies[i].getFile().getName())
            File saveLocationPDB = new File(fileName + ".pdb")
            int minIndex = -1
            if (i == 0) {
                expandedAssemblies[i] = PACCOMFunctions.GenerateBaseShell(assemblies[i], n_mol, original)
                baseAssembly = expandedAssemblies[i]
                currentAssembly = baseAssembly
            } else {
                //TODO: Determine appropriate largeShell value (currently flag)
                expandedAssemblies[i] = PACCOMFunctions.GenerateBaseShell(assemblies[i], largeShell, original)
                currentAssembly = expandedAssemblies[i]
            }
            minIndicies[i] = PACCOMFunctions.CenterMoleculeIndex(currentAssembly)
            logger.info("\n " + currentAssembly.getFile().getName() + " --> "
                    + baseAssembly.getFile().getName())
            logger.info(format(" AssemblySize: %d", currentAssembly.getMolecules().size()))

            /**         SUPERPOSE CENTER MOLECULES             **/
            logger.info("\n Superposing center-most molecule from each shell.\n")

            MSNode centerMolecule = currentAssembly.getMolecules().get(minIndicies[i])
            // Record masses of atoms in center molecule
            double[] centerMasses = new double[compareAtomsSize]
            // Record atom positions of center molecule
            double[] centerCoords = new double[compareAtomsSize * 3]

            ArrayList<Atom> centerMolAtomList = centerMolecule.getAtomList()
            ArrayList<Atom> assemblyAtomList = currentAssembly.getAtomList()
            int assemblyAtomsSize = assemblyAtomList.size()
            assemblyMasses = new double[assemblyAtomsSize]
            // Record coordinates for all atoms within the assembly in a manner the Superpose class will handle.
            assemblyCoords = new double[assemblyAtomsSize * 3]
            massIndex = 0
            coordIndex = 0
            for (Atom atom : assemblyAtomList) {
                assemblyMasses[massIndex++] = atom.getMass()
                assemblyCoords[coordIndex++] = atom.getX()
                assemblyCoords[coordIndex++] = atom.getY()
                assemblyCoords[coordIndex++] = atom.getZ()
            }
            massIndex = 0
            coordIndex = 0
            if (original) { //For conformity with original PACCOM.
                centerCoords = new double[3]
                centerMasses = new double[1]

                Atom atom = centerMolAtomList.get(comparisonAtoms[0])
                centerMasses[0] = 1.00
                centerCoords[0] = atom.getX()
                centerCoords[1] = atom.getY()
                centerCoords[2] = atom.getZ()

                logger.info(format(" Atom coordinates pre translation:"))
                logger.info(format(" N1: " + centerMolAtomList.get(comparisonAtoms[0]).toString()))
                logger.info(format(" N2: " + centerMolAtomList.get(comparisonAtoms[1]).toString()))
                logger.info(format(" N3: " + centerMolAtomList.get(comparisonAtoms[2]).toString()))
            } else {
                for (int index : comparisonAtoms) { //ca update
                    Atom atom = centerMolAtomList.get(index)
                    centerMasses[massIndex++] = atom.getMass()
                    centerCoords[coordIndex++] = atom.getX()
                    centerCoords[coordIndex++] = atom.getY()
                    centerCoords[coordIndex++] = atom.getZ()
                }
            }

            // Determine/apply translation necessary to get CoM for middle molecule to origin
            double[] translation = Superpose.calculateTranslation(centerCoords, centerMasses)
            String translationString = new String()
            for (double valueT : translation) {
                translationString += format("%16.8f", valueT) + "\t"
            }
            logger.info(format(" Translation Vector: (Move to Origin)\n %s", translationString))
            Superpose.applyTranslation(assemblyCoords, translation)

            // Update expanded assemblies to include translation
            coordIndex = 0
            for (Atom atom : assemblyAtomList) {
                double[] newXYZ = new double[]{assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++]}
                atom.setXYZ(newXYZ)
            }

            coordIndex = 0
            if (original) { //PACCOM conformity.
                Atom atom = centerMolAtomList.get(comparisonAtoms[0])
                centerCoords[coordIndex++] = atom.getX()
                centerCoords[coordIndex++] = atom.getY()
                centerCoords[coordIndex] = atom.getZ()
            } else {
                // Update centerCoords with the translation.
                for (int index : comparisonAtoms) { //ca update
                    Atom atom = centerMolAtomList.get(index)
                    centerCoords[coordIndex++] = atom.getX()
                    centerCoords[coordIndex++] = atom.getY()
                    centerCoords[coordIndex++] = atom.getZ()
                }
                double[] coMcenter = centerMolecule.getCenter(true)
                logger.info(format(" Center Mol CoM post Trans: %16.8f %16.8f %16.8f", coMcenter[0], coMcenter[1], coMcenter[2]))
            }

            // Record atom coordinates for center molecule of first assembly (Other assemblies will rotate to match first).
            if (i == 0) {
                centerOfFirst = Arrays.copyOf(centerCoords, centerCoords.length)
            }
            if (original) {
                // PACCOM conformity. Calculation of sigma, phai, and cita angles needed to get n2 on axis and n3 on plane
                centerOfFirst = new double[numMolAlign * 3]
                centerCoords = new double[numMolAlign * 3]
                Atom atomN1 = centerMolAtomList.get(comparisonAtoms[0])
                Atom atomN2 = centerMolAtomList.get(comparisonAtoms[1])
                Atom atomN3 = centerMolAtomList.get(comparisonAtoms[2])
                centerMasses = new double[numMolAlign]
                Arrays.fill(centerMasses, 1.0)

                double p1n2 = atomN2.getX()
                double q1n2 = atomN2.getY()
                double r1n2 = atomN2.getZ()
                double p2n2 = atomN2.getX()
                double q2n2 = atomN2.getY()

                double p1n3 = atomN3.getX()
                double q1n3 = atomN3.getY()
                double r1n3 = atomN3.getZ()
                double p2n3 = atomN3.getX()
                double q2n3 = atomN3.getY()

                logger.fine(format(" START: N1: %f %f %f", atomN1.getX(), atomN1.getY(), atomN1.getZ()))
                logger.fine(format(" START: N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format(" START: N3: %f %f %f", p1n3, q1n3, r1n3))

                double cita0 = acos(p1n2 / sqrt(pow(p1n2, 2) + pow(q1n2, 2)))
                double phai0 = acos(sqrt(pow(p1n2, 2) + pow(q1n2, 2)) / sqrt(pow(p1n2, 2) + pow(q1n2, 2) + pow(r1n2, 2)))
                if (q1n2 < 0.0) {
                    cita0 = -cita0
                }

                logger.fine(format(" cita: %f", cita0))
                logger.fine(format(" cos(cita): %f", cos(cita0)))
                logger.fine(format(" sin(cita): %f", sin(cita0)))

                p2n2 = p1n2 * cos(cita0) + q1n2 * sin(cita0)
                q1n2 = -p1n2 * sin(cita0) + q1n2 * cos(cita0)
                p2n3 = p1n3 * cos(cita0) + q1n3 * sin(cita0)
                q1n3 = -p1n3 * sin(cita0) + q1n3 * cos(cita0)

                p1n2 = p2n2
                p1n3 = p2n3

                logger.fine(format(" Step 1 N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format(" Step 1 N3: %f %f %f", p1n3, q1n3, r1n3))

                if (r1n2 > 0.0) {
                    phai0 = -phai0
                }

                logger.fine(format(" phai: %f", phai0))

                p2n2 = p1n2 * cos(phai0) - r1n2 * sin(phai0)
                r1n2 = p1n2 * sin(phai0) + r1n2 * cos(phai0)
                p2n3 = p1n3 * cos(phai0) - r1n3 * sin(phai0)
                r1n3 = p1n3 * sin(phai0) + r1n3 * cos(phai0)

                p1n2 = p2n2
                p1n3 = p2n3

                logger.fine(format(" Step 2 N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format(" Step 2 N3: %f %f %f", p1n3, q1n3, r1n3))

                double sigma0 = acos(q1n3 / sqrt(pow(q1n3, 2) + pow(r1n3, 2)))
                if (r1n3 < 0.0) {
                    sigma0 = -sigma0
                }

                logger.fine(format(" sigma: %f", sigma0))

                q2n2 = q1n2 * cos(sigma0) + r1n2 * sin(sigma0)
                r1n2 = -q1n2 * sin(sigma0) + r1n2 * cos(sigma0)
                q2n3 = q1n3 * cos(sigma0) + r1n3 * sin(sigma0)
                r1n3 = -q1n3 * sin(sigma0) + r1n3 * cos(sigma0)

                q1n2 = q2n2
                q1n3 = q2n3

                logger.fine(format(" DONE N1: %f %f %f", atomN1.getX(), atomN1.getY(), atomN1.getZ()))
                logger.fine(format(" DONE N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format(" DONE N3: %f %f %f", p1n3, q1n3, r1n3))

                centerOfFirst[0] = atomN1.getX()
                centerOfFirst[1] = atomN1.getY()
                centerOfFirst[2] = atomN1.getZ()
                centerOfFirst[3] = p1n2
                centerOfFirst[4] = q1n2
                centerOfFirst[5] = r1n2
                centerOfFirst[6] = p1n3
                centerOfFirst[7] = q1n3
                centerOfFirst[8] = r1n3

                centerCoords[0] = atomN1.getX()
                centerCoords[1] = atomN1.getY()
                centerCoords[2] = atomN1.getZ()
                centerCoords[3] = atomN2.getX()
                centerCoords[4] = atomN2.getY()
                centerCoords[5] = atomN2.getZ()
                centerCoords[6] = atomN3.getX()
                centerCoords[7] = atomN3.getY()
                centerCoords[8] = atomN3.getZ()
            }

            double[][] rotation = Superpose.calculateRotation(centerOfFirst, centerCoords, centerMasses)
            String rotString = new String()
            for (double[] rotRow : rotation) {
                rotString += "\n  "
                for (double rotCol : rotRow) {
                    rotString += format("%16.8f", rotCol) + "\t"
                }
            }
            logger.info(" First Rotation Matrix:" + rotString)
            Superpose.applyRotation(assemblyCoords, rotation)

            // Update expanded assemblies to include rotation.
            coordIndex = 0
            for (Atom atom : assemblyAtomList) {
                double[] newXYZ = new double[]{assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++]}
                atom.setXYZ(newXYZ)
            }

            if (original) { //PACCOM conformity.  later
                Atom atomN1 = centerMolAtomList.get(comparisonAtoms[0])
                Atom atomN2 = centerMolAtomList.get(comparisonAtoms[1])
                Atom atomN3 = centerMolAtomList.get(comparisonAtoms[2])
                centerCoords[0] = atomN1.getX()
                centerCoords[1] = atomN1.getY()
                centerCoords[2] = atomN1.getZ()
                centerCoords[3] = atomN2.getX()
                centerCoords[4] = atomN2.getY()
                centerCoords[5] = atomN2.getZ()
                centerCoords[6] = atomN3.getX()
                centerCoords[7] = atomN3.getY()
                centerCoords[8] = atomN3.getZ()
            } else {
                // Update centerCoords with the translation.
                coordIndex = 0
                for (int index : comparisonAtoms) { //ca update
                    Atom atom = centerMolAtomList.get(index)
                    centerCoords[coordIndex++] = atom.getX()
                    centerCoords[coordIndex++] = atom.getY()
                    centerCoords[coordIndex++] = atom.getZ()
                }
            }
            double rmsdCenter = Superpose.rmsd(centerOfFirst, centerCoords, centerMasses)
            logger.info(format(" RMSD for Center Molecules After First Rotation: %16.8f\n", rmsdCenter))
            //At this point, all expanded assemblies should have their center molecules at the origin.

            /**     PERFORM numMolAlign MOLECULE SUPERPOSITION     **/
            logger.info(format(" Superposing %d center-most molecules from each shell.", numMolAlign))
            int[] similarMolecules = new int[numMolAlign]
            Arrays.fill(similarMolecules, -1) //Later we check for duplicate molecules and need 0 to be an option
            //Store molecules matching molecules in first
            double[] centerOfFirstMidMol = new double[3]
            HashMap<Integer, Double> baseClosestMols = new HashMap<Integer, Double>()
            if (original) { //PACCOM Conformity
                similarMolecules = new int[numMolAlign]
                // Determine next numMolAlign closest molecules based on n1 atom in first assembly
                baseAssembly.getMolecules().get(minIndicies[0]).getAtomList().get(comparisonAtoms[0]).getXYZ(centerOfFirstMidMol)
                for (int j = 0; j < baseAssembly.getMolecules().size(); j++) {
                    double[] nextCenter = new double[3]
                    baseAssembly.getMolecules().get(j).getAtomList().get(comparisonAtoms[0]).getXYZ(nextCenter)
                    double value = sqrt(pow(centerOfFirstMidMol[0] - nextCenter[0], 2)
                            + pow(centerOfFirstMidMol[1] - nextCenter[1], 2)
                            + pow(centerOfFirstMidMol[2] - nextCenter[2], 2))
                    if (baseClosestMols.size() < numMolAlign) {
                        baseClosestMols.put(j, value)
                    } else {
                        double maxValue = Collections.max(baseClosestMols.values())
                        if (value < maxValue) {
                            Iterator<Map.Entry<Integer, Double>> it = baseClosestMols.entrySet().iterator()
                            while (it.hasNext()) {
                                Map.Entry<Integer, Double> item = it.next()
                                if (item.getValue() == maxValue) {
                                    it.remove()
                                    break
                                }
                            }
                            baseClosestMols.put(j, value)
                        }
                    }
                }
                baseClosestMols = PACCOMFunctions.SortHashMapByValue(baseClosestMols)
                logger.info(format(" baseClosestMol size: %d", baseClosestMols.size()))
                //Determine numMolAlign molecules in other assemblies closest to those in first
                // TODO Right now first come first serve, maybe optimize for closest molecule.

                //Determine numMolAlign molecules in other assemblies closest to those in first
                int baseIndex = 0
                for (int j : baseClosestMols.keySet()) {
                    double[] baseCenter = baseAssembly.getMolecules().get(j).getCenter(true)
                    int index = -1
                    double minDist = Double.MAX_VALUE
                    for (int k = 0; k < currentAssembly.getMolecules().size(); k++) {
                        boolean duplicate = false
                        double[] nextCenter = currentAssembly.getMolecules().get(k).getCenter(true)
                        double value = sqrt(pow(baseCenter[0] - nextCenter[0], 2)
                                + pow(baseCenter[1] - nextCenter[1], 2)
                                + pow(baseCenter[2] - nextCenter[2], 2))
                        if (value < minDist) {
                            for (int prevKey : similarMolecules) {
                                if (prevKey == k) {
                                    duplicate = true;
                                    break;
                                }
                            }
                            logger.fine(format(" Key: " + k + " Value: " + value + " is duplicate? " + duplicate))
                            if (!duplicate) {
                                minDist = value
                                index = k
                            }
                        }
                    }
                    logger.fine(format(" Set: %d", index))
                    similarMolecules[baseIndex++] = index
                }
            } else {
                // Determine next numMolAlign closest molecules in first assembly
                centerOfFirstMidMol = baseAssembly.getMolecules().get(minIndicies[0]).getCenter(true)
                for (int j = 0; j < baseAssembly.getMolecules().size(); j++) {
                    double[] nextCenter = baseAssembly.getMolecules().get(j).getCenter(true)
                    double value = sqrt(pow(centerOfFirstMidMol[0] - nextCenter[0], 2)
                            + pow(centerOfFirstMidMol[1] - nextCenter[1], 2)
                            + pow(centerOfFirstMidMol[2] - nextCenter[2], 2))
                    if (baseClosestMols.size() < numMolAlign) {
                        baseClosestMols.put(j, value)
                    } else {
                        double maxValue = Collections.max(baseClosestMols.values())
                        if (value < maxValue) {
                            Iterator<Map.Entry<Integer, Double>> it = baseClosestMols.entrySet().iterator()
                            while (it.hasNext()) {
                                Map.Entry<Integer, Double> item = it.next()
                                if (item.getValue() == maxValue) {
                                    it.remove()
                                    break;
                                }
                            }
                            baseClosestMols.put(j, value)
                        }
                    }
                }
                baseClosestMols = PACCOMFunctions.SortHashMapByValue(baseClosestMols)
                for (Map.Entry<Integer, Double> item : baseClosestMols) {
                    logger.fine(format(" BaseKey: %d\tValue: %f", item.key, item.value))
                }

                //Determine numMolAlign molecules in other assemblies closest to those in first
                int baseIndex = 0
                for (int j : baseClosestMols.keySet()) {
                    double[] baseCenter = baseAssembly.getMolecules().get(j).getCenter(true)
                    int index = -1
                    double minDist = Double.MAX_VALUE
                    for (int k = 0; k < currentAssembly.getMolecules().size(); k++) {
                        boolean duplicate = false
                        double[] nextCenter = currentAssembly.getMolecules().get(k).getCenter(true)
                        double value = sqrt(pow(baseCenter[0] - nextCenter[0], 2)
                                + pow(baseCenter[1] - nextCenter[1], 2)
                                + pow(baseCenter[2] - nextCenter[2], 2))
                        if (value < minDist) {
                            for (int prevKey : similarMolecules) {
                                if (prevKey == k) {
                                    duplicate = true;
                                    break;
                                }
                            }
                            logger.fine(format(" Key: " + k + " Value: " + value + " is duplicate? " + duplicate))
                            if (!duplicate) {
                                minDist = value
                                index = k
                            }
                        }
                    }
                    logger.fine(format(" Set: %d", index))
                    similarMolecules[baseIndex++] = index
                }
            }
            // Closest numMolAlign molecules have been determined for each shell
            // Prepare coordinates for Superpose algorithm
            centerCoords = new double[compareAtomsSize * numMolAlign * 3]
            massIndex = 0
            coordIndex = 0
            double[] masses = new double[numMolAlign * compareAtomsSize]
            double[] coordsToMove = new double[numMolAlign * compareAtomsSize * 3]
            if (original) {
                centerCoords = new double[compareAtomsSize * 3]
                coordsToMove = new double[compareAtomsSize * 3]

                //Calculation of sigma, phai, and cita angles needed to get n1 atoms of each molecule to desired loci
                Atom atomN1 = assemblyAtomList.get(moleculeSize * similarMolecules[0] + comparisonAtoms[0])
                Atom atomN2 = assemblyAtomList.get(moleculeSize * similarMolecules[0] + comparisonAtoms[1])
                Atom atomN3 = assemblyAtomList.get(moleculeSize * similarMolecules[0] + comparisonAtoms[2])

                double p1n2 = atomN2.getX()
                double q1n2 = atomN2.getY()
                double r1n2 = atomN2.getZ()
                double p2n2 = atomN2.getX()
                double q2n2 = atomN2.getY()

                double p1n3 = atomN3.getX()
                double q1n3 = atomN3.getY()
                double r1n3 = atomN3.getZ()
                double p2n3 = atomN3.getX()
                double q2n3 = atomN3.getY()

                logger.fine(format("START: N1: %f %f %f", atomN1.getX(), atomN1.getY(), atomN1.getZ()))
                logger.fine(format("START: N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format("START: N3: %f %f %f", p1n3, q1n3, r1n3))

                double cita0 = acos(p1n2 / sqrt(pow(p1n2, 2) + pow(q1n2, 2)))
                double phai0 = acos(sqrt(pow(p1n2, 2) + pow(q1n2, 2)) / sqrt(pow(p1n2, 2) + pow(q1n2, 2) + pow(r1n2, 2)))
                if (q1n2 < 0.0) {
                    cita0 = -cita0
                }

                logger.fine(format("cita: %f", cita0))
                logger.fine(format("cos(cita): %f", cos(cita0)))
                logger.fine(format("sin(cita): %f", sin(cita0)))

                p2n2 = p1n2 * cos(cita0) + q1n2 * sin(cita0)
                q1n2 = -p1n2 * sin(cita0) + q1n2 * cos(cita0)
                p2n3 = p1n3 * cos(cita0) + q1n3 * sin(cita0)
                q1n3 = -p1n3 * sin(cita0) + q1n3 * cos(cita0)

                p1n2 = p2n2
                p1n3 = p2n3

                logger.fine(format("Step 1 N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format("Step 1 N3: %f %f %f", p1n3, q1n3, r1n3))

                if (r1n2 > 0.0) {
                    phai0 = -phai0
                }

                logger.fine(format("phai: %f", phai0))

                p2n2 = p1n2 * cos(phai0) - r1n2 * sin(phai0)
                r1n2 = p1n2 * sin(phai0) + r1n2 * cos(phai0)
                p2n3 = p1n3 * cos(phai0) - r1n3 * sin(phai0)
                r1n3 = p1n3 * sin(phai0) + r1n3 * cos(phai0)

                p1n2 = p2n2
                p1n3 = p2n3

                logger.fine(format("Step 2 N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format("Step 2 N3: %f %f %f", p1n3, q1n3, r1n3))

                double sigma0 = acos(q1n3 / sqrt(pow(q1n3, 2) + pow(r1n3, 2)))
                if (r1n3 < 0.0) {
                    sigma0 = -sigma0
                }

                logger.fine(format("sigma: %f", sigma0))

                q2n2 = q1n2 * cos(sigma0) + r1n2 * sin(sigma0)
                r1n2 = -q1n2 * sin(sigma0) + r1n2 * cos(sigma0)
                q2n3 = q1n3 * cos(sigma0) + r1n3 * sin(sigma0)
                r1n3 = -q1n3 * sin(sigma0) + r1n3 * cos(sigma0)

                q1n2 = q2n2
                q1n3 = q2n3

                logger.fine(format("DONE: N1: %f %f %f", atomN1.getX(), atomN1.getY(), atomN1.getZ()))
                logger.fine(format("DONE N2: %f %f %f", p1n2, q1n2, r1n2))
                logger.fine(format("DONE N3: %f %f %f", p1n3, q1n3, r1n3))

                centerCoords[0] = atomN1.getX()
                centerCoords[1] = atomN1.getY()
                centerCoords[2] = atomN1.getZ()
                centerCoords[3] = p1n2
                centerCoords[4] = q1n2
                centerCoords[5] = r1n2
                centerCoords[6] = p1n3
                centerCoords[7] = q1n3
                centerCoords[8] = r1n3

                coordsToMove[0] = atomN1.getX()
                coordsToMove[1] = atomN1.getY()
                coordsToMove[2] = atomN1.getZ()
                coordsToMove[3] = atomN2.getX()
                coordsToMove[4] = atomN2.getY()
                coordsToMove[5] = atomN2.getZ()
                coordsToMove[6] = atomN3.getX()
                coordsToMove[7] = atomN3.getY()
                coordsToMove[8] = atomN3.getZ()

                baseCoords = centerCoords

                coordIndex = 0
                masses = new double[numMolAlign]
                Arrays.fill(masses, 1.0) //Oringial PACCOM does not use masses
                Arrays.fill(assemblyMasses, 1.0) //Oringial PACCOM does not use masses

                for (int j = 0; j < numMolAlign; j++) {
                    int xyzIterate = j * 3
                    double[] centerAtomCoords = new double[compareAtomsSize * 3]
                    double[] moveAtomCoords = new double[compareAtomsSize * 3]
                    double[] tempMasses = new double[compareAtomsSize]
                    Arrays.fill(tempMasses, 1.0) //Original PACCOM does not use masses
                    centerAtomCoords[0] = centerCoords[xyzIterate]
                    moveAtomCoords[0] = coordsToMove[xyzIterate++]

                    centerAtomCoords[1] = centerCoords[xyzIterate]
                    moveAtomCoords[1] = coordsToMove[xyzIterate++]

                    centerAtomCoords[2] = centerCoords[xyzIterate]
                    moveAtomCoords[2] = coordsToMove[xyzIterate++]

                    logger.info(format("x: %16.8f y: %16.8f z: %16.8f\tvs\tx: %16.8f y: %16.8f z: %16.8f",
                            centerAtomCoords[0], centerAtomCoords[1], centerAtomCoords[2],
                            moveAtomCoords[0], moveAtomCoords[1], moveAtomCoords[2]))

                    rmsdCenter = Superpose.rmsd(centerAtomCoords, moveAtomCoords, tempMasses)
                    logger.info(format(" Atom %d Alignment RMSD After First Rotation: %16.8f", j, rmsdCenter))
                }
            } else {
                if (i == 0) { //First assembly: record as center position to move remaining assemblies to
                    for (int j: baseClosestMols.keySet()) {
                        for (int compareAtom : comparisonAtoms) {
                            Atom atom = assemblyAtomList.get(moleculeSize * j + compareAtom)
                            baseCoords[coordIndex++] = atom.getX()
                            baseCoords[coordIndex++] = atom.getY()
                            baseCoords[coordIndex++] = atom.getZ()
                        }
                    }
                }
                massIndex = 0
                coordIndex = 0
                for (int j = 0; j < numMolAlign; j++) {
                    for (int compareAtom : comparisonAtoms) {
                        Atom atom = assemblyAtomList.get(moleculeSize * similarMolecules[j] + compareAtom)
                        masses[massIndex++] = atom.getMass()
                        coordsToMove[coordIndex++] = atom.getX()
                        coordsToMove[coordIndex++] = atom.getY()
                        coordsToMove[coordIndex++] = atom.getZ()
                    }
                }
                // Record coordinates for all atoms within the assembly in a manner the Superpose class will handle.
                coordIndex = 0
                massIndex = 0
                for (Atom atom : assemblyAtomList) {
                    assemblyMasses[massIndex++] = atom.getMass()
                    assemblyCoords[coordIndex++] = atom.getX()
                    assemblyCoords[coordIndex++] = atom.getY()
                    assemblyCoords[coordIndex++] = atom.getZ()
                }
            }

            rotation = Superpose.calculateRotation(baseCoords, coordsToMove, masses)
            rotString = new String()
            for (double[] rotRow : rotation) {
                rotString += "\n  "
                for (double rotCol : rotRow) {
                    rotString += format("%16.8f", rotCol) + "\t"
                }
            }
            logger.info(" Second Rotation Matrix:" + rotString)
            Superpose.applyRotation(assemblyCoords, rotation)

            // Update expanded assemblies to include rotation.
            coordIndex = 0
            for (Atom atom : assemblyAtomList) {
                double[] newXYZ = new double[]{assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++],
                        assemblyCoords[coordIndex++]}
                atom.setXYZ(newXYZ)
            }

            if (original) {
                // TEMP Update coords for after rotation RMSD comparions:
                for (int j = 0; j < numMolAlign; j++) {
                    massIndex = 0
                    coordIndex = 0
                    for (int compareAtom : comparisonAtoms) {
                        Atom atom = assemblyAtomList.get(moleculeSize * similarMolecules[j] + compareAtom)
                        coordsToMove[coordIndex++] = atom.getX()
                        coordsToMove[coordIndex++] = atom.getY()
                        coordsToMove[coordIndex++] = atom.getZ()
                    }
                }
                // Update coordinates for all atoms within the assembly in a manner the Superpose class will handle.
                coordIndex = 0
                massIndex = 0
                for (Atom atom : assemblyAtomList) {
                    assemblyMasses[massIndex++] = atom.getMass()
                    assemblyCoords[coordIndex++] = atom.getX()
                    assemblyCoords[coordIndex++] = atom.getY()
                    assemblyCoords[coordIndex++] = atom.getZ()
                }

                for (int j = 0; j < numMolAlign; j++) {
                    for (int k = 0; k < compareAtomsSize; k++) {
                        int xyzIterate = k * 3
                        double[] centerAtomCoords = new double[3]
                        double[] moveAtomCoords = new double[3]
                        double[] tempMasses = new double[1]
                        tempMasses[0] = masses[k]
                        centerAtomCoords[0] = baseCoords[xyzIterate]
                        moveAtomCoords[0] = coordsToMove[xyzIterate++]

                        centerAtomCoords[1] = baseCoords[xyzIterate]
                        moveAtomCoords[1] = coordsToMove[xyzIterate++]

                        centerAtomCoords[2] = baseCoords[xyzIterate]
                        moveAtomCoords[2] = coordsToMove[xyzIterate]

                        rmsdCenter = Superpose.rmsd(centerAtomCoords, moveAtomCoords, tempMasses)
                        logger.fine(format(" Mol %d Atom %d Alignment RMSD After Second Rotation: %16.8f", j, k, rmsdCenter))
                    }
                }
            } else {
                massIndex = 0
                coordIndex = 0
                for (int j = 0; j < numMolAlign; j++) {
                    for (int compareAtom : comparisonAtoms) {
                        Atom atom = assemblyAtomList.get(moleculeSize * similarMolecules[j] + compareAtom)
                        masses[massIndex++] = atom.getMass()
                        coordsToMove[coordIndex++] = atom.getX()
                        coordsToMove[coordIndex++] = atom.getY()
                        coordsToMove[coordIndex++] = atom.getZ()
                    }
                }
                // Update coordinates for all atoms within the assembly in a manner the Superpose class will handle.
                for (int j = 0; j < numMolAlign; j++) {
                    for (int k = 0; k < compareAtomsSize; k++) {
                        int xyzIterate = k * 3
                        double[] centerAtomCoords = new double[3]
                        double[] moveAtomCoords = new double[3]
                        double[] tempMasses = new double[1]
                        tempMasses[0] = masses[k]
                        centerAtomCoords[0] = baseCoords[xyzIterate]
                        moveAtomCoords[0] = coordsToMove[xyzIterate++]

                        centerAtomCoords[1] = baseCoords[xyzIterate]
                        moveAtomCoords[1] = coordsToMove[xyzIterate++]

                        centerAtomCoords[2] = baseCoords[xyzIterate]
                        moveAtomCoords[2] = coordsToMove[xyzIterate]

                        rmsdCenter = Superpose.rmsd(centerAtomCoords, moveAtomCoords, tempMasses)
                        logger.fine(format(" Mol %d Atom %d Alignment RMSD After Second Rotation: %16.8f", j, k, rmsdCenter))
                    }
                }
            }

            if (saveIntermediates) {
                PACCOMFunctions.SaveAssemblyPDB(currentAssembly, saveLocationPDB)
            }
            expandedAssemblies[i] = PACCOMFunctions.CutSubShell(assemblies[i], baseAssembly, currentAssembly, original)
            currentAssembly = expandedAssemblies[i]
            if (saveFiles) {
                PACCOMFunctions.SaveAssemblyPDB(currentAssembly, saveLocationPDB)
            }
        }

        // Calculate RMSD
        if (allVsAll) {
            logger.info("\n Calculating RMSD between each crystal:\n")
        } else {
            logger.info("\n Calculating RMSDs compared to first crystal:\n")
        }
        double[][] rmsdValues = new double[numAssemblies][numAssemblies]
        String topLabels = new String(" ")
        topLabels += "\t\t"
        String row = new String()
        for (int i = 0; i < numAssemblies; i++) {
            double[] compareCoords = new double[expandedAssemblies[i].getMolecules().size() * compareAtomsSize * 3]
            row += format(" %s\t", filenames.get(i))
            coordIndex = 0
//            logger.fine(format(" RMSD calculation assembly %d", i))
            for (MSNode molecule : expandedAssemblies[i].getMolecules()) {
                for (int compareAtom : comparisonAtoms) {
                    Atom atom = molecule.getAtomList().get(compareAtom)
//                    logger.finest(format(atom.toString()))
                    compareCoords[coordIndex++] = atom.getX()
                    compareCoords[coordIndex++] = atom.getY()
                    compareCoords[coordIndex++] = atom.getZ()
                }
            }
            for (int j = 0; j < numAssemblies; j++) {
                // Record coordinates for all atoms within the assembly in a manner the Superpose class will handle.
                int numMolRMSD = expandedAssemblies[j].getMolecules().size()
                logger.fine(format(" Number of Molecules to Align: %d", numMolAlign))
                logger.fine(format(" Number of Molecules for RMSD: %d", numMolRMSD))
                logger.fine(format(" Number of Atoms to Compare: %d", compareAtomsSize))
                logger.fine(format(" RMSD Denominator: %d\n", numMolRMSD * compareAtomsSize))
                assemblyMasses = new double[
                        expandedAssemblies[j].getMolecules().size() * compareAtomsSize]
                assemblyCoords = new double[
                        expandedAssemblies[j].getMolecules().size() * compareAtomsSize * 3]
                coordIndex = 0
                massIndex = 0
                for (MSNode molecule : expandedAssemblies[j].getMolecules()) {
                    for (int compareAtom : comparisonAtoms) {
                        Atom atom = molecule.getAtomList().get(compareAtom)
                        assemblyMasses[massIndex++] = atom.getMass()
                        assemblyCoords[coordIndex++] = atom.getX()
                        assemblyCoords[coordIndex++] = atom.getY()
                        assemblyCoords[coordIndex++] = atom.getZ()
                    }
                }
                // TODO: Determine normalization (numMolRMSD * compareAtomsSize) include masses?
                rmsdValues[i][j] = Superpose.rmsd(compareCoords, assemblyCoords)/(numMolRMSD * compareAtomsSize)
                if (i == 0) {
                    topLabels += format("%s  ", filenames.get(j))
                }
                row += format("%16.8f\t", rmsdValues[i][j])
                if (!allVsAll) {
                    j = numAssemblies
                }
            }
            row += format("\n")
            topLabels += format("\n")
        }
        logger.info(topLabels + row)
        return this
    }

    @Override
    List<Potential> getPotentials() {
        if (assemblies == null) {
            return new ArrayList<Potential>();
        } else {
            return Arrays.stream(assemblies).
                    filter { a -> a != null
                    }.
                    map { a -> a.getPotentialEnergy()
                    }.
                    filter { e -> e != null
                    }.
                    collect(Collectors.toList());
        }
    }
}