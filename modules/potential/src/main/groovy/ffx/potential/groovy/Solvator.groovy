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

import ffx.crystal.Crystal
import ffx.crystal.ReplicatesCrystal
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Polymer
import ffx.utilities.Constants;

import org.apache.commons.io.FilenameUtils

import java.util.regex.Matcher
import java.util.regex.Pattern
import java.util.stream.Collectors
import java.util.stream.IntStream

import static java.lang.String.format

import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Solvator script puts a box of solvent around a solute.
 * <br>
 * Usage:
 * <br>
 * ffxc Solvator [options] &lt;filename&gt;
 */
@Command(description = " Solvate a system.", name = "ffxc Solvator")
class Solvator extends PotentialScript {

    /**
     * --sFi or --solventFile specifies the file to read for the solvent.
     */
    @Option(names = ['--sFi', '--solventFile'], paramLabel = "water",
            description = 'A file containing the solvent box to be used. A default will eventually be added.')
    private String solventFileName = null;

    /**
     * --iFi or --ionFile specifies a [optional] file to read for ion replacement. Must have a complementary
     * .ions file, with format ReferenceAtomStart ReferenceAtomEnd Concentration [Charge]. See forcefieldx/examples/nacl.ions.
     */
    @Option(names = ['--iFi', '--ionFile'], paramLabel = "[ions]",
            description = "Name of the file containing ions. Must also have a .ions file (e.g. nacl.pdb must also have nacl.ions). Default: no ions.")
    private String ionFileName = null;

    /**
     * -r or --rectangular uses a rectangular prism as the output rather than a cube;
     * this reduces overall box size, but is not recommended for simulations long enough to see solute rotation.
     */
    @Option(names = ['-r', '--rectangular'], paramLabel = "false",
            description = "Use a rectangular prism rather than a cube for solvation. Not recommended due to chance of solute rotations.")
    private boolean rectangular = false;

    /**
     * -p or --padding sets the minimum amount of solvent padding around the solute.
     */
    @Option(names = ['-p', '--padding'], paramLabel = "9.0",
            description = "Include at least this many Angstroms of solvent around the solute in each direction.")
    private double padding = 9.0;

    /**
     * -b or --boundary sets the minimum distance at which to keep solvent molecules around the solute.
     */
    @Option(names = ['-b', '--boundary'], paramLabel = "2.5",
            description = "Delete solvent molecules that infringe closer than this to the solute.")
    private double boundary = 2.5;

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * A domain decomposition scheme for solute atoms intended for rapidly assessing solute-solvent clashes.
     * Instead of a naive double loop over solute & solvent atoms, pre-place solute into this domain decomposition.
     * Then, for each new solvent atom, check which domain it belongs to, and check all atoms in that domain and
     * neighboring domains for possible clashes.
     */
    private class DomainDecomposition {
        // Solute atoms belonging in each domain.
        private final Atom[][][][] domains;
        // Solute atoms in that domain and all surrounding domains.
        private final Atom[][][][] withAdjacent;
        private final nX, nY, nZ;
        private final double domainLength;
        private final Crystal crystal;

        /**
         * Seturn a set of integers corresponding to domains that may neighbor a current domain along an axis.
         *
         * Checks for wraparound, out-of-bounds indices, and the final domain being subsized (and thus always
         * needing to include the penultimate domain).
         *
         * @param i Index to get neighbors for.
         * @param nI Number of domains along this axis.
         * @return All neighboring domain indices in this axis accounting for wraparound.
         */
        private static final Set<Integer> neighborsWithWraparound(int i, int nI) {
            List<Integer> boxesI = new ArrayList<>(3);
            boxesI.add(i);
            // Account for wraparound, and for box[nX-1] possibly being undersized (and thus needing to include box[nX-2])
            if (i == 0) {
                boxesI.add(nI - 1);
                boxesI.add(nI - 2);
                boxesI.add(1);
            } else if (i == (nI - 2)) {
                boxesI.add(nI - 3);
                boxesI.add(nI - 1);
                boxesI.add(0);
            } else if (i == (nI - 1)) {
                boxesI.add(nI - 2);
                boxesI.add(0);
            } else {
                boxesI.add(i-1);
                boxesI.add(i+1);
            }
            // Eliminate out-of-bounds values and duplicate values.
            return boxesI.stream().
                    filter({ it >= 0 && it < nI }).
                    collect(Collectors.toSet());
        }

        DomainDecomposition(double boundary, Atom[] soluteAtoms, Crystal crystal) {
            long time = -System.nanoTime();
            domainLength = 2.1 * boundary;
            nX = (int) (crystal.a / domainLength) + 1;
            nY = (int) (crystal.b / domainLength) + 1;
            nZ = (int) (crystal.c / domainLength) + 1;
            this.crystal = crystal;

            // First, determine how many Atoms per domain.
            // nAts will later be reused as a counter of "how many atoms placed in this domain".
            int[][][] nAts = new int[nX][nY][nZ];
            for (Atom at : soluteAtoms) {
                double[] xyz = new double[3];
                xyz = at.getXYZ(xyz);
                int i = (int) (xyz[0] / domainLength);
                int j = (int) (xyz[1] / domainLength);
                int k = (int) (xyz[2] / domainLength);
                nAts[i][j][k]++;
            }

            // Build the domains.
            domains = new Atom[nX][nY][nZ][];
            for (int i = 0; i < nX; i++) {
                for (int j = 0; j < nY; j++) {
                    for (int k = 0; k < nZ; k++) {
                        domains[i][j][k] = new Atom[nAts[i][j][k]];
                    }
                    // Reset the nAts array, which will be reused for placing atoms.
                    Arrays.fill(nAts[i][j], 0);
                }
            }

            // Add Atoms to domains.
            for (Atom at : soluteAtoms) {
                double[] xyz = new double[3];
                xyz = at.getXYZ(xyz);
                int i = (int) (xyz[0] / domainLength);
                int j = (int) (xyz[1] / domainLength);
                int k = (int) (xyz[2] / domainLength);
                domains[i][j][k][nAts[i][j][k]++] = at;
            }

            // Now construct a faster lookup array that contains all atoms from the domain and its neighbors.
            // I know it's a 6-deep loop, but it's still overall O(N) AFAICT; it does run < 300 ms for a 14000-atom system.
            withAdjacent = new Atom[nX][nY][nZ][];
            for (int i = 0; i < nX; i++) {
                Set<Integer> neighborsI = neighborsWithWraparound(i, nX);
                for (int j = 0; j < nY; j++) {
                    Set<Integer> neighborsJ = neighborsWithWraparound(j, nY);
                    for (int k = 0; k < nZ; k++) {
                        Set<Integer> neighborsK = neighborsWithWraparound(k, nZ);
                        List<Atom> neighborAtoms = new ArrayList<>();
                        for (int l : neighborsI) {
                            for (int m : neighborsJ) {
                                for (int n : neighborsK) {
                                    neighborAtoms.addAll(domains[l][m][n]);
                                }
                            }
                        }
                        logger.fine(format(" Constructing domain %3d-%3d-%3d with %d atoms.", i, j, k, neighborAtoms.size()));
                        logger.fine(format(" Neighbors along X: %s", neighborsI.toListString()));
                        logger.fine(format(" Neighbors along Y: %s", neighborsJ.toListString()));
                        logger.fine(format(" Neighbors along Z: %s\n", neighborsK.toListString()));
                        withAdjacent[i][j][k] = new Atom[neighborAtoms.size()];
                        withAdjacent[i][j][k] = neighborAtoms.toArray(withAdjacent[i][j][k]);
                    }
                }
            }

            int nBoxes = nX * nY * nZ;
            double avgAtoms = ((double) soluteAtoms.length) / nBoxes;
            time += System.nanoTime();
            logger.info(format(" Decomposed the solute into %d domains of side length %11.4g, averaging %10.3g atoms apiece in %8.3g sec.", nBoxes, domainLength, avgAtoms, (time * 1E-9)));
        }

        /**
         * Check if an Atom may clash with something registered with this DomainDecomposition
         * @param a A new solvent atom.
         * @param threshold Solute-solvent boundary.
         * @return If a clash is detected over all solute atoms and symops.
         */
        boolean checkClashes(Atom a, double threshold) {
            double[] xyz = new double[3];
            xyz = a.getXYZ(xyz);
            int i = (int) (xyz[0] / domainLength);
            int j = (int) (xyz[1] / domainLength);
            int k = (int) (xyz[2] / domainLength);

            return Arrays.stream(withAdjacent[i][j][k]).
                    anyMatch({ Atom s ->
                        double[] xyzS = new double[3];
                        xyzS = s.getXYZ(xyzS);
                        double dist = crystal.minDistOverSymOps(xyz, xyzS);
                        return dist < threshold;
                    });
        }
    }

    private class IonAddition {
        // Once again, modestly bad practice by exposing these fields.
        final double conc;
        final boolean toNeutralize;
        final Atom[] atoms;
        private final int nAts;
        final double[][] atomOffsets;
        final double charge;

        IonAddition(Atom[] atoms, double conc, boolean toNeutralize) {
            this.atoms = Arrays.copyOf(atoms, atoms.length);
            this.conc = conc;
            this.toNeutralize = toNeutralize;
            charge = Arrays.stream(atoms).mapToDouble({ it.getMultipoleType().getCharge() }).sum();
            nAts = atoms.length;

            if (nAts > 1) {
                double[] com = new double[3];
                atomOffsets = new double[nAts][3];
                Arrays.fill(com, 0);
                double sumMass = 0;

                // Calculate ionic center of mass.
                for (Atom atom : atoms) {
                    double mass = atom.getMass();
                    sumMass += mass;
                    double[] xyz = new double[3];
                    xyz = atom.getXYZ(xyz);

                    for (int i = 0; i < 3; i++) {
                        xyz[i] *= mass;
                        com[i] += xyz[i];
                    }
                }

                for (int i = 0; i < 3; i++) {
                    com[i] /= sumMass;
                }

                // Calculate positions as an offset from CoM.
                for (int i = 0; i < nAts; i++) {
                    Atom ai = atoms[i];
                    double[] xyz = new double[3];
                    xyz = ai.getXYZ(xyz);

                    for (int j = 0; j < 3; j++) {
                        atomOffsets[i][j] = xyz[j] - com[j];
                    }
                }
            } else {
                atomOffsets = new double[1][3];
                Arrays.fill(atomOffsets[0], 0.0);
            }
        }

        /**
         * Creates a new ion.
         *
         * @param com          Center of mass to place at.
         * @param currXYZIndex XYZ index of the atom directly preceding the ones to add.
         * @param chain        Chain to add ion to.
         * @param resSeq       Residue number to assign to the ion.
         * @return             Newly created ion atoms.
         */
        Atom[] createIon(double[] com, int currXYZIndex, char chain, int resSeq) {
            Atom[] newIonAts = new Atom[nAts];
            for (int i = 0; i < nAts; i++) {
                Atom fromAtom = atoms[i];
                double[] newXYZ = new double[3];
                System.arraycopy(com, 0, newXYZ, 0, 3);
                for (int j = 0 ; j < 3; j++) {
                    newXYZ[j] += atomOffsets[i][j];
                }
                Atom newAtom = new Atom(++currXYZIndex, fromAtom, newXYZ, resSeq, chain, Character.toString(chain));
                logger.fine(format(" New ion atom %s at chain %c on ion molecule %s-%d", newAtom, newAtom.getChainID(), newAtom.getResidueName(), newAtom.getResidueNumber()));
                newAtom.setHetero(true);
                newAtom.setResName(fromAtom.getResidueName());
                if (newAtom.getAltLoc() == null) {
                    newAtom.setAltLoc(' ' as char);
                }
                newIonAts[i] = newAtom;
            }
            return newIonAts;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(format(" Ion addition with %d atoms per ion, concentration %8.2g mM, charge %6.2f, used to neutralize %b", atoms.length, conc, charge, toNeutralize));
            for (Atom atom : atoms) {
                sb.append("\n").append(" Includes atom ").append(atom.toString());
            }
            return sb.toString();
        }
    }

    /**
     * Calculates the center of mass for an array of Atoms.
     *
     * @param atoms
     * @return Center of mass
     */
    private static double[] getCOM(Atom[] atoms) {
        double totMass = 0;
        double[] com = new double[3];
        double[] xyz = new double[3];
        for (Atom atom : atoms) {
            xyz = atom.getXYZ(xyz);
            double mass = atom.getMass();
            totMass += mass
            double invTotMass = 1.0 / totMass;
            for (int i = 0; i < 3; i++) {
                xyz[i] -= com[i];
                xyz[i] *= (mass * invTotMass);
                com[i] += xyz[i];
            }
        }
        return com;
    }

    /**
     * Pops a solvent molecule off the end of the list, and returns an ion at its center of mass.
     *
     * @param solvent      Source of solvent molecules (last member will be removed).
     * @param ia           Ion to replace it with.
     * @param currXYZIndex Index of the atom directly preceding the ion to be placed.
     * @param chain        Chain to place the ion into.
     * @param resSeq       Residue number to assign.
     * @return             Newly created set of ion atoms.
     */
    private static Atom[] swapIon(List<Atom[]> solvent, IonAddition ia, int currXYZIndex, char chain, int resSeq) {
        int nSolv = solvent.size() - 1;
        Atom[] lastSolvent = solvent.get(nSolv);
        solvent.remove(nSolv);
        double[] com = getCOM(lastSolvent);
        return ia.createIon(com, currXYZIndex, chain, resSeq);
    }

    /**
     * Execute the script.
     */
    @Override
    Solvator run() {
        if (!init()) {
            return this
        }

        if (!solventFileName) {
            logger.info(helpString());
            return this;
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        ForceFieldEnergy soluteEnergy = activeAssembly.getPotentialEnergy()
        Atom[] soluteAtoms = activeAssembly.getAtomArray()

        File solventFi = new File(solventFileName);
        MolecularAssembly solvent = potentialFunctions.open(solventFi);
        Atom[] baseSolventAtoms = solvent.getActiveAtomArray();

        Atom[] ionAtoms = null;
        File ionsFile = null;
        if (ionFileName) {
            MolecularAssembly ions = potentialFunctions.open(new File(ionFileName));
            ionAtoms = ions.getActiveAtomArray();
            String ionsName = FilenameUtils.removeExtension(ionFileName);
            ionsName = ionsName + ".ions";
            ionsFile = new File(ionsName);
            assert ionsFile != null && ionsFile.exists();
        }

        Crystal soluteCrystal = activeAssembly.getCrystal();
        Crystal solventCrystal = solvent.getCrystal();
        if (solventCrystal instanceof ReplicatesCrystal) {
            solventCrystal = ((ReplicatesCrystal) solventCrystal).getUnitCell();
        }

        if (!soluteCrystal.aperiodic()) {
            logger.severe(" Solute must be aperiodic to start!");
        }
        if (solventCrystal.aperiodic() || !soluteCrystal.spaceGroup.shortName.equalsIgnoreCase("P1")) {
            logger.severe(" Solvent box must be periodic (and P1)!");
        }

        int nSolute = soluteAtoms.length;
        double[][] soluteCoordinates = new double[nSolute][3];

        for (int i = 0; i < nSolute; i++) {
            Atom ati = soluteAtoms[i];
            double[] xyzi = new double[3];
            xyzi = ati.getXYZ(xyzi);
            System.arraycopy(xyzi, 0, soluteCoordinates[i], 0, 3);
        }

        double[] minSoluteXYZ = new double[3];
        double[] maxSoluteXYZ = new double[3];
        double[] soluteBoundingBox = new double[3];
        for (int i = 0; i < 3; i++) {
            minSoluteXYZ[i] = Arrays.stream(soluteCoordinates).mapToDouble({double[] xyz -> return xyz[i]}).min().getAsDouble();
            maxSoluteXYZ[i] = Arrays.stream(soluteCoordinates).mapToDouble({double[] xyz -> return xyz[i]}).max().getAsDouble();
            soluteBoundingBox[i] = maxSoluteXYZ[i] - minSoluteXYZ[i];
        }

        double soluteMass = Arrays.stream(soluteAtoms).mapToDouble({it.getMass()}).sum();
        double invSolMass = 1.0 / soluteMass;
        double[] origCoM = new double[3];
        for (int i = 0; i < nSolute; i++) {
            Atom ati = soluteAtoms[i];
            double massi = ati.getMass();
            massi *= invSolMass;
            double[] xyzi = new double[3];
            xyzi = ati.getXYZ(xyzi);
            for (int j = 0; j < 3; j++) {
                origCoM[j] += (xyzi[j] * massi);
            }
        }

        logger.fine(format(" Original CoM: %s", Arrays.toString(origCoM)));

        double[] newBox = new double[3];
        if (rectangular) {
            newBox = Arrays.stream(soluteBoundingBox).map({it + (2.0 * padding)}).toArray();
        } else {
            double soluteLinearSize = IntStream.range(0, nSolute - 1).parallel().mapToDouble({int i ->
                double[] xyzi = soluteCoordinates[i];
                return IntStream.range(i + 1, nSolute).mapToDouble({int j ->
                    double[] xyzj = soluteCoordinates[j];
                    double dist2 = 0;
                    for (int k = 0; k < 3; k++) {
                        double dx = xyzi[k] - xyzj[k];
                        dx *= dx;
                        dist2 += dx;
                    }
                    return dist2;
                }).max().getAsDouble();
            }).max().getAsDouble();
            soluteLinearSize = Math.sqrt(soluteLinearSize);
            soluteLinearSize += (2.0 * padding);
            Arrays.fill(newBox, soluteLinearSize);
        }

        logger.info(format(" Molecule will be solvated in a periodic box of size %10.4g, %10.4g, %10.4g", newBox[0], newBox[1], newBox[2]));

        double[] soluteTranslate = new double[3];
        for (int i = 0; i < 3; i++) {
            soluteTranslate[i] = 0.5 * newBox[i];
            soluteTranslate[i] -= origCoM[i];
        }
        for (Atom atom : soluteAtoms) {
            atom.move(soluteTranslate);
        }

        logger.fine(format(" Solute translated by %s", Arrays.toString(soluteTranslate)));

        solvent.moveAllIntoUnitCell();

        int nSolvent = baseSolventAtoms.length;
        double[][] baseSolventCoordinates = new double[nSolvent][3];
        for (int i = 0; i < nSolvent; i++) {
            Atom ati = baseSolventAtoms[i];
            double[] xyzi = new double[3];
            xyzi = ati.getXYZ(xyzi);
            System.arraycopy(xyzi, 0, baseSolventCoordinates[i], 0, 3);
        }

        logger.info(format(" Solute size: %12.4g, %12.4g, %12.4g", soluteBoundingBox[0], soluteBoundingBox[1], soluteBoundingBox[2]));

        int[] solventReplicas = new int[3];
        double[] solventBoxVectors = new double[3];
        solventBoxVectors[0] = solventCrystal.a;
        solventBoxVectors[1] = solventCrystal.b;
        solventBoxVectors[2] = solventCrystal.c;

        for (int i = 0; i < 3; i++) {
            solventReplicas[i] = (int) Math.ceil((newBox[i] / solventBoxVectors[i]));
        }

        Crystal newCrystal = new Crystal(newBox[0], newBox[1], newBox[2], 90, 90, 90, "P1");
        soluteEnergy.setCrystal(newCrystal);

        List<MSNode> bondedNodes = solvent.getAllBondedEntities();
        MSNode[] solventEntities = bondedNodes.toArray(new MSNode[bondedNodes.size()]);

        int nSolventMols = solventEntities.length;
        double[][] solventCOMs = new double[nSolventMols][];
        for (int i = 0; i < nSolventMols; i++) {
            List<Atom> moleculeAtoms = solventEntities[i].getAtomList();
            double totMass = moleculeAtoms.stream().mapToDouble({it.getAtomType().atomicWeight}).sum();
            double invMass = 1.0 / totMass;
            double[] com = new double[3];
            Arrays.fill(com, 0);

            for (Atom atom : moleculeAtoms) {
                double[] xyz = new double[3];
                xyz = atom.getXYZ(xyz);
                double mass = atom.getAtomType().atomicWeight;
                for (int j = 0; j < 3; j++) {
                    com[j] += (mass * xyz[j] * invMass);
                }
            }
            solventCrystal.toPrimaryCell(com, com);
            solventCOMs[i] = com;
        }

        double[] xyzOffset = new double[3];
        int currXYZIndex = Arrays.stream(soluteAtoms).
                mapToInt({it.getXyzIndex()})
                .max().getAsInt();
        int currResSeq = 1;

        char[] possibleChains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'] as char[];
        char solventChain = ' ';
        Set<Character> soluteChains = Arrays.stream(soluteAtoms).map({it.getChainID()}).collect(Collectors.toSet());
        for (char solvChainOpt : possibleChains) {
            if (!soluteChains.contains(solvChainOpt)) {
                solventChain = solvChainOpt;
                break;
            }
        }
        char ionChain = ' ';
        for (char solvChainOpt : possibleChains) {
            if (solvChainOpt != solventChain && !soluteChains.contains(solvChainOpt)) {
                ionChain = solvChainOpt;
                break;
            }
        }

        // Decompose the solute into domains so solute-solvent clashes can be quickly checked.
        DomainDecomposition ddc = new DomainDecomposition(boundary, soluteAtoms, newCrystal);

        // Unfortunately, Groovy treats ' ' as a String, not as a char.
        if (solventChain == ' '.charAt(0)) {
            logger.severe(" Could not find an unused character A-Z for the new solvent!");
        }
        logger.info(" New solvent molecules will be placed in chain ${solventChain}");

        // Accumulator for new solvent molecules.
        // Currently implemented as an ArrayList with the notion of removing from the end of the List.
        List<Atom[]> newMolecules = new ArrayList<>();

        for (int ai = 0; ai < solventReplicas[0]; ai++) {
            xyzOffset[0] = ai * solventBoxVectors[0];

            for (int bi = 0; bi < solventReplicas[1]; bi++) {
                xyzOffset[1] = bi * solventBoxVectors[1];

                for (int ci = 0; ci < solventReplicas[2]; ci++) {
                    xyzOffset[2] = ci * solventBoxVectors[2];

                    logger.info(format(" Tiling solvent replica %d,%d,%d", ai+1, bi+1, ci+1));

                    MoleculePlace: for (int i = 0; i < nSolventMols; i++) {
                        MSNode moli = solventEntities[i];
                        double[] comi = new double[3];
                        for (int j = 0; j < 3; j++) {
                            comi[j] = xyzOffset[j] + solventCOMs[i][j];
                            if (comi[j] < 0) {
                                logger.warning(format(" Skipping a copy of solvent molecule %d for violating minimum boundary 0,0,0. This should not occur!", i));
                                continue MoleculePlace;
                            } else if (comi[j] > newBox[j]) {
                                logger.fine(format(" Skipping a copy of solvent molecule %d for violating maximum boundary %12.4g,%12.4g,%12.4g", i, newBox[0], newBox[1], newBox[2]));
                                continue MoleculePlace;
                            }
                        }

                        logger.fine(format(" Placing a molecule at %f,%f,%f", comi[0], comi[1], comi[2]));

                        List<Atom> parentAtoms = moli.getAtomList();
                        int nMolAtoms = parentAtoms.size();
                        Atom[] newAtomArray = new Atom[nMolAtoms];

                        for (int atI = 0; atI < nMolAtoms; atI++) {
                            Atom parentAtom = parentAtoms.get(atI);
                            double[] newXYZ = new double[3];
                            newXYZ = parentAtom.getXYZ(newXYZ);
                            for (int j = 0; j < 3; j++) {
                                newXYZ[j] += xyzOffset[j];
                            }
                            Atom newAtom = new Atom(++currXYZIndex, parentAtom, newXYZ, currResSeq, solventChain, Character.toString(solventChain));
                            logger.fine(format(" New atom %s at chain %c on residue %s-%d", newAtom, newAtom.getChainID(), newAtom.getResidueName(), newAtom.getResidueNumber()));

                            newAtom.setHetero(!moli instanceof Polymer);
                            newAtom.setResName(moli.getName());
                            if (newAtom.getAltLoc() == null) {
                                newAtom.setAltLoc(' ' as char);
                            }
                            newAtomArray[atI] = newAtom;
                        }

                        boolean overlapFound = Arrays.stream(newAtomArray).
                                anyMatch({ Atom a -> ddc.checkClashes(a, boundary)});

                        if (overlapFound) {
                            logger.fine(format(" Skipping a copy of molecule %d for overlapping with the solute.", i));
                            continue MoleculePlace;
                        }

                        newMolecules.add(newAtomArray);
                        ++currResSeq;
                    }
                }
            }
        }

        Collections.shuffle(newMolecules);

        if (ionFileName) {
            logger.info(" Ions will be placed into chain ${ionChain}");
            double volume = newBox[0] * newBox[1] * newBox[2];
            // newCrystal.volume may also work.
            double ionsPermM = volume * Constants.LITERS_PER_CUBIC_ANGSTROM * Constants.AVOGADRO;

            List<IonAddition> byConc = new ArrayList<>();
            IonAddition neutAnion = null;
            IonAddition neutCation = null;

            BufferedReader br;
            Pattern ionicPattern = Pattern.compile("^\\s*([0-9]+) +([0-9]+) +([0-9]+(?:\\.[0-9]*)?|NEUT\\S*) *([-+]?[0-9]+)?");
            Pattern concPatt = Pattern.compile("^[0-9]+(?:\\.[0-9]*)?");
            // Parse .ions file to figure out which ions need to be added.
            try {
                br = new BufferedReader(new FileReader(ionsFile));
                String line = br.readLine();
                while (line != null) {
                    Matcher m = ionicPattern.matcher(line);
                    if (m.matches()) {
                        int start = Integer.parseInt(m.group(1));
                        int end = Integer.parseInt(m.group(2));
                        if (end < start) {
                            throw new IllegalArgumentException(" end < start");
                        }
                        int nAts = (end - start) + 1;

                        Atom[] atoms = new Atom[nAts];
                        for (int i = 0; i < nAts; i++) {
                            int atI = start + i - 1;
                            atoms[i] = ionAtoms[atI];
                        }

                        double conc = 0;
                        boolean toNeutralize;
                        if (concPatt.matcher(m.group(3)).matches()) {
                            conc = Double.parseDouble(m.group(3));
                            toNeutralize = false;
                        } else {
                            toNeutralize = true;
                        }
                        IonAddition ia = new IonAddition(atoms, conc, toNeutralize);
                        if (toNeutralize) {
                            if (ia.charge > 0) {
                                neutCation = ia;
                            } else if (ia.charge < 0) {
                                neutAnion = ia;
                            } else {
                                logger.severe(format(" Specified a neutralizing ion %s with no net charge!", ia.toString()));
                            }
                        } else {
                            byConc.add(ia);
                        }
                    }
                    line = br.readLine();
                }
            } finally {
                br?.close();
            }

            List<Atom[]> addedIons = new ArrayList<>();
            int ionResSeq = 0;

            // Begin swapping waters for ions.
            double initialCharge = activeAssembly.getCharge(false);
            logger.info(" Charge before addition of ions is ${initialCharge}");

            for (IonAddition ia : byConc) {
                logger.info(ia.toString());
                int nIons = Math.ceil(ionsPermM * ia.conc);
                if (nIons > newMolecules.size()) {
                    logger.severe(" Insufficient solvent molecules remain (${newMolecules.size()}) to add ${nIons} ions!");
                }
                logger.info(format(" Number of ions to place: %d\n", nIons));
                for (int i = 0; i < nIons; i++) {
                    Atom[] newIon = swapIon(newMolecules, ia, currXYZIndex, ionChain, ionResSeq++);
                    currXYZIndex += newIon.length;
                    addedIons.add(newIon);
                }
                initialCharge += nIons * ia.charge;
            }

            logger.info(" Charge before neutralization is ${initialCharge}")

            if (initialCharge > 0) {
                if (neutAnion == null) {
                    logger.info(" No counter-anion specified; system will be cationic at charge ${initialCharge}");
                } else {
                    logger.info(" Neutralizing system with ${neutAnion.toString()}");
                    double charge = neutAnion.charge;
                    int nAnions = Math.round(-1.0 * (initialCharge / charge));
                    double netCharge = initialCharge + (nAnions * charge);
                    if (nAnions > newMolecules.size()) {
                        logger.severe(" Insufficient solvent molecules remain (${newMolecules.size()}) to add ${nAnions} counter-anions!");
                    }
                    for (int i = 0; i < nAnions; i++) {
                        Atom[] newIon = swapIon(newMolecules, neutAnion, currXYZIndex, ionChain, ionResSeq++);
                        currXYZIndex += newIon.length;
                        addedIons.add(newIon);
                    }
                    logger.info(format(" System neutralized to %8.4g charge with %d counter-anions", netCharge, nAnions));
                }
            } else if (initialCharge < 0) {
                if (neutCation == null) {
                    logger.info(" No counter-cation specified; system will be anionic at charge ${initialCharge}");
                } else {
                    logger.info(" Neutralizing system with ${neutCation.toString()}");
                    double charge = neutCation.charge;
                    int nCations = Math.round(-1.0 * (initialCharge / charge));
                    double netCharge = initialCharge + (nCations * charge);
                    if (nCations > newMolecules.size()) {
                        logger.severe(" Insufficient solvent molecules remain (${newMolecules.size()}) to add ${nCations} counter-cations!");
                    }
                    for (int i = 0; i < nCations; i++) {
                        Atom[] newIon = swapIon(newMolecules, neutCation, currXYZIndex, ionChain, ionResSeq++);
                        currXYZIndex += newIon.length;
                        addedIons.add(newIon);
                    }
                    logger.info(format(" System neutralized to %8.4g charge with %d counter-cations", netCharge, nCations));
                }
            } else {
                logger.info(" System is neutral; no neutralizing ions needed.");
            }
            newMolecules.addAll(addedIons);
        }

        for (Atom[] atoms : newMolecules) {
            for (Atom atom : atoms) {
                activeAssembly.addMSNode(atom);
            }
        }

        String solvatedName = activeAssembly.getFile().getName().replaceFirst(~/\.[^.]+$/, ".pdb");
        potentialFunctions.saveAsPDB(activeAssembly, new File(solvatedName));

        return this
    }
}
