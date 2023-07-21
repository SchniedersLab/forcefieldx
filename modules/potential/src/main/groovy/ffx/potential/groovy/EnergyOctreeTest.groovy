//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import com.google.common.collect.MinMaxPriorityQueue
import com.github.quickhull3d.QuickHull3D
import edu.rit.pj.ParallelTeam
import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.numerics.atomic.AtomicDoubleArray
import ffx.potential.nonbonded.GeneralizedKirkwood
import ffx.potential.nonbonded.octree.*
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.nonbonded.NeighborList
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.utils.ConvexHullOps
import ffx.potential.utils.StructureMetrics
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.potential.utils.StructureMetrics.momentsOfInertia
import static ffx.potential.utils.StructureMetrics.radiusOfGyration
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc EnergyOctreeTest &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "EnergyOctreeTest")
class EnergyOctreeTest extends PotentialScript {

    @Mixin
    AtomSelectionOptions atomSelectionOptions

    /**
     * -m or --moments print out electrostatic moments.
     */
    @Option(names = ['-m', '--moments'], paramLabel = "false", defaultValue = "false",
            description = 'Print out electrostatic moments.')
    private boolean moments = false

    /**
     * -t or --theta set theta for Octree.
     */
    @Option(names = ['-t', '--theta'], paramLabel = "0.5", defaultValue = "0.5",
            description = 'Set theta for Octree.')
    private double theta = 0.5

    /**
     * --rg or --gyrate Print out the radius of gyration.
     */
    @Option(names = ['--rg', '--gyrate'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the radius of gyration.')
    private boolean gyrate = false

    /**
     * --in or --inertia Print out the moments of inertia.
     */
    @Option(names = ['--in', '--inertia'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the moments of inertia.')
    private boolean inertia = false

    /**
     * -g or --gradient to print out gradients.
     */
    @Option(names = ['-g', '--gradient'], paramLabel = "false", defaultValue = "false",
            description = 'Compute the atomic gradient as well as energy.')
    private boolean gradient = false

    /**
     * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
     */
    @Option(names = ['--fl', '--findLowest'], paramLabel = "0", defaultValue = "0",
            description = 'Return the n lowest energies from an ARC/PDB file.')
    private int fl = 0

    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @Option(names = ['-v', '--verbose'], paramLabel = "false", defaultValue = "false",
            description = "Print out all energy components for each snapshot.")
    private boolean verbose = false

    /**
     * --dc or --densityCutoff Collect structures above a specified density.
     */
    @Option(names = ['--dc', '--densityCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures above a specified density.")
    private double dCutoff = 0.0

    /**
     * --ec or --energyCutOff Collect structures below a specified energy range from the minimum energy.
     */
    @Option(names = ['--ec', '--energyCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
    private double eCutoff = 0.0

    /**
     * The final argument is a PDB or XYZ coordinate file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private String filename = null

    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null
//    public NeighborList.DomainDecomposition domainDecomposition = null
    private AssemblyState assemblyState = null

    // Parallel Team Set-up
    int cores = Runtime.getRuntime().availableProcessors()
    ParallelTeam parallelTeam = new ParallelTeam(cores)

    /**
     * Energy constructor.
     */
    EnergyOctreeTest() {
        this(new Binding())
    }

    /**
     * Energy constructor.
     * @param binding The Groovy Binding to use.
     */
    EnergyOctreeTest(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    EnergyOctreeTest run() {
        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // Load the MolecularAssembly.
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        // Set the filename.
        filename = activeAssembly.getFile().getAbsolutePath()

        logger.info("\n Running Energy on " + filename)

        // Apply atom selections
        atomSelectionOptions.setActiveAtoms(activeAssembly)

        // Selections for used atoms
        if (atomSelectionOptions.isUsedAtomSelectionSet()) {
//      logger.info("Used atom flag set")
            atomSelectionOptions.setUsed(activeAssembly)
        }

        // Get Radius of Gyration
//        double radgy = radiusOfGyration(activeAssembly.getAtomArray())
        // Define Crystal
//        Crystal octreeCrystal = new Crystal(3 * radgy, 3 * radgy, 3 * radgy, 90.0, 90.0, 90.0, "P1")
//        octreeCrystal.setAperiodic(true)
//        activeAssembly.setCrystal(octreeCrystal)

        // Set up Neighbor List
//        Atom[] nlAtoms = activeAssembly.getAtomArray()
//        NeighborList neighborList = new NeighborList(null, octreeCrystal,nlAtoms,2.0*radgy,2.0,parallelTeam)
//        int nAtoms = nlAtoms.length
//        int[][][] newLists = new int[1][nAtoms][]
//        double[][] coords = new double[1][nAtoms * 3]
//        for (int i = 0; i < nAtoms; i++) {
//            coords[0][i * 3] = nlAtoms[i].getX()
//            coords[0][i * 3 + 1] = nlAtoms[i].getY()
//            coords[0][i * 3 + 2] = nlAtoms[i].getZ()
//        }
//        boolean[] use = new boolean[nAtoms]
//        for (int i = 0; i < nAtoms; i++){
//            use[i] = true
//        }
//        neighborList.buildList(coords, newLists, use, true, true)

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        double[] collectGradient = new double[nVars]

        if (gradient) {
            double[] g = new double[nVars]
            int nAts = (int) (nVars / 3)
            energy = forceFieldEnergy.energyAndGradient(x, g, true)
//      logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
            for (int i = 0; i < nAts; i++) {
                int i3 = 3 * i
                collectGradient[i3] += g[i3]
                collectGradient[i3 + 1] += g[i3 + 1]
                collectGradient[i3 + 2] += g[i3 + 2]
//        logger.info(format("GRADIENT  %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
            }
        } else {
            energy = forceFieldEnergy.energy(x, true)
        }

        // Add Octree building
        // Use Convex Hull method to determine the maximum pairwise distance between any atoms
//        QuickHull3D quickHull3D = ConvexHullOps.constructHull(activeAssembly.getAtomArray())
//        double maxDist = ConvexHullOps.maxDist(quickHull3D)
//        NeighborList neighborList = forceFieldEnergy.getNeighborList()
//        int nAtoms = activeAssembly.getAtomArray().size()
//        logger.info(format("** maxDist %4.8f nAtoms %d",maxDist,nAtoms))
//        neighborList.buildOctree(nAtoms, forceFieldEnergy.getCrystal(), maxDist+2.0,8)

        // Testing other Octree Method
        // Create Octree particles from Atoms
        Atom[] atoms = activeAssembly.getActiveAtomArray()
        GeneralizedKirkwood generalizedKirkwood = forceFieldEnergy.getGK()
        AtomicDoubleArray selfEnergies = generalizedKirkwood.getSelfEnergy()
        double[] selfEnergiesArray = new double[selfEnergies.size()]

        for (int i = 0; i < selfEnergies.size(); i++) {
            selfEnergiesArray[i] = selfEnergies.get(i)
        }

        ArrayList<OctreeParticle> particles = new ArrayList<>()
        for (Atom atom : atoms){
            double[] coords = new double[]{atom.getX(),atom.getY(),atom.getZ()}
            OctreeParticle particle = new OctreeParticle(coords,1.0,atom.getCharge())
            particles.add(particle)
        }

        // Use Particles to build a tree
        int nCritical = atoms.size() * 0.1
        logger.info(format("nCritical = %d",nCritical))
        logger.info(format("theta = %4.4f",theta))
//        double theta = 0.0001
        Octree octree = new Octree(nCritical, atoms, theta, forceFieldEnergy.getPmeNode().globalMultipole,activeAssembly,gradient,selfEnergiesArray)

        // Determine maximum separation dis0tance between any two atoms
        QuickHull3D quickHull3D = ConvexHullOps.constructHull(atoms)
        double maxDist = ConvexHullOps.maxDist(quickHull3D)

        // Compute geometric center of atoms
//        double[] center = activeAssembly.computeGeometricCenter(atoms)
        double[] center = activeAssembly.computeCenterOfMass(atoms)

        // Define root OctreeCell
        OctreeCell root = new OctreeCell(nCritical)
        root.setX(center[0])
        root.setY(center[1])
        root.setZ(center[2])
        root.setR(maxDist*0.5)
        logger.info(format("Center of octree root: %4.3f %4.3f %4.3f",root.getX(),root.getY(),root.getZ()))
        logger.info(format("Sidelength of root cell = %4.3f",maxDist));

        // Build tree from root OctreeCell
        octree.buildTree(root)
//        logger.info(format("Printing Cells from Octree"))
//        octree.printCells()

        // Generate NF and FF lists
        octree.neighborCount(1)

        // Compute Multipoles for all leaf cells
        octree.getMultipole(0)

        // Perform the upward sweep to calculate parent cells' multipoles based on child cells' multipoles
        octree.upwardSweep()

        // Calculate the Potential at each atom
        octree.evalPotential()
        octree.printPhi()
/*
        if (moments) {
            logger.info("** Moments being calculated")
            NeighborList.Cell[][][] cells = neighborList.getCells()
            int[] numDivisions = neighborList.getNumDivisions()
            forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false)
            for (int i = 0; i < numDivisions[0]; i++) {
                for (int j = 0; j < numDivisions[1]; j++) {
                    for (int k = 0; k < numDivisions[2]; k++) {
                        double sideLength = maxDist+2.0
                        logger.info(format("** Side length passed to computeMomentsGeometric %2.4f",sideLength))
                        forceFieldEnergy.getPmeNode().computeMomentsGeometric(activeAssembly.getAtomArray(), false, cells[i][j][k], sideLength)
                    }
                }
            }
        }
*/
        if (gyrate) {
            double rg = radiusOfGyration(activeAssembly.getActiveAtomArray())
            logger.info(format(" Radius of gyration:           %10.5f A", rg))
        }

        if (inertia) {
            double[][] inertiaValue =
                    momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true)
        }

        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {

            int numSnaps = fl
            double lowestEnergy = energy
            assemblyState = new AssemblyState(activeAssembly)
            int index = 1

            // Making the MinMax priority queue that will expel the largest entry when it reaches its maximum size N/
            MinMaxPriorityQueue<StateContainer> lowestEnergyQueue = null
            if (fl > 0) {
                lowestEnergyQueue = MinMaxPriorityQueue.maximumSize(numSnaps).create()
                lowestEnergyQueue.add(new StateContainer(assemblyState, lowestEnergy))
            }

            int numModels = systemFilter.countNumModels()
            //Store densities in ordered encountered (used in density cutoff implementation).
            double[] densities = new double[numModels]
            //Store energies in ordered encountered (used in energy cutoff implementation).
            double[] energies = new double[numModels]
            energies[0] = energy

            while (systemFilter.readNext()) {
                index++
                Crystal crystal = activeAssembly.getCrystal()
                densities[index - 1] = crystal.getDensity(activeAssembly.getMass())
                forceFieldEnergy.setCrystal(crystal)
                forceFieldEnergy.getCoordinates(x)
                if (verbose) {
                    logger.info(format(" Snapshot %4d", index))
                    if (!crystal.aperiodic()) {
                        logger.info(format("\n Density:                                %6.3f (g/cc)",
                                crystal.getDensity(activeAssembly.getMass())))
                    }
                    if (gradient) {
                        double[] g = new double[nVars]
                        int nAts = (int) (nVars / 3)
                        energy = forceFieldEnergy.energyAndGradient(x, g, true)
//            logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
                        for (int i = 0; i < nAts; i++) {
                            int i3 = 3 * i
                            collectGradient[i3] += g[i3]
                            collectGradient[i3 + 1] += g[i3 + 1]
                            collectGradient[i3 + 2] += g[i3 + 2]
//              logger.info(format("GRADIENT  %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
                        }
                    } else {
                        energy = forceFieldEnergy.energy(x, true)
                    }
                } else {
                    energy = forceFieldEnergy.energy(x, false)
                    logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))
                }

                energies[index - 1] = energy

                //Update lowest encountered energy
                if (energy < lowestEnergy) {
                    lowestEnergy = energy
                }

                if (fl > 0) {
                    lowestEnergyQueue.add(new StateContainer(new AssemblyState(activeAssembly), energy))
                }

                if (moments) {
                    forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false)
                }

                if (gyrate) {
                    double rg = radiusOfGyration(activeAssembly.getActiveAtomArray())
                    logger.info(format(" Radius of gyration:          %10.5f A", rg))
                }

                if (inertia) {
                    double[][] inertiaValue =
                            momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true)
                }
            }

            if (gradient) {
                double inverseIndex = 1.0 / (double) index
                int nAts = (int) (nVars / 3)
                logger.info(format("    Atom       X, Y and Z Mean Gradient Components (kcal/mol/A)"))
                for (int i = 0; i < nAts; i++) {
                    int i3 = 3 * i
                    logger.info(format("GRADIENT  %7d %16.8f %16.8f %16.8f", i + 1, collectGradient[i3] * inverseIndex,
                            collectGradient[i3 + 1] * inverseIndex, collectGradient[i3 + 2] * inverseIndex))
                }
            }


            // Use the current base directory, or update if necessary based on the given filename.
            String dirString = getBaseDirString(filename)

            // If cutoffs have been selected create an ARC or PDB to store structures that satisfy cutoff.
            if ((eCutoff > 0.0 || dCutoff > 0.0) && numModels > 1) {
                systemFilter.readNext(true)
                activeAssembly = systemFilter.getActiveMolecularSystem()

                String name = getName(filename)
                String ext = getExtension(name)
                name = removeExtension(name)

                File saveFile
                SystemFilter writeFilter
                if (ext.toUpperCase().contains("XYZ")) {
                    saveFile = new File(dirString + name + ".xyz")
                    writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                            activeAssembly.getProperties())
                    potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
                } else if (ext.toUpperCase().contains("ARC")) {
                    saveFile = new File(dirString + name + ".arc")
                    saveFile = potentialFunctions.versionFile(saveFile)
                    writeFilter = new XYZFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                            activeAssembly.getProperties())
                    logger.info(" Saving to " + saveFile.getAbsolutePath())
                    saveFile.createNewFile()
                } else {
                    saveFile = new File(dirString + name + ".pdb")
                    saveFile = potentialFunctions.versionFile(saveFile)
                    writeFilter = new PDBFilter(saveFile, activeAssembly, activeAssembly.getForceField(),
                            activeAssembly.getProperties())
                    if (numModels > 1) {
                        writeFilter.setModelNumbering(0)
                    }
                    saveFile.createNewFile()
                }

                // Determine if each structure meets the cutoff condition
                if (systemFilter instanceof XYZFilter || systemFilter instanceof PDBFilter) {
                    int structNum = 0
                    if (eCutoff > 0.0) {
                        logger.info(
                                format("Lowest Energy of: %16.4f\n Saving structures with energy below: %16.4f",
                                        lowestEnergy,
                                        lowestEnergy + eCutoff))
                    }

                    do {
                        if (dCutoff > 0.0 && eCutoff > 0.0) {
                            if (energies[structNum] < lowestEnergy + eCutoff && densities[structNum] > dCutoff) {
                                if (systemFilter instanceof PDBFilter) {
                                    PDBFilter pdbFilter = (PDBFilter) systemFilter
                                    pdbFilter.writeFile(saveFile, true, false, false)
                                    saveFile.append("ENDMDL\n")
                                } else if (systemFilter instanceof XYZFilter) {
                                    writeFilter.writeFile(saveFile, true)
                                }
                            }
                        } else if (dCutoff > 0.0) {
                            if (densities[structNum] > dCutoff) {
                                if (systemFilter instanceof PDBFilter) {
                                    PDBFilter pdbFilter = (PDBFilter) systemFilter
                                    pdbFilter.writeFile(saveFile, true, false, false)
                                    saveFile.append("ENDMDL\n")
                                } else if (systemFilter instanceof XYZFilter) {
                                    writeFilter.writeFile(saveFile, true)
                                }
                            }
                        } else if (eCutoff > 0.0) {
                            if (energies[structNum] < lowestEnergy + eCutoff) {
                                if (systemFilter instanceof PDBFilter) {
                                    PDBFilter pdbFilter = (PDBFilter) systemFilter
                                    pdbFilter.writeFile(saveFile, true, false, false)
                                    saveFile.append("ENDMDL\n")
                                } else if (systemFilter instanceof XYZFilter) {
                                    writeFilter.writeFile(saveFile, true)
                                }
                            }
                        }
                        structNum++
                    } while (systemFilter.readNext())
                    if (systemFilter instanceof PDBFilter) {
                        saveFile.append("END\n")
                    }
                }
            }

            if (fl > 0) {
                if (numSnaps > index) {
                    logger.warning(format(
                            " Requested %d snapshots, but file %s has only %d snapshots. All %d energies will be reported",
                            numSnaps, filename, index, index))
                    numSnaps = index
                }

                String name = getName(filename)

                for (int i = 0; i < numSnaps - 1; i++) {
                    StateContainer savedState = lowestEnergyQueue.removeLast()
                    AssemblyState finalAssembly = savedState.getState()
                    finalAssembly.revertState()
                    double finalEnergy = savedState.getEnergy()
                    logger.info(format(" The potential energy found is %16.8f (kcal/mol)", finalEnergy))
                    File saveFile = potentialFunctions.versionFile(new File(dirString + name))
                    MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                    potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
                }

                StateContainer savedState = lowestEnergyQueue.removeLast()
                AssemblyState lowestAssembly = savedState.getState()
                lowestEnergy = savedState.getEnergy()

                assemblyState.revertState()
                logger.info(format(" The lowest potential energy found is %16.8f (kcal/mol)", lowestEnergy))

                // Prints our final energy (which will be the lowest energy
                File saveFile = potentialFunctions.versionFile(new File(dirString + name))
                MolecularAssembly molecularAssembly = assemblyState.getMolecularAssembly()
                potentialFunctions.saveAsPDB(molecularAssembly, saveFile)
            }
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (forceFieldEnergy == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList((Potential) forceFieldEnergy)
        }
        return potentials
    }

    private class StateContainer implements Comparable<StateContainer> {

        private final AssemblyState state
        private final double e

        StateContainer(AssemblyState state, double e) {
            this.state = state
            this.e = e
        }

        AssemblyState getState() {
            return state
        }

        double getEnergy() {
            return e
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy())
        }
    }

    /**
     * This entry point is being used to test GraalVM ahead-of-time compilation.
     * @param args Command line arguments.
     */
    public static void main(String... args) {
        Binding binding = new Binding(args)
        Energy energyScript = new Energy(binding)
        energyScript.run()
        System.exit(0)
    }
}

