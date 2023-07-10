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

import com.google.common.collect.MinMaxPriorityQueue
import com.ibm.icu.util.CodePointTrie
import edu.rit.pj.Comm;
import ffx.algorithms.cli.AlgorithmsScript;
import ffx.numerics.math.Double3;
import ffx.numerics.math.HilbertCurveTransforms;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly
import ffx.potential.utils.EnergyException;
import ffx.potential.Utilities
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils
import org.apache.commons.lang.ArrayUtils
import org.apache.commons.math3.util.FastMath
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters

import java.util.logging.Level

import static ffx.potential.utils.Superpose.applyRotation;
import static ffx.potential.utils.Superpose.applyTranslation;
import static ffx.potential.utils.Superpose.calculateTranslation;
import static ffx.potential.utils.Superpose.translate;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.potential.utils.Superpose.rotate;
import static java.lang.String.format;
import static org.apache.commons.lang3.exception.ExceptionUtils.getStackTrace;
import static org.apache.commons.lang.ArrayUtils.toPrimitive;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cos
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * The TorsionScan command enumerates conformations of a molecule using torsional scans around rotatable bonds.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc TorsionScan &lt;filename&gt;
 */
@Command(description = " The TorsionScan command enumerates conformations of a molecule using torsional scans around rotatable bonds.",
        name = "TorsionScan")
class TorsionScan extends AlgorithmsScript {

    /**
     * --th or --theta Step size for bond rotations.
     */
    @Option(names = ['--th', '--theta'], paramLabel = '60.0', defaultValue = '60.0',
            description = "Step size for bond rotations (in Degrees).")
    private double increment = 60.0

    /**
     * --saveNumStates
     */
    @Option(names = ['--saveNumStates', "--sns"], paramLabel = '10', defaultValue = '10',
            description = 'Save this many of the lowest energy states per worker. This is the default.')
    private int saveNumStates = 10

    /**
     * --elimMax
     */
    @Option(names = ['--elimMax', "--em"], paramLabel = '2', defaultValue = '2',
            description = 'Eliminate bonds where one torsion causes high energies. Reduces the complexity of the search.')
    private int elimMax = 2

    /**
     * --saveAll
     */
    @Option(names = ['--saveAll'], paramLabel = 'false', defaultValue = 'false',
            description = 'Save out all states.')
    private boolean saveAll = false

    /**
     * --printEnergies
     */
    @Option(names = ['--noPrint'], paramLabel = 'true', defaultValue = 'true',
            description = 'Do not store or print energies and hilbert indicies at the end.')
    private boolean printOut = true

    /**
     * --sc or --staticComparison Hold angles fixed.
     */
    @Option(names = ['--sc', '--staticComparison'], paramLabel = "false", defaultValue = "false",
            description = 'If set, each bond is rotated independently (faster, but fewer permutations).')
    private static boolean staticCompare

    /**
     * --eliminationThreshold
     */
    @Option(names = ['--eliminationThreshold', '--et'], paramLabel = '50.0', defaultValue = '50.0',
            description = "Remove bonds that cause > this energy change during static analysis (kcal/mol).")
    private double eliminationThreshold = 50.0

    /**
     * The final argument should be the structure filename.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'Atomic coordinate file(s) to permute in XYZ format.')
    String filename = null

    /** Counter to determine how many rotations have been performed so far.*/
    public int progress = 0
    /** Number of bonds in the system. */
    public int numBonds = -1
    private double[] energies
    private long[] hilbertIndices
    /** Minimum permutaiton energy encountered. */
    private double minEnergy = Double.MAX_VALUE;

    /**
     * CrystalSuperpose Constructor.
     */
    TorsionScan() {
        this(new Binding())
    }

    /**
     * CrystalSuperpose Constructor.
     * @param binding Groovy Binding to use.
     */
    TorsionScan(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    @Override
    TorsionScan run() {
        // Direct is sufficient to scan for atom clashes.
        System.setProperty("polarization", "direct");

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // Ensure file exists.
        if (filename == null) {
            logger.info(helpString())
            return this
        }
        activeAssembly = getActiveAssembly(filename)
        SystemFilter systemFilter = algorithmFunctions.getFilter()
        MolecularAssembly ma = systemFilter.getActiveMolecularSystem();
        logger.info("")
        ForceFieldEnergy forceFieldEnergy = activeAssembly.potentialEnergy
        double[] x = new double[forceFieldEnergy.getNumberOfVariables()]

        // Get the bonds and rotation groups
        List<Bond> bonds = getTorsionalBonds(ma)
        logger.info(" " + bonds.size() + " torsional bonds found.")
        List<Atom[]> rotationGroups = getRotationGroups(bonds)
        int turns = (increment < 360.0) ? (int) (360.0/increment) : 1.0
        int bits = (int) FastMath.ceil(FastMath.log(2, turns))

        MinMaxPriorityQueue<StateContainer> queue
        if(saveAll){
            queue = MinMaxPriorityQueue.maximumSize(bonds.size()*turns as int).create()
        }
        else{
            queue = MinMaxPriorityQueue.maximumSize(saveNumStates).create()
        }

        // Potentially eliminate bonds
        eliminateBonds(elimMax, turns, bonds.size(), bonds, rotationGroups,
                queue,
                activeAssembly,
                forceFieldEnergy,
                x)

        // Calculate the maximum index of number of configurations using BigInteger
        BigInteger maxIndex = BigInteger.valueOf(2).pow(bits*bonds.size())
        BigInteger numConfigs = BigInteger.valueOf(turns).pow(bonds.size())
        logger.info("\n Performing " + turns + " turns about each bond.")
        logger.info(" Maximum number of configurations: " + numConfigs)
        logger.info(" Maximum hilbert index: " + maxIndex)
        long start = 0
        long end = maxIndex.longValue()
        long samples = numConfigs.longValue()

        // Parallel workers
        Comm world = Comm.world()
        long[] startIndices = new long[10]
        long jobsPerWorkerSplit = 0
        // Make worker assignment adjustments based on total number of configurations
        if(world.size() > 1 && !staticCompare){ // Workers not required for static comparison
            startIndices = new long[world.size()]
            // Assign each worker the same number of indices
            long jobsPerWorker = (long) (maxIndex.longValue() / world.size())
            logger.info("\n Jobs per worker: " + jobsPerWorker)
            jobsPerWorkerSplit = (long) (jobsPerWorker / world.size())
            logger.info(" Jobs per worker split: " + jobsPerWorkerSplit)
            // Assign starting indexes to each worker
            for(int i = 0; i < world.size(); i++){
                startIndices[i] = i * jobsPerWorker + world.rank() * jobsPerWorkerSplit
            }
            logger.info(" Worker " + world.rank() + " assigned indices: " + startIndices)
            if(printOut) {
                energies = new double[jobsPerWorker]
                hilbertIndices = new double[jobsPerWorker]
            }
            if(saveAll){
                MinMaxPriorityQueue<StateContainer> newQ = MinMaxPriorityQueue.maximumSize(jobsPerWorker as int).create()
                while(!queue.isEmpty()){
                    newQ.add(queue.removeFirst())
                }
                queue = newQ
            }
        }
        else if (!staticCompare && printOut){
            energies = new double[samples]
            hilbertIndices = new double[samples]
        }

        forceFieldEnergy.getCoordinates(x)
        double energy = forceFieldEnergy.energy(x, false)
        minEnergy = energy
        logger.info(" Initial energy: " + energy)

        if(!staticCompare) {
            if (world.size() > 1) {
                for (int i = 0; i < world.size(); i++) {
                    spinTorsions(bonds.size(), turns, bonds, rotationGroups,
                            startIndices[i], startIndices[i] + jobsPerWorkerSplit - 1,
                            queue,
                            activeAssembly,
                            forceFieldEnergy,
                            x)
                }
            } else {
                spinTorsions(bonds.size(), turns, bonds, rotationGroups, start, end,
                        queue,
                        activeAssembly,
                        forceFieldEnergy,
                        x)
            }
        }

        if(printOut && !staticCompare){
            logger.info("\n In order energies:")
            logger.info(" " + energies.toString())
            logger.info(" In order hilbert indices:")
            logger.info(" " + hilbertIndices.toString())
        }

        int count = 0
        logger.info("\n Saving " + queue.size() + " states.")
        String extension = "_rot.arc"
        if(world.size() > 1){
            extension = "_rank" + world.rank() + extension
        }
        File saveLocation = new File(FilenameUtils.removeExtension(filename) + extension)
        logger.info(" Logging structures into: " + saveLocation)
        logger.info(" Negative hilbert indicies values are from static analysis.")
        XYZFilter xyzFilter = new XYZFilter(saveLocation,
                activeAssembly,
                activeAssembly.getForceField(),
                activeAssembly.properties)
        while(!queue.empty){
            StateContainer toBeSaved = queue.removeFirst()
            AssemblyState assembly = toBeSaved.getState()
            assembly.revertState()
            forceFieldEnergy.getCoordinates(x)
            energy = toBeSaved.energy
            logger.info(format(" Writing to file. Configuration #%-5d energy: %-8.3f Hilbert index: %-10d",
                    (count+1),
                    energy,
                    toBeSaved.getIndex()))
            if(count == 0){
                xyzFilter.writeFile(saveLocation, false)
            }
            else{
                xyzFilter.writeFile(saveLocation, true)
            }
            count++
        }

        // Create properties/key file
        File key = new File(FilenameUtils.removeExtension(filename) + ".key");
        File properties = new File(FilenameUtils.removeExtension(filename) + ".properties");
        try {
            if (key.exists()) {
                File keyComparison = new File(FilenameUtils.removeExtension(filename) + "_rot.key");
                if (keyComparison.createNewFile()) {
                    FileUtils.copyFile(key, keyComparison);
                }
            } else if (properties.exists()) {
                File propertiesComparison = new File(FilenameUtils.removeExtension(filename) + "_rot.properties");
                if (propertiesComparison.createNewFile()) {
                    FileUtils.copyFile(properties, propertiesComparison);
                }
            } else{
                logger.info(" No key or properties file found.")
            }
        } catch (Exception e) {
            e.printStackTrace()
        }

        return this
    }

    static List<Bond> getTorsionalBonds(MolecularAssembly molecularAssembly) {
        List<Bond> torsionalBonds = new ArrayList<>();
        for (Bond bond : molecularAssembly.getBondList()) {
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            List<Bond> bond1 = a1.getBonds();
            int b1 = bond1.size();
            List<Bond> bond2 = a2.getBonds();
            int b2 = bond2.size();
            // Ignore single bonds
            if (b1 > 1 && b2 > 1){
                // Ignore carbon w 3 hydrogens
                if(a1.getAtomicNumber() == 6 && a1.getNumberOfBondedHydrogen() == 3 || a2.getAtomicNumber() == 6 && a2.getNumberOfBondedHydrogen() == 3){
                    continue
                }
                // Ignore ring atoms
                if (a1.isRing(a2)) {
                    continue
                }
                torsionalBonds.add(bond)
                logger.info(" Bond " + bond + " is a torsional bond.")
            }
        }
        return torsionalBonds;
    }

    List<Atom[]> getRotationGroups(List<Bond> bonds) {
        List<Atom[]> rotationGroups = new ArrayList<>();
        for(Bond bond : bonds){
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            // We should have two atoms with a spinnable bond
            List<Atom> a1List = new ArrayList<>();
            List<Atom> a2List = new ArrayList<>();
            // Search for rotation groups
            searchTorsions(a1, a1List, a2);
            searchTorsions(a2, a2List, a1);
            // Add smaller list to rotation groups
            if(a1List.size() > a2List.size()){
                rotationGroups.add(a2List.toArray() as Atom[])
            }
            else{
                rotationGroups.add(a1List.toArray() as Atom[])
            }
        }
        return rotationGroups
    }

    /**
     * Scanning through nTorsions^nTorsions (nTorsions x nTorsions x nTorsions x ...) array.
     * Uses a hilbert curve to get to each point in the discrete space with minimal rotations
     * between each state without recursion. Hilbert curve implementation based on OpenMM's
     * c++ implementation.
     * @param torsionalBonds
     * @param atomGroups
     * @param nTorsions
     * @param start
     * @param end
     * @param stateContainers
     * @param molecularAssembly
     * @param forceFieldEnergy
     */
    void spinTorsions(int nTorsionalBonds, int nTorsions, List<Bond> bonds, List<Atom[]> atomGroups,
                      long start, long end,
                      MinMaxPriorityQueue<StateContainer> stateContainers,
                      MolecularAssembly molecularAssembly,
                      ForceFieldEnergy forceFieldEnergy,
                      double[] x) {
        long numConfigs = (long) Math.pow(nTorsions, nTorsionalBonds)
        int nBits = (int) FastMath.ceil(FastMath.log(2, nTorsions))
        int[] currentState = new int[nTorsionalBonds]
        while(start <= end && progress < numConfigs){
            int[] newState = HilbertCurveTransforms.hilbertIndexToCoordinates(nTorsionalBonds, nBits, start)
            //logger.info(" Hilbert index: " + start + "; Coordinates: " + Arrays.toString(newState))
            if(ArrayUtils.contains(newState, nTorsions)){
                //logger.info(" Skipping state.")
                start++
                continue
            }
            // Permute from currentState to newState
            changeState(currentState, newState, nTorsions, bonds, atomGroups)
            // Update coordinates
            forceFieldEnergy.getCoordinates(x)
            // Calculate energy
            double energy = forceFieldEnergy.energy(x)
            // Log minimum energy
            if(energy < minEnergy){
                minEnergy = energy
                logger.info(format("\n New minimum energy: %12.5f", minEnergy))
                // Print out the coordinates
                StringBuilder sb = new StringBuilder(format(" Hilbert Index: %d; Coordinates: (", start))
                for (int i = 0; i < nTorsionalBonds; i++) {
                    sb.append(format("%d", newState[i]))
                    if (i < nTorsionalBonds - 1) {
                        sb.append(", ")
                    }
                }
                sb.append(")")
                logger.info(sb.toString())
            }
            // Add to queue
            stateContainers.add(new StateContainer(new AssemblyState(molecularAssembly), energy, start))
            currentState = newState
            if(printOut) {
                energies[progress] = energy
                hilbertIndices[progress] = start
            }
            start++
            progress++
        }
        if(progress == numConfigs){
            logger.info("\n Completed all configurations before end index.")
        }
    }

    static void changeState(int[] oldState, int[] newState, int nTorsions,
                               List<Bond> bonds, List<Atom[]> atoms) {
        for(int i = 0; i < oldState.length; i++) {
            if (oldState[i] != newState[i]) {
                // Add change required to go from oldState to newState to change array
                int change = newState[i] - oldState[i]
                // Apply change at this index to atoms
                int turnDegrees = (int) (change * (360 / nTorsions))
                // Get vector from atom to atom of bond i
                double[] u = new double[3]
                double[] translation = new double[3]
                double[] a1 = bonds.get(i).getAtom(0).getXYZ(new double[3])
                double[] a2 = bonds.get(i).getAtom(1).getXYZ(new double[3])
                if(atoms[i][0] == bonds.get(i).getAtom(0)){ // if a1 had the smaller rotation group
                    for(int j = 0; j < 3; j++){
                        // Rotation about this vector is same as the same rotation about the opposite vector?
                        u[j] = a1[j]-a2[j] // Vector defined as line segment from a2 to a1
                        translation[j] = -a1[j]
                    }
                    // Unit vector
                    double unit = 1/Math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
                    for (int j = 0; j < 3; j++) {
                        u[j] *= unit
                    }
                    for(int j = 0; j < atoms[i].length; j++){
                        atoms[i][j].move(translation)
                        rotateAbout(u, atoms[i][j], turnDegrees)
                        atoms[i][j].move(a1)
                    }
                }
                else{
                    for(int j = 0; j < 3; j++){
                        u[j] = a2[j]-a1[j] // Vector defined as line segment from a1 to a2
                        translation[j] = -a2[j]
                    }
                    double unit = 1/Math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
                    for (int j = 0; j < 3; j++) {
                        u[j] *= unit
                    }
                    for(int j = 0; j < atoms[i].length; j++){
                        atoms[i][j].move(translation)
                        rotateAbout(u, atoms[i][j], turnDegrees)
                        atoms[i][j].move(a2) // This line different then above
                    }
                }
            }
        }
    }

    /**
     * Rotate an atom about another.
     *
     * @param u a unit vector to rotate {@link ffx.potential.bonded.Atom} object about.
     * @param a2 a {@link ffx.potential.bonded.Atom} object to create axis of rotation.
     * @param theta Amount to rotate by in degrees.
     */
    static void rotateAbout(double[] u, Atom a2, double theta) {

        theta = toRadians(theta);

        double[][] rotation = new double[3][3];
        rotation[0][0] = cos(theta) + (u[0] * u[0]) * (1-cos(theta));
        rotation[0][1] = u[0]*u[1] * (1-cos(theta)) - u[2]*sin(theta);
        rotation[0][2] = u[0]*u[2] * (1-cos(theta)) + u[1]*sin(theta);
        rotation[1][0] = u[1]*u[0] * (1-cos(theta)) + u[2]*sin(theta);
        rotation[1][1] = cos(theta) + (u[1] * u[1]) * (1-cos(theta));
        rotation[1][2] = u[1]*u[2] * (1-cos(theta)) - u[0]*sin(theta);
        rotation[2][0] = u[2]*u[0] * (1-cos(theta)) - u[1]*sin(theta);
        rotation[2][1] = u[2]*u[1] * (1-cos(theta)) + u[0]*sin(theta);
        rotation[2][2] = cos(theta) + (u[2] * u[2]) * (1-cos(theta));

        double[] a2XYZ = new double[3];
        a2.getXYZ(a2XYZ);
        applyRotation(a2XYZ, rotation);
        a2.setXYZ(a2XYZ);
    }

    /**
     * Identify atoms that should be rotated.
     *
     * @param seed an {@link ffx.potential.bonded.Atom} object to rotate about.
     * @param atoms a list of {@link ffx.potential.bonded.Atom} objects to rotate.
     * @param notAtom Avoid this atom (wrong side of bond).
     */
    void searchTorsions(Atom seed, List<Atom> atoms, Atom notAtom){
        if (seed == null) {
            return;
        }
        atoms.add(seed);
        for (Bond b : seed.getBonds()) {
            Atom nextAtom = b.get1_2(seed);
            if (nextAtom == notAtom || atoms.contains(nextAtom)) { //nextAtom.getParent() != null ||
                continue;
            }
            searchTorsions(nextAtom, atoms, notAtom);
        }
    }

    void eliminateBonds(int numRemove, int nTorsions, int nBonds, List<Bond> bonds, List<Atom[]> atoms,
                         MinMaxPriorityQueue<StateContainer> states,
                         MolecularAssembly molecularAssembly,
                         ForceFieldEnergy forceFieldEnergy,
                         double[] x) {
        // Update coordinates
        forceFieldEnergy.getCoordinates(x)
        // Calculate energy
        double initialE = forceFieldEnergy.energy(x)
        ArrayList<Integer> remove = new ArrayList<>()
        int[] state = new int[nBonds]
        int[] oldState = new int[nBonds]
        for(int i = 0; i < nBonds; i++){
            for(int j = 1; j < nTorsions-1; j++){
                state[i] = j
                changeState(oldState, state, nTorsions, bonds, atoms)
                // Update coordinates
                forceFieldEnergy.getCoordinates(x)
                // Calculate energy
                double newEnergy = forceFieldEnergy.energy(x)
                if(newEnergy - initialE > eliminationThreshold){
                    remove.add(i)
                    break
                } else {
                    states.add(new StateContainer(new AssemblyState(molecularAssembly), newEnergy, -1))
                }
                state[i] = 0
            }
        }
        changeState(oldState, state, nTorsions, bonds, atoms)
        Arrays.sort(remove)
        logger.info("\n " + remove.size() + " bonds that cause large energy increase: " + remove)
        if(remove.size() > numRemove){
            remove = remove.subList(0, numRemove)
        }
        for(int i = remove.size()-1; i >= 0; i--){
            logger.info(" Removing bond: " + bonds.get(remove.get(i)))
            logger.info(" Bond index: " + remove.get(i))
            bonds.set(remove.get(i), null)
            atoms.set(remove.get(i), null)
        }
        bonds.removeAll(Collections.singleton(null))
        atoms.removeAll(Collections.singleton(null))
    }

    /**
     * Implements StateContainer to store the coordinates of a state and its energy
     */
    private class StateContainer implements Comparable<StateContainer> {

        private final AssemblyState state
        private final double e
        private final long hilbertIndex

        StateContainer(AssemblyState state, double e, long hilbertIndex) {
            this.state = state
            this.e = e
            this.hilbertIndex = hilbertIndex
        }

        AssemblyState getState() {
            return state
        }

        double getEnergy() {
            return e
        }

        long getIndex() {
            return hilbertIndex
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy())
        }
    }
}