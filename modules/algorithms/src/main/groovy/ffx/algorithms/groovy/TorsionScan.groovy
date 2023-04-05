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
package ffx.algorithms.groovy;

import ffx.algorithms.cli.AlgorithmsScript;
import ffx.numerics.math.Double3;
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
import static org.apache.commons.math3.util.FastMath.cos;
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
    private double increment

    /**
     * --ec or --energyCutoff Exclude conformations above this energy.
     */
    @Option(names = ['--ec', '--energyCutoff'], paramLabel = '30.0', defaultValue = '30.0',
            description = "Exclude conformations with relative energies more than this cutoff (in kcal/mol).")
    private double energyCutoff

    /**
     * --st or --saveTolerance Save conformations greater than this value.
     */
    @Option(names = ['--st', '--saveTolerance'], paramLabel = '-1.0', defaultValue = '-1.0',
            description = "Save out conformations with an RMSD larger than this value.")
    private double saveCutoff

    /**
     * --sc or --staticComparison Hold angles fixed.
     */
    @Option(names = ['--sc', '--staticComparison'], paramLabel = "false", defaultValue = "false",
            description = 'If set, each bond is rotated independently (faster, but fewer permutations).')
    private static boolean staticCompare

    /** Minimum permutaiton energy encountered. */
    private double minEnergy;

    /**
     * The final argument should be the structure filename.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'Atomic coordinate file(s) to permute in XYZ format.')
    List<String> filenames = null

    /** List to store conformations determined unique. */
    List<double[]> uniqueConformations = new ArrayList<>();

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
        if (filenames == null) {
            logger.info(helpString())
            return this
        }
        algorithmFunctions.openAll(filenames.get(0))
        SystemFilter systemFilter = algorithmFunctions.getFilter()

        // Define the filename to use for the RMSD values.
        String filename = filenames.get(0)
        MolecularAssembly ma = systemFilter.getActiveMolecularSystem();
        int size = ma.getAtomArray().size();
        double[] reference = new double[size * 3];
        int index = 0;
        Atom[] atoms = ma.getAtomArray();
        for(int i = 0; i < size; i++){
            Atom atom = atoms[i];
            reference[index++] = atom.getX();
            reference[index++] = atom.getY();
            reference[index++] = atom.getZ();
        }

        uniqueConformations.add(reference);

        List<Bond> bonds = ma.getBondList();
        int end = bonds.size();
        scanTorsions(ma, bonds, filename, 0, end);

        if(saveCutoff >= 0.0){
            for(double[] coords : uniqueConformations) {
                saveCoordinatesAsAssembly(ma, coords, filename)
            }
        }
        return this
    }

    /**
     * Recursive method to spin torsions.
     * @param ma Molecular assembly of interest
     * @param bonds Bonds in the molecule
     * @param filename Name of original file (altered before saving)
     * @param start Starting bond index
     * @param end Ending bond index
     */
    void scanTorsions(MolecularAssembly ma, List<Bond> bonds, String filename, int start, int end) {
        for(int i = start; i < end; i++) {
            Bond bond = bonds.get(i);
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            List<Bond> bond1 = a1.getBonds();
            int b1 = bond1.size();
            List<Bond> bond2 = a2.getBonds();
            int b2 = bond2.size();
            if (a1.isHydrogen() || a2.isHydrogen()) {
                continue;
            }
            if(a1.getAtomicNumber() == 6 && a1.getNumberOfBondedHydrogen() == 3 || a2.getAtomicNumber() == 6 && a2.getNumberOfBondedHydrogen() == 3){
                continue;
            }
            // No need to spin torsion if only bond...
            if (b1 > 1 && b2 > 1) {
                if (a1.isRing(a2)) {
                    // Don't spin rings... seems messy
                    if(logger.is(Level.FINER)){
                        logger.finer("Ring detected");
                    }

                    continue;
                }else{
                    if(logger.is(Level.FINER)){
                        logger.finer(" No rings detected.");
                    }

                }
                // We should have two atoms with a spinnable bond
                List<Atom> a1List = new ArrayList<>();
                List<Atom> a2List = new ArrayList<>();
                searchTorsions(a1, a1List, a2);
                searchTorsions(a2, a2List, a1);
                int listSize1 = a1List.size();
                int listSize2 = a2List.size();
                Atom[] aArray;
                Atom[] otherAtoms;
                if(listSize1 > listSize2){
                    aArray = a2List.toArray() as Atom[];
                    otherAtoms = a1List.toArray() as Atom[];
                }else{
                    aArray = a1List.toArray() as Atom[];
                    otherAtoms = a2List.toArray() as Atom[];
                }
                Double3 a1XYZ = a1.getXYZ()
                Double3 a2XYZ = a2.getXYZ()
                double[] x = new double[]{a1XYZ.get(0), a1XYZ.get(1), a1XYZ.get(2), a2XYZ.get(0), a2XYZ.get(1), a2XYZ.get(2)}
                double[] mass = new double[]{1.0,1.0}
                double[] translation = calculateTranslation(x, mass);
                applyTranslation(x, translation);
                // Unit vector to rotate about.
                double[] u = a2XYZ.sub(a1XYZ).normalize().get();
                double turns = (increment < 360.0) ? 360.0/increment : 1.0;
                for (int j = 0; j < turns; j++) {
                    for (Atom a : aArray) {
                        if (a == aArray[0]) {
                            //First atom is what we are rotating about...
                            continue;
                        }
                        a.move(translation);
                        rotateAbout(u, a, increment);
                        a.move(new double[]{-translation[0], -translation[1], -translation[2]});
                    }
                    // TODO Apply rotations above to molecule and add to MolecularAssembly to save out...
                    // All atoms should be rotated... Save structure.
                    Atom[] atoms = new Atom[listSize1 + listSize2];
                    int targetIndex = 1;
                    while(targetIndex <= (listSize1 + listSize2)){
                        boolean found = false;
                        for(Atom a: aArray){
                            if(a.getIndex() == targetIndex){
                                atoms[targetIndex++ - 1] = a;
                                found = true;
                                break;
                            }
                        }
                        if(!found){
                            for(Atom a: otherAtoms){
                                if(a.getIndex() == targetIndex){
                                    atoms[targetIndex++ - 1] = a;
                                    break;
                                }
                            }
                        }
                    }
                    if(saveCutoff >= 0.0){
                        addUnique(atoms, listSize1+listSize2);
                    }else{
                        // Save all torsions.
                        saveCoordinatesAsAssembly(ma, atoms, filename);
                    }

                    if (!staticCompare) {
                        scanTorsions(ma, bonds, filename, i+1, end);
                        if(i == end){
                            return;
                        }
                    }
                }
            }else{
                if(logger.is(Level.FINER)){
                    logger.finer(" One bond.");
                }
            }
        }
    }

    /**
     * Add new torsinal conformations to list of acquired structures.
     * @param atoms Atoms of the structure to be added.
     * @param size Number of atoms.
     * @return
     */
    boolean addUnique(final Atom[] atoms, int size){
        // Set up potentially new system.
        double[] test = new double[size * 3];
        double[] mass = new double[size];
        int index = 0;
        for(int i = 0; i < size; i++){
            Atom atom = atoms[i];
            mass[i] = 1.0
            test[index++] = atom.getX();
            test[index++] = atom.getY();
            test[index++] = atom.getZ();
        }
        double[] translation = calculateTranslation(test, mass);
        applyTranslation(test, translation);
        // Minimum difference between new structure and old list.
        double minDiff = Double.MAX_VALUE;
        for(Double[] reference0 : uniqueConformations){
            double[] reference = toPrimitive(reference0)
            translate(reference, mass);
            // Reference is rotated to match test
            rotate(test, reference, mass);

            double diff = rmsd(test, reference, mass);
            if(diff < minDiff){
                minDiff = diff;
            }
        }
        if(minDiff > saveCutoff){
            applyTranslation(test, new double[]{-translation[0], -translation[1], -translation[2]});
            uniqueConformations.add(test);
        }
        return true;
    }

    /**
     * Rotate an atom about another.
     *
     * @param a0 a {@link ffx.potential.bonded.Atom} object to rotate about.
     * @param a1 a {@link ffx.potential.bonded.Atom} object to create axis of rotation.
     * @param theta Amount to rotate by (degrees).
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
     * @param seed a {@link ffx.potential.bonded.Atom} object to rotate about.
     * @param atoms a list of {@link ffx.potential.bonded.Atom} objects to rotate.
     * @param searchDisulfide Whether to cross disulfides
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

    /**
     * Save File from atom array.
     * @param ma Starting MolecularAssembly object.
     * @param aArray Atoms to save out.
     * @param filename Location to save File.
     * @return
     */
    boolean saveCoordinatesAsAssembly(MolecularAssembly ma, Atom[] aArray, String filename) {
        // Add atoms from moved original to atom list.
        ArrayList<Atom> atomList = new ArrayList<>();
        for (Atom a : aArray) {
            Atom newAtom = new Atom(a.getIndex(), a.getName(), a.getAtomType(), a.getXYZ(null));
            atomList.add(newAtom);
        }
        saveAssemblyandFiles(ma, atomList, filename);
    }

    /**
     * Save File from XYZ coordinates.
     * @param ma Starting MolecularAssembly object.
     * @param newXYZ Desired coordiantes for atoms.
     * @param filename Location to save File.
     * @return
     */
    boolean saveCoordinatesAsAssembly(MolecularAssembly ma, double[] newXYZ, String filename) {
        // Add atoms from moved original to atom list.
        ArrayList<Atom> atomList = new ArrayList<>();
        Atom[] aArray = ma.getAtomArray();
        int size = aArray.length;
        for (int i = 0; i < size; i++) {
            Atom a = aArray[i];
            Atom newAtom = new Atom(a.getIndex(), a.getName(), a.getAtomType(), new double[]{newXYZ[i * 3], newXYZ[i * 3 + 1], newXYZ[i * 3 + 2]});
            atomList.add(newAtom);
        }
        saveAssemblyandFiles(ma, atomList, filename);
    }

    boolean saveAssemblyandFiles(MolecularAssembly ma, ArrayList<Atom> atomList, String filename){
        String fileName = FilenameUtils.removeExtension(filename);
        String description = "_rot";
        ArrayList<Bond> bondList = ma.getBondList();
        for (Bond bond : bondList) {
            Atom a1 = bond.getAtom(0);
            Atom a2 = bond.getAtom(1);
            //Indices stored as human-readable.
            int a1Ind = a1.getIndex() - 1;
            int a2Ind = a2.getIndex() - 1;
            Atom newA1 = atomList.get(a1Ind);
            Atom newA2 = atomList.get(a2Ind);
            Bond b = new Bond(newA1, newA2);
            b.setBondType(bond.getBondType());
        }

        File key = new File(fileName + ".key");
        File properties = new File(fileName + ".properties");
        if (key.exists()) {
            File keyComparison = new File(fileName + description + ".key");
            try {
                if (keyComparison.createNewFile()) {
                    FileUtils.copyFile(key, keyComparison);
                } else {
                    if(logger.is(Level.FINER)){
                        logger.finer(" Could not create properties file.");
                    }
                }
            } catch (Exception ex) {
                // Likely using properties file.
                if (logger.isLoggable(Level.FINER)) {
                    logger.finer(ex.toString());
                }
            }
        } else if (properties.exists()) {
            File propertiesComparison = new File(fileName + description + ".properties");
            try {
                if (propertiesComparison.createNewFile()) {
                    FileUtils.copyFile(properties, propertiesComparison);
                } else {
                    if (logger.isLoggable(Level.FINER)) {
                        logger.finer(" Could not create properties file.");
                    }
                }
            } catch (Exception ex) {
                // Likely not using a key/properties file... so PDB?
                logger.info(" Neither key nor properties file detected therefore only creating XYZ.");
                if (logger.isLoggable(Level.FINER)) {
                    logger.finer(ex.toString());
                    logger.finer(getStackTrace(ex));
                }
            }
        }

        MolecularAssembly saveAssembly = new MolecularAssembly(fileName);

        // Construct the force field for the expanded set of molecules
        saveAssembly.setForceField(ma.getForceField());
        saveAssembly.setCrystal(ma.getCrystal());

        // The biochemistry method is designed to load chemical entities into the
        // Polymer, Molecule, Water and Ion data structure.
        Utilities.biochemistry(saveAssembly, atomList);

        try {
            saveAssembly.setPotential(ForceFieldEnergy.energyFactory(saveAssembly));
            double energy = saveAssembly.potentialEnergy.energy();
            if(energy < minEnergy){
                minEnergy = energy;
            }
            double relativeEnergy = abs(energy - minEnergy)
            if( relativeEnergy > energyCutoff){
                logger.info(format(" Conformation energy (%9.4f) is greater than cutoff (%9.4f > %9.4f).", energy, relativeEnergy, energyCutoff));
                saveAssembly.destroy();
                return false;
            }else{
                logger.info(format(" Energy: %9.4f", energy));
                // superpose initial object for RMSD?
                File saveLocation = new File(fileName + description + ".arc");
                saveAssembly.setFile(saveLocation);
                XYZFilter xyzFilter = new XYZFilter(saveLocation, saveAssembly, saveAssembly.getForceField(), saveAssembly.getProperties());
                boolean success = xyzFilter.writeFile(saveLocation, true);

                saveAssembly.destroy();
                return success;
            }
        }catch(EnergyException eex){
            logger.info(format(" Unstable conformation skipped."));
            saveAssembly.destroy();
            return false
        }catch(Exception ex){
            logger.severe(getStackTrace(ex))
            saveAssembly.destroy();
            return false;
        }
    }
}
