// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
// ******************************************************************************
package ffx.potential.utils;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.pow;
import static java.lang.String.format;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXScript;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import java.io.File;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * PACCOMFunctions contains methods utilized in the PACCOM.groovy file.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 */
public class PACCOMFunctions {

    /**
     * The logger for this class.
     */
    protected static final Logger logger = Logger.getLogger(FFXScript.class.getName());

    /**
     * Distance calculator from Apache Commons.
     */
    protected static final EuclideanDistance distance = new EuclideanDistance();

    /**
     * Generate a molecular assembly of a sphere containing a specified number of molecules.
     *
     * @param molecularAssembly Base molecular assembly to expand.
     * @param nMolecules        Number of molecules to include in sphere.
     * @param original          If following the original version of PACCOM (by Okimasa Okada)
     * @return MolecularAssembly of the desired sphere size.
     */
    public static MolecularAssembly generateBaseSphere(MolecularAssembly molecularAssembly,
                                                       int nMolecules,
                                                       boolean original) {

        // Move molecules into the unit cell.
        molecularAssembly.moveAllIntoUnitCell();

        Crystal crystal = molecularAssembly.getCrystal();
        double asymmetricUnitVolume = crystal.volume / crystal.getNumSymOps();

        // Add wiggle room for boundary cutoffs
        double inflationFactor = 2.0;

        // Estimate a radius that will include "nMolecules".
        double radius = cbrt((3.0 / (4.0 * PI) * nMolecules * asymmetricUnitVolume)) + inflationFactor;

        // Atoms in the asymmetric unit.
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        //TODO update l, m, and n with meaningful values
        Crystal replicatesCrystal = new ReplicatesCrystal(molecularAssembly.getCrystal(), 12, 12, 12, 12*molecularAssembly.getCrystal().c);
                //ReplicatesCrystal.replicatesCrystalFactory(molecularAssembly.getCrystal(), radius * 2.0);
        molecularAssembly.setCrystal(replicatesCrystal);

        // Collect asymmetric unit atomic coordinates.
        double[] x = new double[nAtoms];
        double[] y = new double[nAtoms];
        double[] z = new double[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            x[i] = atom.getX();
            y[i] = atom.getY();
            z[i] = atom.getZ();
        }

        int nSymm = replicatesCrystal.spaceGroup.symOps.size();
        // Symmetry coordinates for each molecule in replicates crystal
        double[][] xS = new double[nSymm][nAtoms];
        double[][] yS = new double[nSymm][nAtoms];
        double[][] zS = new double[nSymm][nAtoms];
        double[][] centerMolsFrac = new double[nSymm][3]; //Fractional center of each molecule
        double[][] centerMolsCart = new double[nSymm][3]; //Cartesian center of each molecule

        // Loop over replicate crystal SymOps
        for (int iSym = 0; iSym < nSymm; iSym++) {
            SymOp symOp = replicatesCrystal.spaceGroup.getSymOp(iSym);
            // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
            replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
            // Compute center-of-mass (CoM) for Cartesian coordinates
            double[] centerOfMass = new double[3];
            int index = 0;
            if (original) {
                // PACCOM conformity -- original algorithm did not use masses.
                for (int i = 0; i < nAtoms; i++) {
                    centerOfMass[0] += xS[iSym][index];
                    centerOfMass[1] += yS[iSym][index];
                    centerOfMass[2] += zS[iSym][index++];
                }
                centerOfMass[0] /= nAtoms;
                centerOfMass[1] /= nAtoms;
                centerOfMass[2] /= nAtoms;
            } else {
                double totalMass = 0.0;
                for (Atom atom : atoms) {
                    double m = atom.getMass();
                    centerOfMass[0] += xS[iSym][index] * m;
                    centerOfMass[1] += yS[iSym][index] * m;
                    centerOfMass[2] += zS[iSym][index++] * m;
                    totalMass += m;
                }
                centerOfMass[0] /= totalMass;
                centerOfMass[1] /= totalMass;
                centerOfMass[2] /= totalMass;
            }

            // Convert CoM to fractional coordinates
            replicatesCrystal.toFractionalCoordinates(centerOfMass, centerOfMass);

            // Save CoM fractional coordinates
            centerMolsFrac[iSym] = centerOfMass;
        }

        int[] desiredMols = new int[nMolecules];

        // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
        // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
        double[] fracCenter = {0.5, 0.5, 0.5};
        double[] cartCenter = new double[3];
        replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

        // Molecule index closest to center of crystal.
        int minIndex = -1;
        double minDist = Double.MAX_VALUE;
        for (int i = 0; i < nSymm; i++) {
            // Convert each CoM to cartesian coordinates
            replicatesCrystal.toCartesianCoordinates(centerMolsFrac[i], centerMolsCart[i]);

            // Then compute euclidean distance from cartesian center of the replicates cell
            double value = distance.compute(cartCenter, centerMolsCart[i]);

            // Save index of molecule closest to replicates cell center
            if (value < minDist) {
                minDist = value;
                minIndex = i;
            }
        }

        // Hashmap with index and distance
        HashMap<Integer, Double> closestMols = new HashMap<>();

        // Loop over saved CoM coordinates to find the nMolecules closest to the marked molecule
        for (int i = 0; i < nSymm; i++) {
            double value = distance.compute(centerMolsCart[minIndex], centerMolsCart[i]);

            if (closestMols.size() < nMolecules) {
                closestMols.put(i, value);
            } else {
                double maxValue = Collections.max(closestMols.values());
                if (value < maxValue) {
                    // Remove most distant molecule
                    Iterator<Entry<Integer, Double>> it = closestMols.entrySet().iterator();
                    while (it.hasNext()) {
                        Entry<Integer, Double> item = it.next();
                        // Hopefully precision issues do not cause any issues with a direct comparison of double values.
                        if (item.getValue() == maxValue) {
                            it.remove();
                            break;
                        }
                    }
                    closestMols.put(i, value);
                }
            }
        }

        // Final collection of nMolcules closest to the marked molecule.
        int index = 0;
        for (Integer key : closestMols.keySet()) {
            desiredMols[index++] = key;
        }

        // Create expanded assemblies containing nMolecules
        MolecularAssembly expandedAssembly = new MolecularAssembly(molecularAssembly.getName());
        List<Bond> bondList = molecularAssembly.getBondList();
        ArrayList<Atom> newAtomList = new ArrayList<>();
        int atomIndex = 0;
        for (int iSym : desiredMols) {
            ArrayList<Atom> atomList = new ArrayList<>();
            // Create a new set of Atoms for each SymOp molecule
            for (int i = 0; i < nAtoms; i++) {
                Atom a = atoms[i];
                double[] xyz = new double[3];
                xyz[0] = xS[iSym][i];
                xyz[1] = yS[iSym][i];
                xyz[2] = zS[iSym][i];
                Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
                atomList.add(atom);
            }
            // Create a new set of Bonds for each SymOp molecule
            for (Bond bond : bondList) {
                Atom a1 = bond.getAtom(0);
                Atom a2 = bond.getAtom(1);
                Atom newA1 = atomList.get(a1.getIndex() - 1);
                Atom newA2 = atomList.get(a2.getIndex() - 1);
                Bond b = new Bond(newA1, newA2);
                b.setBondType(bond.getBondType());
            }
            newAtomList.addAll(atomList);
        }

        // Construct the force field for the expanded set of molecules
        ForceField forceField = molecularAssembly.getForceField();

        // Clear all periodic boundary keywords to achieve an aperiodic system.
        forceField.clearProperty("a-axis");
        forceField.clearProperty("b-axis");
        forceField.clearProperty("c-axis");
        forceField.clearProperty("alpha");
        forceField.clearProperty("beta");
        forceField.clearProperty("gamma");
        forceField.clearProperty("spacegroup");

        expandedAssembly.setForceField(forceField);
        Utilities.biochemistry(expandedAssembly, newAtomList);
        expandedAssembly.finalize(true, forceField);
        ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
        expandedAssembly.setPotential(energy);

        // This energy evaluation could be removed in the future. Is it useful as a consistency check?
        //energy.energy(true, true);

        expandedAssembly.setFile(molecularAssembly.getFile());

        return expandedAssembly;
    }

    /**
     * Cut a molecular sub-sphere from a larger sphere (molecularAssembly) that matches a smaller molecular sphere
     * (baseAssembly).
     *
     * @param asymmetricUnitAssembly Assembly of the original molecule.
     * @param inflatedAssembly       Assembly of the current sphere.
     * @param similarMolecules       Indices for desired molecules.
     * @return MolecularAssembly cut from molecularAssembly that is similar to baseAssembly.
     */
    public static MolecularAssembly cutSubSphere(MolecularAssembly asymmetricUnitAssembly,
                                                 MolecularAssembly inflatedAssembly, int[] similarMolecules) {

        int nAtoms = asymmetricUnitAssembly.getAtomArray().length;
        int nInflatedMolecules = inflatedAssembly.getMolecules().size();

        Atom[] inflatedAtoms = inflatedAssembly.getAtomArray();
        double[][] xInflated = new double[nInflatedMolecules][nAtoms];
        double[][] yInflated = new double[nInflatedMolecules][nAtoms];
        double[][] zInflated = new double[nInflatedMolecules][nAtoms];
//    // This assumes the atoms of each molecule are contiguous. This is not enforced for XYZ files.
        for (int i = 0; i < nInflatedMolecules; i++) {
            for (int j = 0; j < nAtoms; j++) {
                Atom atom = inflatedAtoms[nAtoms * i + j];
                xInflated[i][j] = atom.getX();
                yInflated[i][j] = atom.getY();
                zInflated[i][j] = atom.getZ();
            }
        }

        // Create assembly for second system based on optimal nMolecule matches
        MolecularAssembly expandedAssembly = new MolecularAssembly(asymmetricUnitAssembly.getName());
        List<Bond> bondList = asymmetricUnitAssembly.getBondList();
        ArrayList<Atom> newAtomList = new ArrayList<>();
        int atomIndex = 0;
        for (int molecule : similarMolecules) {
            ArrayList<Atom> atomList = new ArrayList<>();
            // Create a new set of Atoms for each SymOp molecule
            for (int i = 0; i < nAtoms; i++) {
                Atom a = inflatedAtoms[molecule * nAtoms + i];
                double[] xyz = new double[3];
                xyz[0] = xInflated[molecule][i];
                xyz[1] = yInflated[molecule][i];
                xyz[2] = zInflated[molecule][i];
                Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
                atomList.add(atom);
            }
            // Create a new set of Bonds for each SymOp molecule
            for (Bond bond : bondList) {
                Atom a1 = bond.getAtom(0);
                Atom a2 = bond.getAtom(1);
                Atom newA1 = atomList.get(a1.getIndex() - 1);
                Atom newA2 = atomList.get(a2.getIndex() - 1);
                Bond b = new Bond(newA1, newA2);
                b.setBondType(bond.getBondType());
            }
            newAtomList.addAll(atomList);
        }

        // Construct the force field for the expanded set of molecules
        ForceField forceField = inflatedAssembly.getForceField();

        // Clear all periodic boundary keywords to achieve an aperiodic system.
        forceField.clearProperty("a-axis");
        forceField.clearProperty("b-axis");
        forceField.clearProperty("c-axis");
        forceField.clearProperty("alpha");
        forceField.clearProperty("beta");
        forceField.clearProperty("gamma");
        forceField.clearProperty("spacegroup");

        expandedAssembly.setForceField(forceField);
        Utilities.biochemistry(expandedAssembly, newAtomList);
        expandedAssembly.finalize(true, forceField);
        ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
        expandedAssembly.setPotential(energy);

        // This energy evaluation could be removed in the future. Is it useful as a consistency check?
        //energy.energy(true, true);

        expandedAssembly.setFile(inflatedAssembly.getFile());
        return expandedAssembly;
    }

    /**
     * Determine the index for the central molecule based on center of mass for an assembly.
     *
     * @param assembly   Assembly of which the center is desired.
     * @param nMolecules Number of central molecules to find.
     * @param original   Is the script using the original PACCOM algorithm?
     * @return int index of the molecule closest to the center of mass.
     */
    public static int[] centerMoleculesIndices(MolecularAssembly assembly, int nMolecules, boolean original) {
        int nMolAssembly = assembly.getMolecules().size();
        int nAtoms = assembly.getAtomArray().length;
        double[] centerOfMass = new double[3];
        double totalMass = 0.0;

        //Determine center of assembly
        if (original) {
            for (Atom atom : assembly.getAtomArray()) {
                centerOfMass[0] += atom.getX();
                centerOfMass[1] += atom.getY();
                centerOfMass[2] += atom.getZ();
            }
            centerOfMass[0] /= nAtoms;
            centerOfMass[1] /= nAtoms;
            centerOfMass[2] /= nAtoms;
        } else {
            for (Atom atom : assembly.getAtomArray()) {
                double m = atom.getMass();
                centerOfMass[0] += atom.getX() * m;
                centerOfMass[1] += atom.getY() * m;
                centerOfMass[2] += atom.getZ() * m;
                totalMass += m;
            }
            centerOfMass[0] /= totalMass;
            centerOfMass[1] /= totalMass;
            centerOfMass[2] /= totalMass;
        }

        logger.info(" Center of shell:");
        logger.info(format(" %.4f %.4f %.4f", centerOfMass[0], centerOfMass[1], centerOfMass[2]));
        // Find closest nMolecules to center.
        HashMap<Integer, Double> targetClosestMols = new HashMap<>();
        for (int j = 0; j < nMolAssembly; j++) {
            double[] nextCenter = assembly.getMolecules().get(j).getCenter(false);

            double value = distance.compute(centerOfMass, nextCenter);
            if (targetClosestMols.size() < nMolecules) {
                targetClosestMols.put(j, value);
            } else {
                double maxValue = Collections.max(targetClosestMols.values());
                if (value < maxValue) {
                    Iterator<Map.Entry<Integer, Double>> it = targetClosestMols.entrySet().iterator();
                    while (it.hasNext()) {
                        Map.Entry<Integer, Double> item = it.next();
                        if (item.getValue() == maxValue) {
                            it.remove();
                            break;
                        }
                    }
                    targetClosestMols.put(j, value);
                }
            }
        }
        targetClosestMols = PACCOMFunctions.sortHashMapByValue(targetClosestMols);
        for (Integer key : targetClosestMols.keySet()) {
            logger.finest(format("Key: %d Dist: %f", key, targetClosestMols.get(key)));
        }
        Integer[] values = targetClosestMols.keySet().toArray(targetClosestMols.keySet().toArray(new Integer[nMolecules]));
        // Convert Integer to int
        int[] returnValues = new int[nMolecules];
        int index = 0;
        for (Integer value : values) {
            returnValues[index++] = value;
        }

        return returnValues;
    }

    /**
     * Determine the indices for the molecules near target molecule in an assembly.
     *
     * @param assembly   Assembly of which the center is desired.
     * @param targetMol Index of desired center molecule.
     * @param nMolecules Number of central molecules to find.
     * @param original   Is the script using the original PACCOM algorithm?
     * @return int index of the molecule closest to the center of mass.
     */
    public static int[] moleculeIndicesNear(MolecularAssembly assembly, int targetMol, int currAssembly, int[][] minIndicies, int nMolecules, boolean original) {
        int nMolAssembly = assembly.getMolecules().size();
        int nAtoms = assembly.getAtomArray().length;
        double[] centerOfMass = assembly.getMolecules().get(minIndicies[currAssembly][targetMol]).getCenter(!original);

        // Find closest nMolecules to center.
        HashMap<Integer, Double> targetClosestMols = new HashMap<>();
        for (int j = 0; j < nMolAssembly; j++) {
            double[] nextCenter = assembly.getMolecules().get(j).getCenter(!original);

            double value = distance.compute(centerOfMass, nextCenter);
            if (targetClosestMols.size() < nMolecules) {
                targetClosestMols.put(j, value);
            } else {
                double maxValue = Collections.max(targetClosestMols.values());
                if (value < maxValue) {
                    Iterator<Map.Entry<Integer, Double>> it = targetClosestMols.entrySet().iterator();
                    while (it.hasNext()) {
                        Map.Entry<Integer, Double> item = it.next();
                        if (item.getValue() == maxValue) {
                            it.remove();
                            break;
                        }
                    }
                    targetClosestMols.put(j, value);
                }
            }
        }
        targetClosestMols = PACCOMFunctions.sortHashMapByValue(targetClosestMols);
        for (Integer key : targetClosestMols.keySet()) {
            logger.finest(format("Key: %d Dist: %f", key, targetClosestMols.get(key)));
        }
        Integer[] values = targetClosestMols.keySet().toArray(targetClosestMols.keySet().toArray(new Integer[nMolecules]));
        // Convert Integer to int
        int[] returnValues = new int[nMolecules];
        int index = 0;
        for (Integer value : values) {
            returnValues[index++] = value;
        }

        return returnValues;
    }

    /**
     * Sort the given HashMap (Integer, Double) by values (Double).
     *
     * @param hm Hashmap to be sorted.
     * @return HashMap(Integer, Double) that has been sorted.
     */
    public static HashMap<Integer, Double> sortHashMapByValue(HashMap<Integer, Double> hm) {
        // Create a list from elements of HashMap
        List<Entry<Integer, Double>> list = new LinkedList<>(hm.entrySet());

        // Sort the list
        list.sort(Entry.comparingByValue());

        // put data from sorted list to hashMap
        HashMap<Integer, Double> temp = new LinkedHashMap<>();
        for (Entry<Integer, Double> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }

    /**
     * Reorients the given atoms so n2 and n3 are on the X axis and X-Y plane. n1 should be at origin.
     *
     * @param n1 Atom located at the origin
     * @param n2 Atom to go to X axis
     * @param n3 Atom to go to X-Y plane
     * @return double[] New coordinates for the atoms
     * @author Okimasa OKADA
     * @author Aaron J. Nessler
     */
    static double[] standardOrientation(Atom n1, Atom n2, Atom n3) {
        double[] orientedCoords = new double[9];
        double p1n2 = n2.getX();
        double q1n2 = n2.getY();
        double r1n2 = n2.getZ();
        double p2n2;
        double q2n2;

        double p1n3 = n3.getX();
        double q1n3 = n3.getY();
        double r1n3 = n3.getZ();
        double p2n3;
        double q2n3;

        logger.fine(format("START: N1: %16.8f %16.8f %16.8f", n1.getX(), n1.getY(), n1.getZ()));
        logger.fine(format("START: N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
        logger.fine(format("START: N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));
        // Calculation of sigma, phai, and cita angles needed to get specified atoms to desired loci
        double cita0 = acos(p1n2 / sqrt(pow(p1n2, 2) + pow(q1n2, 2)));
        double phai0 = acos(sqrt(pow(p1n2, 2) + pow(q1n2, 2)) / sqrt(pow(p1n2, 2) + pow(q1n2, 2) + pow(r1n2, 2)));
        if (q1n2 < 0.0) {
            cita0 = -cita0;
        }

        logger.fine(format("cita: %16.8f", cita0));
        logger.fine(format("cos(cita): %16.8f", cos(cita0)));
        logger.fine(format("sin(cita): %16.8f", sin(cita0)));

        p2n2 = p1n2 * cos(cita0) + q1n2 * sin(cita0);
        q1n2 = -p1n2 * sin(cita0) + q1n2 * cos(cita0);
        p2n3 = p1n3 * cos(cita0) + q1n3 * sin(cita0);
        q1n3 = -p1n3 * sin(cita0) + q1n3 * cos(cita0);

        p1n2 = p2n2;
        p1n3 = p2n3;

        logger.fine(format("Step 1 N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
        logger.fine(format("Step 1 N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

        if (r1n2 > 0.0) {
            phai0 = -phai0;
        }

        logger.fine(format("phai: %16.8f", phai0));

        p2n2 = p1n2 * cos(phai0) - r1n2 * sin(phai0);
        r1n2 = p1n2 * sin(phai0) + r1n2 * cos(phai0);
        p2n3 = p1n3 * cos(phai0) - r1n3 * sin(phai0);
        r1n3 = p1n3 * sin(phai0) + r1n3 * cos(phai0);

        p1n2 = p2n2;
        p1n3 = p2n3;

        logger.fine(format("Step 2 N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
        logger.fine(format("Step 2 N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

        double sigma0 = acos(q1n3 / sqrt(q1n3 * q1n3 + r1n3 * r1n3));
        if (r1n3 < 0.0) {
            sigma0 = -sigma0;
        }

        logger.fine(format("sigma: %16.8f", sigma0));

        q2n2 = q1n2 * cos(sigma0) + r1n2 * sin(sigma0);
        r1n2 = -q1n2 * sin(sigma0) + r1n2 * cos(sigma0);
        q2n3 = q1n3 * cos(sigma0) + r1n3 * sin(sigma0);
        r1n3 = -q1n3 * sin(sigma0) + r1n3 * cos(sigma0);

        q1n2 = q2n2;
        q1n3 = q2n3;

        logger.fine(format("DONE: N1: %16.8f %16.8f %16.8f", n1.getX(), n1.getY(), n1.getZ()));
        logger.fine(format("DONE N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
        logger.fine(format("DONE N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

        orientedCoords[0] = n1.getX();
        orientedCoords[1] = n1.getY();
        orientedCoords[2] = n1.getZ();
        orientedCoords[3] = p1n2;
        orientedCoords[4] = q1n2;
        orientedCoords[5] = r1n2;
        orientedCoords[6] = p1n3;
        orientedCoords[7] = q1n3;
        orientedCoords[8] = r1n3;

        return orientedCoords;
    }

    /**
     * Save the current assembly as a PDB file.
     *
     * @param currentAssembly Current assembly of which to save a PDB.
     * @param fileName        Prefix for the file being saved.
     * @param i               Number of file being saved.
     */
    public static void saveAssemblyPDB(MolecularAssembly currentAssembly, String fileName, int i) {
        String name = Integer.toString(i);
        File saveLocationPDB = new File(fileName + "_" + name + ".pdb");
        //Save aperiodic system of n_mol closest atoms for visualization
        logger.info(" Saving " + saveLocationPDB.getName());
        PDBFilter pdbfilter = new PDBFilter(saveLocationPDB, currentAssembly, null, null);
        pdbfilter.writeFile(saveLocationPDB, false);
    }
}