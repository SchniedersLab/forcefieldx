package ffx.potential.utils;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSNode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import org.apache.commons.io.FilenameUtils;

import java.util.logging.Logger;

import java.io.File;
import java.util.*;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.*;
import static org.apache.commons.math3.util.FastMath.pow;

public class PACCOMFunctions {

    static MolecularAssembly GenerateBaseShell(MolecularAssembly molecularAssembly, int n_mol,
                                               boolean original) {
        molecularAssembly.moveAllIntoUnitCell();

        Crystal crystal = molecularAssembly.getCrystal();
        double asymmetricUnitVolume = crystal.volume / crystal.getNumSymOps();
        double inflationFactor = 2.0; //Add wiggle room for boundary cutoffs
        double radius = cbrt((3.0 / (4.0 * PI) * n_mol * asymmetricUnitVolume)) + inflationFactor;
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length; //# of atoms in one molecule

        Crystal replicatesCrystal =
                ReplicatesCrystal.replicatesCrystalFactory(molecularAssembly.getCrystal(), radius * 2.0);
        molecularAssembly.setCrystal(replicatesCrystal);

        double[] x = new double[nAtoms];
        double[] y = new double[nAtoms];
        double[] z = new double[nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            x[i] = atom.getX();
            y[i] = atom.getY();
            z[i] = atom.getZ();
        }

        int numSymOps = replicatesCrystal.spaceGroup.symOps.size();
        //Symmetry coordinates for each molecule in replicates crystal
        double[][] xS = new double[numSymOps][nAtoms];
        double[][] yS = new double[numSymOps][nAtoms];
        double[][] zS = new double[numSymOps][nAtoms];
        double[][] centerMolsFrac = new double[numSymOps][3]; //Fractional center of each molecule
        double[][] centerMolsCart = new double[numSymOps][3]; //Cartesian center of each molecule
        // Loop over replicate crystal SymOps
        int numWarn = 0; //Number of molecules outside replicates cell
        boolean moleOutofCell = false;

        for (int iSym = 0; iSym < numSymOps; iSym++) {
            SymOp symOp = replicatesCrystal.spaceGroup.getSymOp(iSym);
            // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
            replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
            // Compute center-of-mass (CoM) for Cartesian coordinates
            double[] centerOfMass = new double[3];
            double totalMass = 0.0;
            int index = 0;
            if (original) { //PACCOM conformity. Original PACCOM did not use masses.
                for (Atom atom : atoms) {
                    centerOfMass[0] += xS[iSym][index];
                    centerOfMass[1] += yS[iSym][index];
                    centerOfMass[2] += zS[iSym][index++];
                }
                centerOfMass[0] /= nAtoms;
                centerOfMass[1] /= nAtoms;
                centerOfMass[2] /= nAtoms;
            } else {
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
            // Check that CoM fractional coordinates are within replicates cell.
            for (int i = 0; i < centerOfMass.length; i++) {
                if (centerOfMass[i] <= 0.0 || centerOfMass[i] >= 1.0) {
                    moleOutofCell = true;
                    numWarn++;
                    i = centerOfMass.length; //Exit loop if one of the coordinates are out of the cell
                }
            }
            // ToDo: If CoM fractional coordinates are outside the replicates cell, move the molecule back inside
            // Save CoM fractional coordinates
            centerMolsFrac[iSym] = centerOfMass;
        }
        int[] desiredMols = new int[n_mol];
        // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
        // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
        double[] fracCenter = new double[3];
        double[] cartCenter = new double[3];
        fracCenter[0] = 0.5;
        fracCenter[1] = 0.5;
        fracCenter[2] = 0.5;
        replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

        // Convert each CoM to cartesian coordinates
        // Then compute euclidean distance from cartesian center of the replicates cell
        // Save index of molecule closest to replicates cell center
        int minIndex = -1; // Molecule index closest to center of crystal.
        double minDist = Double.MAX_VALUE;
        for (int i = 0; i < numSymOps; i++) {
            replicatesCrystal.toCartesianCoordinates(centerMolsFrac[i], centerMolsCart[i]);
            double value = sqrt(pow(cartCenter[0] - centerMolsCart[i][0], 2)
                    + pow(cartCenter[1] - centerMolsCart[i][1], 2)
                    + pow(cartCenter[2] - centerMolsCart[i][2], 2));
            if (value < minDist) {
                minDist = value;
                minIndex = i;
            }
        }
        //Hashmap with index and distance
        HashMap<Integer, Double> closestMols = new HashMap<Integer, Double>();
        // Loop over saved CoM coordinates to find the n_mol closest molecules to the marked molecule
        for (int i = 0; i < numSymOps; i++) {
            double value = sqrt(pow(centerMolsCart[minIndex][0] - centerMolsCart[i][0], 2)
                    + pow(centerMolsCart[minIndex][1] - centerMolsCart[i][1], 2)
                    + pow(centerMolsCart[minIndex][2] - centerMolsCart[i][2], 2));
            if (closestMols.size() < n_mol) {
                closestMols.put(i, value);
            } else {
                double maxValue = Collections.max(closestMols.values());
                if (value < maxValue) {
                    Iterator<Map.Entry<Integer, Double>> it = closestMols.entrySet().iterator();
                    while (it.hasNext() && closestMols.size() == n_mol) {
                        Map.Entry<Integer, Double> item = it.next();
                        if (item.getValue() == maxValue) {
                            it.remove();
                        }
                    }
                    closestMols.put(i, value);
                }
            }
        }
        int index = 0;
        for (Integer key : closestMols.keySet()) {
            desiredMols[index++] = key;
        }

        // CREATE EXPANDED ASSEMBLIES CONTAINING n_mol MOLECULES
        MolecularAssembly expandedAssembly = new MolecularAssembly(molecularAssembly.getName());
        List<Bond> bondList = molecularAssembly.getBondList();
        ArrayList<Bond> newBondList = new ArrayList<>();
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
                newBondList.add(b);
            }
            newAtomList.addAll(atomList);
        }
        // Construct the force field for the expanded set of molecules
        ForceField forceField = molecularAssembly.getForceField();
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
        energy.energy(true, true);

        expandedAssembly.setFile(molecularAssembly.getFile());

        return expandedAssembly;
    }

    static MolecularAssembly CutSubShell(MolecularAssembly originalAssembly, MolecularAssembly baseAssembly,
                                         MolecularAssembly molecularAssembly, boolean original) {
        molecularAssembly.moveAllIntoUnitCell();

        Atom[] baseAtoms = baseAssembly.getAtomArray();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        int nBaseAtoms = baseAtoms.length;
        int nAtomsPerMolecule = molecularAssembly.getMolecules().get(0).getAtomList().size();
        int nBaseMolecules = baseAssembly.getMolecules().size();
        int nMolecules = molecularAssembly.getMolecules().size();

        double[][] xBase = new double[nBaseMolecules][nAtomsPerMolecule];
        double[][] yBase = new double[nBaseMolecules][nAtomsPerMolecule];
        double[][] zBase = new double[nBaseMolecules][nAtomsPerMolecule];
        for (int i = 0; i < nBaseMolecules; i++) {
            for (int j = 0; j < nAtomsPerMolecule; j++) {
                Atom atom = baseAtoms[nAtomsPerMolecule * i + j];
                xBase[i][j] = atom.getX();
                yBase[i][j] = atom.getY();
                zBase[i][j] = atom.getZ();
            }
        }

        double[][] x2 = new double[nMolecules][nAtomsPerMolecule];
        double[][] y2 = new double[nMolecules][nAtomsPerMolecule];
        double[][] z2 = new double[nMolecules][nAtomsPerMolecule];
        for (int i = 0; i < nMolecules; i++) {
            for (int j = 0; j < nAtomsPerMolecule; j++) {
                Atom atom = atoms[nAtomsPerMolecule * i + j];
                x2[i][j] = atom.getX();
                y2[i][j] = atom.getY();
                z2[i][j] = atom.getZ();
            }
        }

        double[][] centerOfMols = new double[nMolecules][3];
        // Compute center-of-mass (CoM) for Cartesian coordinates


        if (original) { //PACCOM conformity. Original PACCOM did not use masses.
            for (int molIndex = 0; molIndex < nMolecules; molIndex++) {
                double[] centerOfMass = new double[3];
                double totalMass = 0.0;
                for (int atomIndex = 0; atomIndex < nAtomsPerMolecule; atomIndex++) {
                    centerOfMass[0] += x2[molIndex][atomIndex];
                    centerOfMass[1] += y2[molIndex][atomIndex];
                    centerOfMass[2] += z2[molIndex][atomIndex];
                }
                centerOfMass[0] /= atoms.length;
                centerOfMass[1] /= atoms.length;
                centerOfMass[2] /= atoms.length;
                centerOfMols[molIndex] = centerOfMass;
            }
        } else {
            for (int molIndex = 0; molIndex < nMolecules; molIndex++) {
                double[] centerOfMass = new double[3];
                double totalMass = 0.0;
                for (int atomIndex = 0; atomIndex < nAtomsPerMolecule; atomIndex++) {
                    double m = atoms[molIndex * nAtomsPerMolecule + atomIndex].getMass();
                    centerOfMass[0] += x2[molIndex][atomIndex] * m;
                    centerOfMass[1] += y2[molIndex][atomIndex] * m;
                    centerOfMass[2] += z2[molIndex][atomIndex] * m;
                    totalMass += m;
                }
                centerOfMass[0] /= totalMass;
                centerOfMass[1] /= totalMass;
                centerOfMass[2] /= totalMass;
                centerOfMols[molIndex] = centerOfMass;
            }
            // ToDo: If CoM fractional coordinates are outside the replicates cell, move the molecule back inside
            // Save CoM fractional coordinates
        }

        int[] desiredMols = new int[nBaseMolecules];
        Arrays.fill(desiredMols, -1);
        int index = 0;
        for (MSNode baseMolecule : baseAssembly.getMolecules()) {
            double[] baseCenter = baseMolecule.getCenter(true);
            int key = -1;
            double minValue = Double.MAX_VALUE;
            for (int i = 0; i < nMolecules; i++) {
                boolean duplicate = false;
                double value = sqrt(pow(baseCenter[0] - centerOfMols[i][0], 2)
                        + pow(baseCenter[1] - centerOfMols[i][1], 2)
                        + pow(baseCenter[2] - centerOfMols[i][2], 2));
                if (value < minValue) {
                    for (int prevKey : desiredMols) {
                        if (prevKey == i) {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate) {
                        minValue = value;
                        key = i;
                    }
                }
            }
            desiredMols[index++] = key;
        }

        // CREATE EXPANDED ASSEMBLIES CONTAINING n_mol MOLECULES
        MolecularAssembly expandedAssembly = new MolecularAssembly(originalAssembly.getName());
        List<Bond> bondList = originalAssembly.getBondList();
        ArrayList<Bond> newBondList = new ArrayList<>();
        ArrayList<Atom> newAtomList = new ArrayList<>();
        int atomIndex = 0;
        for (int molecule : desiredMols) {
            ArrayList<Atom> atomList = new ArrayList<>();
            // Create a new set of Atoms for each SymOp molecule
            for (int i = 0; i < nAtomsPerMolecule; i++) {
                Atom a = atoms[molecule*nAtomsPerMolecule + i];
                double[] xyz = new double[3];
                xyz[0] = x2[molecule][i];
                xyz[1] = y2[molecule][i];
                xyz[2] = z2[molecule][i];
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
                newBondList.add(b);
            }
            newAtomList.addAll(atomList);
        }
        // Construct the force field for the expanded set of molecules
        ForceField forceField = molecularAssembly.getForceField();
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
        energy.energy(true, true);

        expandedAssembly.setFile(molecularAssembly.getFile());

        return expandedAssembly;
    }

    static int CenterMoleculeIndex(MolecularAssembly assembly){
        int minIndex = -1; // Molecule index closest to center of crystal.
        double minDist = Double.MAX_VALUE;
        int nMolecules = assembly.getMolecules().size();
        double[] centerOfMass = new double[3];
        double totalMass = 0.0;

        for(Atom atom: assembly.getAtomArray()){
            double m = atom.getMass();
            centerOfMass[0]+=atom.getX() * m;
            centerOfMass[1]+=atom.getY() * m;
            centerOfMass[2]+=atom.getZ() * m;
            totalMass += m;
        }
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;

        assembly.getMolecules().get(0).getCenter(true);
        for (int i=0; i<nMolecules;i++) {
            double[] centerMol = assembly.getMolecules().get(i).getCenter(true);
            double value = sqrt(pow(centerOfMass[0] - centerMol[0], 2)
                    + pow(centerOfMass[1] - centerMol[1], 2)
                    + pow(centerOfMass[2] - centerMol[2], 2));
            if (value < minDist) {
                minDist = value;
                minIndex = i;
            }
        }
        return minIndex;
    }

    /**
     * Sort the given HashMap by values (Doubles).
     *
     * @param hm Hashmap to be sorted.
     */
    static HashMap<Integer, Double> SortHashMapByValue(HashMap<Integer, Double> hm) {
        // Create a list from elements of HashMap
        List<Map.Entry<Integer, Double>> list =
                new LinkedList<Map.Entry<Integer, Double>>(hm.entrySet());
        // Sort the list
        list.sort(new Comparator<Map.Entry<Integer, Double>>() {
            public int compare(Map.Entry<Integer, Double> o1,
                               Map.Entry<Integer, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });

        // put data from sorted list to hashMap
        HashMap<Integer, Double> temp = new LinkedHashMap<Integer, Double>();
        for (Map.Entry<Integer, Double> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }

    /**
     * Save the current assembly as a PDB file.
     *
     * @param logger          Logger used to display information.
     * @param currentAssembly Current assembly of which to save a PDB.
     * @param saveLocationPDB A pre-made file to save the PDB.
     */
    static void SaveAssemblyPDB(Logger logger, MolecularAssembly currentAssembly, File saveLocationPDB) {
        //Save aperiodic system of n_mol closest atoms for visualization
        logger.info(" Saving " + saveLocationPDB.getName());
        PDBFilter pdbfilter = new PDBFilter(saveLocationPDB, currentAssembly, null, null);
        pdbfilter.writeFile(saveLocationPDB, false);
    }
}
