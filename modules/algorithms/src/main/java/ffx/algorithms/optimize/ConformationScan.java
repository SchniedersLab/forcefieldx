package ffx.algorithms.optimize;

import com.google.common.collect.MinMaxPriorityQueue;
import ffx.algorithms.AlgorithmListener;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.parameters.BondType;
import ffx.potential.parsers.XYZFilter;

import java.io.File;
import java.util.AbstractQueue;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.logging.Logger;

import static ffx.potential.utils.Superpose.applyRotation;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.*;

/**
 * This class is for a configuration optimization of two small monomers. The search
 * consists of aligning heavy (except carbon) and hydrogen atoms at a h-bond distance
 * with their center of masses behind them all on the z-axis in hopes that a global
 * conformational minima will be found by allowing for strong h-bonds to occur. After
 * aligning, an optional static torsion scan and NONE,DIRECT,MUTUAL minimizations are
 * performed. A binding energy can be calculated with this by:
 *
 * ddEbind = dEbind (of a dimer)
 *         - dEbind (of one monomer w/ itself)
 *         - dEbind (of the other monomer w/ itself)
 */
public class ConformationScan {

    private static final Logger logger = Logger.getLogger(ConformationScan.class.getName());

    private final MolecularAssembly mola;
    private final ForceFieldEnergy forceFieldEnergy;
    private double[] x;
    private AlgorithmListener algorithmListener;

    private final Molecule m1;
    private final Molecule m2;
    private final Atom[] m1Atoms;
    private final Atom[] m2Atoms;
    private final ArrayList<Atom> m1TargetAtoms;
    private final ArrayList<Atom> m2TargetAtoms;

    private AbstractQueue<StateContainer> statesQueue;

    private final double eps;
    private final int maxIter;
    private final double hBondDist;
    private final double flatBottomRadius;
    private final boolean tScan;
    private final boolean excludeH;
    private final boolean minimize;

    private double m1MinEnergy;
    private double m2MinEnergy;
    private double totalMonomerMinimizedEnergy;

    private double minimumEnergy;
    private double averageEnergy;
    private double averageEnergyNoOutlier;
    private double stdOfEnergies;
    private double stdOfEnergiesNoOutlier;

    public ConformationScan(
            MolecularAssembly molecularAssembly,
            Molecule monomerOne,
            Molecule monomerTwo,
            double minimizeEps,
            int minimizeMaxIterations,
            double hydrogenBondDistance,
            double flatBottomRadius,
            boolean staticTorsionScan,
            boolean excludeExtraHydrogen,
            boolean minimize
    )
    {
        mola = molecularAssembly;
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        m1 = monomerOne;
        m2 = monomerTwo;
        eps = minimizeEps;
        maxIter = minimizeMaxIterations;
        hBondDist = hydrogenBondDistance;
        this.flatBottomRadius = flatBottomRadius;
        tScan = staticTorsionScan;
        excludeH = excludeExtraHydrogen;
        this.minimize = minimize;
        m1Atoms = m1.getAtomList().toArray(new Atom[0]);
        m2Atoms = m2.getAtomList().toArray(new Atom[0]);
        m1TargetAtoms = new ArrayList<>();
        m2TargetAtoms = new ArrayList<>();
        setTargetAtoms(mola.getAtomArray());
        statesQueue = MinMaxPriorityQueue.maximumSize(m1TargetAtoms.size() * m2TargetAtoms.size()).create();
        if(!minimize) { // If minimization is not performed, then the queue will be in order of the loop
            statesQueue = new ArrayBlockingQueue<>(m1TargetAtoms.size() * m2TargetAtoms.size() + 1);
        }
    }

    public void scan(){
        minimizeEachMolecule(minimize);
        double[] zAxis = new double[]{0,0,1};
        // Loop through interactions between the two molecules --> Not necessarily symmetric but close?
        int loopCounter = 1;
        for(Atom a: m1TargetAtoms){
            alignMoleculeCOMtoAtomVecWithAxis(m1, a, zAxis, m1Atoms);
            // Align with opposite side of axis --> Reset after loop
            zAxis[2] = -1;
            for(Atom b: m2TargetAtoms){
                logger.info("\n ----- Trial " + loopCounter + " out of " +
                        (m2TargetAtoms.size() * m1TargetAtoms.size()) + " -----");
                loopCounter++;
                alignMoleculeCOMtoAtomVecWithAxis(m2, b, zAxis, m2Atoms);
                // Move the second molecule to the perfect h-bond distance
                double[] hBondVector = new double[]{0, 0, a.getZ() - b.getZ() + hBondDist};
                for (Atom m2Atom : m2Atoms) {
                    m2Atom.move(hBondVector);
                }
                // Minimize the energy of the system subject to a harmonic restraint on the distance
                // between the two atoms. Keep the state if minimization works.
                try {
                    if(minimize) {
                        minimizeSystem(a, b);
                    }
                    forceFieldEnergy.getCoordinates(x);
                    double e = forceFieldEnergy.energy(x, false) - totalMonomerMinimizedEnergy;
                    logger.info(" Binding energy of trial " + (loopCounter-1) + ": " + e);
                    statesQueue.add(new StateContainer(new AssemblyState(mola), e));
                } catch (Exception ignored) {
                    logger.warning(" Minimization failed. No state will be saved.");
                    //e.printStackTrace()
                }
            }
            // Reset z-axis for mol 1 alignment
            zAxis[2] = 1;
        }
        logger.info("\n ------------------------- End of Trials -------------------------");
        if(statesQueue.peek() != null) {
            minimumEnergy = statesQueue.peek().getEnergy();
        } else{
            logger.warning(" No states were saved. No minimum energy found.");
            minimumEnergy = Double.NaN;
            return;
        }
        calculateMeanStd();
    }
    public void logAllEnergyInformation(){
        logger.info(format(" Minimum energy of monomer 1:                 %12.5f", m1MinEnergy));
        logger.info(format(" Minimum energy of monomer 2:                 %12.5f", m2MinEnergy));
        logger.info(format(" Sum of minimum monomers:                     %12.5f", totalMonomerMinimizedEnergy));
        logger.info(format(" Minimum binding energy of system:            %12.5f", minimumEnergy));
        logger.info(format(" Average binding energy:                      %12.5f +- %6.3f", averageEnergy, stdOfEnergies));
        logger.info(format(" Average binding energy (no outlier):         %12.5f +- %6.3f", averageEnergyNoOutlier, stdOfEnergiesNoOutlier));
    }
    public void writeStructuresToXYZ(File outputFile){
        XYZFilter xyzFilter = new XYZFilter(outputFile, mola, mola.getForceField(), mola.getForceField().getProperties());
        logger.info("\n Writing structures to " + outputFile.getAbsolutePath());
        int count = 1;
        StateContainer[] states = statesQueue.toArray(new StateContainer[0]);
        for(StateContainer s: states){
            AssemblyState state = s.getState();
            double energy = s.getEnergy();
            logger.info(" Writing configuration #" + count + " with energy: " + energy);
            state.revertState();
            xyzFilter.writeFile(outputFile, count != 1);
            count++;
        }
    }

    public static ArrayList<Double> logBindingEnergyCalculation(ConformationScan m1, ConformationScan m2, ConformationScan dimer){
        ArrayList<Double> bindingEnergies = new ArrayList<>();
        logger.info("\n ------------------------- Binding Energy -------------------------");
        logger.info(" ddE (of binding) = dE (of coformer dimer binding) - dE (of both monomer dimer bindings summed)");
        logger.info("\n Calculation using minimized energies: ");
        logger.info(format(" Minimal Structure    ddE = %8.3f - %8.3f = %12.5f",
                dimer.getMinimumEnergy(),
                m1.getMinimumEnergy() + m2.getMinimumEnergy(),
                dimer.getMinimumEnergy() - m1.getMinimumEnergy() - m2.getMinimumEnergy()));
        bindingEnergies.add(dimer.getMinimumEnergy() - m1.getMinimumEnergy() - m2.getMinimumEnergy());
        logger.info("\n Calculation using average energies: ");
        logger.info(format(" Average              ddE = %8.3f - %8.3f = %12.5f",
                dimer.getAverageEnergy(),
                m1.getAverageEnergy() + m2.getAverageEnergy(),
                dimer.getAverageEnergy() - m1.getAverageEnergy() - m2.getAverageEnergy()));
        bindingEnergies.add(dimer.getAverageEnergy() - m1.getAverageEnergy() - m2.getAverageEnergy());
        logger.info("\n Calculation using average energies (no outliers): ");
        logger.info(format(" Average (no outlier) ddE = %8.3f - %8.3f = %12.5f",
                dimer.getAverageEnergyNoOutlier(),
                m1.getAverageEnergyNoOutlier() + m2.getAverageEnergyNoOutlier(),
                dimer.getAverageEnergyNoOutlier() - m1.getAverageEnergyNoOutlier() - m2.getAverageEnergyNoOutlier()));
        bindingEnergies.add(dimer.getAverageEnergyNoOutlier() - m1.getAverageEnergyNoOutlier() - m2.getAverageEnergyNoOutlier());
        logger.info("\n N.B. -- Monomer energies are calculated using the minimized structures in all cases.");
        logger.info(" ------------------------- End of Binding Energy -------------------------\n");
        return bindingEnergies;
    }

    public ArrayList<AssemblyState> getStates(){
        ArrayList<AssemblyState> states = new ArrayList<>();
        for(StateContainer s: statesQueue){
            states.add(s.getState());
        }
        return states;
    }
    public ArrayList<AssemblyState> getStatesWithinEnergy(double energy){
        ArrayList<AssemblyState> states = new ArrayList<>();
        StateContainer minState = statesQueue.peek();
        if(minState != null){
            double minEnergy = minState.getEnergy();
            for(StateContainer s: statesQueue){
                if(s.getEnergy() - minEnergy < energy){
                    states.add(s.getState());
                }
            }
        }
        return states;
    }
    public ArrayList<AssemblyState> getStatesFilteredByRMSD(double rmsdCutoff){return null;}
    public ArrayList<AssemblyState> getStatesAroundAverage(double averageE, double radialCutoff){return null;}
    public ArrayList<Double> getEnergies(){
        ArrayList<Double> energies = new ArrayList<>();
        for(StateContainer s: statesQueue){
            energies.add(s.getEnergy());
        }
        return energies;
    }
    public ArrayList<Double> getEnergiesWithinEnergy(double energy){
        ArrayList<Double> energies = new ArrayList<>();
        StateContainer minState = statesQueue.peek();
        if(minState != null){
            double minEnergy = minState.getEnergy();
            for(StateContainer s: statesQueue){
                if(s.getEnergy() - minEnergy < energy){
                    energies.add(s.getEnergy());
                }
            }
        }
        return energies;
    }
    public double getM1MinEnergy() {return m1MinEnergy;}
    public double getM2MinEnergy() {return m2MinEnergy;}
    public double getTotalMonomerMinimizedEnergy() {return totalMonomerMinimizedEnergy;}
    public double getMinimumEnergy() {return minimumEnergy;}
    public double getAverageEnergy() {return averageEnergy;}
    public double getAverageEnergyNoOutlier() {return averageEnergyNoOutlier;}
    public double getStdOfEnergies() {return stdOfEnergies;}
    public double getStdOfEnergiesNoOutlier() {return stdOfEnergiesNoOutlier;}

    private void setTargetAtoms(Atom[] atoms){
        for(Atom a: atoms){
            if(a.getAtomType().atomicNumber == 7|| a.getAtomType().atomicNumber == 8
                    || a.getAtomType().atomicNumber == 9 || a.getAtomType().atomicNumber == 15
                    || a.getAtomType().atomicNumber == 16 || a.getAtomType().atomicNumber == 17){ // N,O,F,P,S,Cl
                if(a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                    m1TargetAtoms.add(a);
                } else{
                    m2TargetAtoms.add(a);
                }

                // Searching for bonded h's only if we are excluding H's that aren't bonded to electronegative atoms
                if(excludeH) {
                    for (Bond b : a.getBonds()) {
                        int num = b.get1_2(a).getAtomType().atomicNumber;
                        if (num == 1) {
                            if (a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                                m1TargetAtoms.add(a);
                            } else {
                                m2TargetAtoms.add(a);
                            }
                        }
                    }
                }
            }
            else if(a.getAtomType().atomicNumber == 1 && !excludeH){ // Add all H's
                if(a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                    m1TargetAtoms.add(a);
                } else{
                    m2TargetAtoms.add(a);
                }
            }
        }

    }
    private void minimizeEachMolecule(boolean minimize){
        // Monomer one energy
        Minimize monomerMinEngine;
        for(Atom a: m2Atoms){ a.setUse(false); }
        if(minimize) {
            logger.info("\n --------- Minimize Monomer 1 --------- ");
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE);
            logger.info("\n --------- Monomer 1 Static Torsion Scan --------- ");
            TorsionSearch m1TorsionSearch = new TorsionSearch(mola, mola.getMoleculeArray()[0], 32, 1);
            m1TorsionSearch.staticAnalysis(0, 100);
            if (!m1TorsionSearch.getStates().isEmpty()) {
                AssemblyState minState = m1TorsionSearch.getStates().get(0);
                minState.revertState();
            }
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x);
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT);
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x);
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL);
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(eps, maxIter).getCoordinates(x);
        }
        logger.info("\n --------- Monomer 1 Energy Breakdown --------- ");
        double monomerEnergy = forceFieldEnergy.energy(x, true);
        for(Atom a: m2Atoms){ a.setUse(true); }

        // Monomer two energy
        for(Atom a: m1Atoms){ a.setUse(false); }
        if(minimize) {
            logger.info("\n --------- Minimize Monomer 2 --------- ");
            logger.info("\n --------- Monomer 2 Static Torsion Scan --------- ");
            TorsionSearch m2TorsionSearch = new TorsionSearch(mola, mola.getMoleculeArray()[1], 32, 1);
            m2TorsionSearch.staticAnalysis(0, 100);
            if (!m2TorsionSearch.getStates().isEmpty()) {
                AssemblyState minState = m2TorsionSearch.getStates().get(0);
                minState.revertState();
            }
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE);
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x);
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT);
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(1.0, maxIter).getCoordinates(x);
            forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL);
            monomerMinEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
            monomerMinEngine.minimize(eps, maxIter).getCoordinates(x);
        }
        logger.info("\n --------- Monomer 2 Energy Breakdown --------- ");
        double monomerEnergy2 = forceFieldEnergy.energy(x, true);
        for(Atom a: m1Atoms){ a.setUse(true); }

        // Log potentials
        logger.info(format("\n %-29s%12.7f kcal/mol", "Monomer energy 1:", monomerEnergy));
        logger.info(format(" %-29s%12.7f kcal/mol", "Monomer energy 2:", monomerEnergy2));

        m1MinEnergy = monomerEnergy;
        m2MinEnergy = monomerEnergy2;
        totalMonomerMinimizedEnergy = monomerEnergy + monomerEnergy2;
    }
    private void minimizeSystem(Atom a, Atom b) throws Exception{
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE);
        if(tScan){ // Molecules feel each other
            logger.info("\n --------- Monomer 1 Static Torsion Scan --------- ");
            forceFieldEnergy.getCoordinates(x);
            double tscanE = forceFieldEnergy.energy(x, false);
            TorsionSearch m1TorsionSearch = new TorsionSearch(mola, mola.getMoleculeArray()[0], 32, 1);
            m1TorsionSearch.staticAnalysis(0, 100);
            if(!m1TorsionSearch.getStates().isEmpty()) {
                AssemblyState minState = m1TorsionSearch.getStates().get(0);
                minState.revertState();
            }
            forceFieldEnergy.getCoordinates(x);
            double tscanEAfter = forceFieldEnergy.energy(x, false);
            logger.info("\n Energy before static torsion scan of monomer 1: " + tscanE);
            logger.info(" Energy after static torsion scan of monomer 1: " + tscanEAfter);

            logger.info("\n --------- Monomer 2 Static Torsion Scan --------- ");
            forceFieldEnergy.getCoordinates(x);
            tscanE = forceFieldEnergy.energy(x, false);
            TorsionSearch m2TorsionSearch = new TorsionSearch(mola, mola.getMoleculeArray()[1], 32, 1);
            m2TorsionSearch.staticAnalysis(0, 100);
            if(!m2TorsionSearch.getStates().isEmpty()) {
                AssemblyState minState = m2TorsionSearch.getStates().get(0);
                minState.revertState();
            }
            forceFieldEnergy.getCoordinates(x);
            tscanEAfter = forceFieldEnergy.energy(x, false);
            logger.info("\n Energy before static torsion scan of monomer 2: " + tscanE);
            logger.info(" Energy after static torsion scan of monomer 2: " + tscanEAfter);
        }
        forceFieldEnergy.getCoordinates(x);
        double e = forceFieldEnergy.energy(x, false);
        if (e > 100000){
            throw new Exception("Energy too high to minimize.");
        }
        // Set up restraintBond
        BondType restraint = new BondType(new int[]{a.getAtomicNumber(), b.getAtomicNumber()},
                1000.0,
                this.hBondDist,
                BondType.BondFunction.FLAT_BOTTOM_QUARTIC,
                this.flatBottomRadius);
        RestraintBond restraintBond = new RestraintBond(a, b,
                null,
                false,
                0.0, 0.0,
                null);
        restraintBond.setBondType(restraint);
        //Minimize NONE and DIRECT and MUTUAL
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE);
        Minimize minEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
        minEngine.minimize(1.0, this.maxIter);
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.DIRECT);
        minEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
        minEngine.minimize(1.0, this.maxIter);
        forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL);
        minEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
        minEngine.minimize(this.eps, this.maxIter);
        // Delete restraintBond
        a.getBonds().remove(restraintBond);
        b.getBonds().remove(restraintBond);
        a.update();
        b.update();
        mola.getBondList().remove(restraintBond);
        mola.update();
    }
    private void calculateMeanStd(){
        ArrayList<Double> energies = this.getEnergies();
        double highOutlierCutoff = Double.MAX_VALUE;
        double lowOutlierCutoff = Double.MIN_VALUE;
        // Calculate the interquartile range
        if(energies.size() > 4) {
            double iqr = energies.get(energies.size() / 4) - energies.get(3 * energies.size() / 4);
            double q3 = energies.get(3 * energies.size() / 4);
            double q1 = energies.get(energies.size() / 4);
            highOutlierCutoff = q3 + 1.5 * iqr;
            lowOutlierCutoff = q1 - 1.5 * iqr;
        }

        // Calculate the means
        int count = 0;
        double sum = 0.0;
        double sumNoOutlier = 0.0;
        for(double e: energies){
            if(e < highOutlierCutoff && e > lowOutlierCutoff){
                sumNoOutlier += e;
                count++;
            }
            sum += e;
        }
        averageEnergy = sum / energies.size();
        averageEnergyNoOutlier = sumNoOutlier / count;

        // Calculate stds
        double sumOfSquares = 0.0;
        double sumOfSquaresNoOutlier = 0.0;
        for(double e: energies){
            sumOfSquares += (e - averageEnergy) * (e - averageEnergy);
            if(e < highOutlierCutoff && e > lowOutlierCutoff){
                sumOfSquaresNoOutlier += (e - averageEnergyNoOutlier) * (e - averageEnergyNoOutlier);
            }
        }
        stdOfEnergies = Math.sqrt(sumOfSquares / energies.size());
        stdOfEnergiesNoOutlier = Math.sqrt(sumOfSquaresNoOutlier / count);
    }

    private static void alignMoleculeCOMtoAtomVecWithAxis(Molecule m, Atom a, double[] axis, Atom[] mAtoms){
        // Get center of mass of moleculeOneAtoms
        double[] moleculeOneCOM = getCOM(m);
        for(int i = 0; i < mAtoms.length; i++){
            mAtoms[i].move(new double[] {-moleculeOneCOM[0], -moleculeOneCOM[1], -moleculeOneCOM[2]});
        }
        // Get coordinates of a
        double[] aCoords = a.getXYZ().copy().get();
        // Get rotation matrix to align dipole moment with z-axis
        double[][] rotation = getRotationBetween(aCoords, axis);
        // Get a reference to the moleculeOneAtoms
        double[] moleculeAtomicPositions = new double[mAtoms.length * 3];
        for(int i = 0; i < mAtoms.length; i++){
            moleculeAtomicPositions[i*3] = mAtoms[i].getX();
            moleculeAtomicPositions[i*3 + 1] = mAtoms[i].getY();
            moleculeAtomicPositions[i*3 + 2] = mAtoms[i].getZ();
        }
        applyRotation(moleculeAtomicPositions, rotation);
        // Move atoms into rotated positions
        for(int i = 0; i < mAtoms.length; i++){
            mAtoms[i].setXYZ(new double[] {moleculeAtomicPositions[i*3],
                    moleculeAtomicPositions[i*3 + 1],
                    moleculeAtomicPositions[i*3 + 2]});
        }
    }

    /**
     * Gets the center of mass of a molecule
     * @param m
     * @return x,y,z coordinates of center of mass
     */
    private static double[] getCOM(Molecule m){
        // Get center of mass of moleculeOneAtoms
        double[] moleculeOneCOM = new double[3];
        double totalMass = 0.0;
        for(Atom s: m.getAtomList()){
            double[] pos = s.getXYZ().get();
            moleculeOneCOM[0] += pos[0] * s.getMass();
            moleculeOneCOM[1] += pos[1] * s.getMass();
            moleculeOneCOM[2] += pos[2] * s.getMass();
            totalMass += s.getMass();
        }
        totalMass = 1 / totalMass;
        moleculeOneCOM[0] *= totalMass;
        moleculeOneCOM[1] *= totalMass;
        moleculeOneCOM[2] *= totalMass;

        return moleculeOneCOM;
    }

    /**
     * Gets the rotation matrix between two vectors
     * @param v1
     * @param v2
     * @return rotation matrix
     */
    private static double[][] getRotationBetween(double[] v1, double[] v2){
        // Normalize v1 and v2
        double v1Norm = 1/ Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        double v2Norm = 1 / Math.sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
        for(int i = 0; i < 3; i++){
            v1[i] *= v1Norm;
            v2[i] *= v2Norm;
        }
        // Cross product between targetVector and z-axis
        double[] crossProduct = new double[3];
        crossProduct[0] = v1[1] * v2[2] - v1[2] * v2[1];
        crossProduct[1] = v1[2] * v2[0] - v1[0] * v2[2];
        crossProduct[2] = v1[0] * v2[1] - v1[1] * v2[0];
        // Normalize cross product
        double crossProductNorm = 1 / Math.sqrt(crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1] + crossProduct[2] * crossProduct[2]);
        for(int i = 0; i < 3; i++){
            crossProduct[i] *= crossProductNorm;
        }
        // Dot product between v1 and z-axis
        double dotProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        // Angle between v1 and z-axis
        double theta = Math.acos(dotProduct);
        double[] u = crossProduct;
        // Define quaternion from axis-angle
        double[] quaternion = new double[4];
        quaternion[0] = cos(theta/2);
        quaternion[1] = u[0] * sin(theta/2);
        quaternion[2] = u[1] * sin(theta/2);
        quaternion[3] = u[2] * sin(theta/2);
        // Normalize quaternion
        double quaternionNorm = 1 / Math.sqrt(quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]
                + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);
        for(int i = 0; i < 4; i++){ quaternion[i] *= quaternionNorm; }
        // Useful storage
        double q1q1 = quaternion[1] * quaternion[1];
        double q2q2 = quaternion[2] * quaternion[2];
        double q3q3 = quaternion[3] * quaternion[3];
        double q0q1 = quaternion[0] * quaternion[1];
        double q0q2 = quaternion[0] * quaternion[2];
        double q0q3 = quaternion[0] * quaternion[3];
        double q1q2 = quaternion[1] * quaternion[2];
        double q1q3 = quaternion[1] * quaternion[3];
        double q2q3 = quaternion[2] * quaternion[3];
        // Quaternion rotation matrix
        double[][] rotation = new double[3][3];
        rotation[0][0] = 1 - 2 * (q2q2 + q3q3);
        rotation[0][1] = 2 * (q1q2 - q0q3);
        rotation[0][2] = 2 * (q1q3 + q0q2);
        rotation[1][0] = 2 * (q1q2 + q0q3);
        rotation[1][1] = 1 - 2 * (q1q1 + q3q3);
        rotation[1][2] = 2 * (q2q3 - q0q1);
        rotation[2][0] = 2 * (q1q3 - q0q2);
        rotation[2][1] = 2 * (q2q3 + q0q1);
        rotation[2][2] = 1 - 2 * (q1q1 + q2q2);

        return rotation;
    }


    /**
     * Implements StateContainer to store the coordinates of a state and its energy
     */
    private record StateContainer(AssemblyState state, double e) implements Comparable<StateContainer>
    {
        AssemblyState getState() {
            return state;
        }

        double getEnergy() {
            return e;
        }

        @Override
        public int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy());
        }
    }
}
