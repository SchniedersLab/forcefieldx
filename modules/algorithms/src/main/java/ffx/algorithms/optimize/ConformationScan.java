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
 * This class is for a configuration optimization of two small systems. The search
 * consists of aligning heavy (except carbon) and hydrogen atoms at a h-bond distance
 * with their center of masses behind them all on the z-axis in hopes that a global
 * conformational minima will be found by allowing for strong h-bonds to occur. After
 * aligning, an optional static torsion scan and NONE,DIRECT,MUTUAL minimizations are
 * performed. A binding energy can be calculated with this by:
 *
 * ddEbind = dEbind (of two systems)
 *         - dEbind (of one system w/ itself)
 *         - dEbind (of the other system w/ itself)
 */
public class ConformationScan {

    private static final Logger logger = Logger.getLogger(ConformationScan.class.getName());

    private final MolecularAssembly mola;
    private final ForceFieldEnergy forceFieldEnergy;
    private final AssemblyState initState;
    private double[] x;
    private AlgorithmListener algorithmListener;

    private final Molecule[] s1;
    private final Molecule[] s2;
    private final Atom[] s1Atoms;
    private final Atom[] s2Atoms;
    private final ArrayList<Atom> s1TargetAtoms;
    private final ArrayList<Atom> s2TargetAtoms;

    private AbstractQueue<StateContainer> statesQueue;

    private final double eps;
    private final int maxIter;
    private double hBondDist;
    private double flatBottomRadius;
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
            Molecule[] systemOne,
            Molecule[] systemTwo,
            double minimizeEps,
            int minimizeMaxIterations,
            boolean staticTorsionScan,
            boolean excludeExtraHydrogen,
            boolean minimize
    )
    {
        mola = molecularAssembly;
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        initState = new AssemblyState(mola);
        s1 = systemOne;
        s2 = systemTwo;
        eps = minimizeEps;
        maxIter = minimizeMaxIterations;
        tScan = staticTorsionScan;
        excludeH = excludeExtraHydrogen;
        this.minimize = minimize;
        s1Atoms = getAtomsFromMoleculeArray(s1);
        s2Atoms = getAtomsFromMoleculeArray(s2);
        s1TargetAtoms = new ArrayList<>();
        s2TargetAtoms = new ArrayList<>();
        setTargetAtoms(mola.getAtomArray()); // TODO: Adjust this to handle systems not just molecules
        statesQueue = MinMaxPriorityQueue.maximumSize(s1TargetAtoms.size() * s2TargetAtoms.size())
                .create();
        if(!minimize) { // If minimization is not performed, then the queue will be in order of the loop
            statesQueue = new ArrayBlockingQueue<>(s1TargetAtoms.size() * s2TargetAtoms.size() + 1);
        }
    }

    public void scan(){
        systemEnergies(); // get energy of each system by itself for binding e calculation
        double[] zAxis = new double[]{0,0,1};
        // Loop through interactions between the two molecules --> Not necessarily symmetric but close?
        int loopCounter = 1;
        for(Atom a: s1TargetAtoms){
            for(Atom b: s2TargetAtoms){
                zAxis[2] = 1;
                initState.revertState(); // all trials need to be initialized from the initial states
                alignSystemCOMtoAtomVecWithAxis(a, zAxis, s1Atoms);
                zAxis[2] = -1;
                logger.info("\n ----- Trial " + loopCounter + " out of " +
                        (s2TargetAtoms.size() * s1TargetAtoms.size()) + " -----");
                loopCounter++;
                alignSystemCOMtoAtomVecWithAxis(b, zAxis, s2Atoms);
                // Move the second molecule to the perfect h-bond distance
                hBondDist = 2.0;
                double[] hBondVector = new double[]{0, 0, a.getZ() - b.getZ() + hBondDist};
                logger.info(" Initial H-bond distance: " + hBondVector[2]);
                // Finds and moves to minimial vector (no polarization to avoid failures)
                forceFieldEnergy.getPmeNode().setPolarization(Polarization.NONE);
                hBondVector[2] += minimizeVector(hBondVector, -15, 15)[2];
                forceFieldEnergy.getPmeNode().setPolarization(Polarization.MUTUAL);
                logger.info(" Best H-bond distance: " + hBondVector[2]);
                hBondDist = hBondVector[2] - (a.getZ() - b.getZ());
                flatBottomRadius = Math.abs(hBondDist / 2.0);
                logger.info(" Flat bottom radius: " + flatBottomRadius);
                // Minimize the energy of the system subject to a harmonic restraint on the distance
                // between the two atoms. Keep the state if minimization works.
                try {
                    mola.update();
                    forceFieldEnergy.getCoordinates(x);
                    int status = 0;
                    if(minimize) {
                        status = minimizeSystem(a, b);
                    }
                    if(status != -1) {
                        forceFieldEnergy.getCoordinates(x);
                        double e = forceFieldEnergy.energy(x, true) - totalMonomerMinimizedEnergy;
                        logger.info("\n Binding energy of trial " + (loopCounter - 1) + ": " + e);
                        statesQueue.add(new StateContainer(new AssemblyState(mola), e));
                    } else {
                        logger.warning(" Minimization failed. No state will be saved.");
                        //statesQueue.add(new StateContainer(new AssemblyState(mola), -1)); Use for debugging
                    }
                } catch (Exception ignored) {
                    logger.warning(" Minimization failed. No state will be saved.");
                    //statesQueue.add(new StateContainer(new AssemblyState(mola), -1)); Use for debugging
                    //e.printStackTrace()
                }
            }
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
    public boolean writeStructuresToXYZ(File outputFile){
        XYZFilter xyzFilter = new XYZFilter(outputFile, mola, mola.getForceField(), mola.getForceField().getProperties());
        logger.info("\n Writing structures to " + outputFile.getAbsolutePath());
        int count = 1;
        StateContainer[] states = new StateContainer[statesQueue.size()];
        if(minimize){
            // Pull states out in order
            int index = 0;
            while(!statesQueue.isEmpty() && index < states.length){
                states[index] = statesQueue.poll();
                index++;
            }
            if(index < states.length){
                logger.warning(" Not all states will be saved. ");
                return false;
            }
        }
        else{
            statesQueue.toArray(states);
        }
        if(states.length == 0){
            logger.warning(" No states were saved. No structures will be written.");
            return false;
        }
        for(StateContainer s: states){
            AssemblyState state = s.getState();
            double energy = s.getEnergy();
            logger.info(" Writing configuration #" + count + " with energy: " + energy);
            state.revertState();
            xyzFilter.writeFile(outputFile, count != 1);
            count++;
        }
        return true;
    }

    public static ArrayList<Double> logBindingEnergyCalculation(ConformationScan m1, ConformationScan m2, ConformationScan dimer){
        ArrayList<Double> bindingEnergies = new ArrayList<>();
        logger.info("\n ------------------------- Binding Energy -------------------------");
        logger.info(" ddE (of binding) = dE (of coformer dimer binding) - dE (average of both monomer dimer bindings summed)");
        logger.info("\n Calculation using minimized energies: ");
        double averageMonomerEnergy = (m1.getMinimumEnergy() + m2.getMinimumEnergy()) / 2;
        double bindingE = dimer.getMinimumEnergy() - averageMonomerEnergy;
        logger.info(format(" Minimal Structure    ddE = %8.3f - %8.3f = %12.5f",
                dimer.getMinimumEnergy(),
                averageMonomerEnergy,
                bindingE));
        bindingEnergies.add(bindingE);
        logger.info("\n Calculation using average energies: ");
        averageMonomerEnergy = (m1.getAverageEnergy() + m2.getAverageEnergy()) / 2;
        bindingE = dimer.getAverageEnergy() - averageMonomerEnergy;
        logger.info(format(" Average              ddE = %8.3f - %8.3f = %12.5f",
                dimer.getAverageEnergy(),
                averageMonomerEnergy,
                bindingE));
        bindingEnergies.add(bindingE);
        logger.info("\n Calculation using average energies (no outliers): ");
        averageMonomerEnergy = (m1.getAverageEnergyNoOutlier() + m2.getAverageEnergyNoOutlier()) / 2;
        bindingE = dimer.getAverageEnergyNoOutlier() - averageMonomerEnergy;
        logger.info(format(" Average (no outlier) ddE = %8.3f - %8.3f = %12.5f",
                dimer.getAverageEnergyNoOutlier(),
                averageMonomerEnergy,
                bindingE));
        bindingEnergies.add(bindingE);
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
    private double[] minimizeVector(double[] hBondVector, int lowBound, int highBound) {
        highBound += 1; // To include the highBound
        // Grid search from lowBound to highBound w/ 1 ang steps
        double[] coarsePotentialSurface = new double[(int) Math.abs(highBound - lowBound)];
        double[] zSearched = new double[(int) Math.abs(highBound - lowBound)]; // Relative to hBondVector[2]
        double[] coarseVector = new double[3];
        int minIndex = -1;
        double minE = Double.MAX_VALUE;
        coarseVector[2] = hBondVector[2] + lowBound; // Bounds depend on the hBondVector[2] value
        for(int i = 0; i < Math.abs(highBound - lowBound); i++){
            zSearched[i] = i == 0 ? lowBound : zSearched[i - 1] + 1;
            // Move mol2 to new position
            for(Atom a: s2Atoms){
                a.move(coarseVector);
            }
            // Calculate the energy of the system
            forceFieldEnergy.getCoordinates(x);
            coarsePotentialSurface[i] = forceFieldEnergy.energy(x, false);
            if(coarsePotentialSurface[i] < minE){
                minE = coarsePotentialSurface[i];
                minIndex = i;
            }
            coarseVector[2] = 1;
        }
        // Return to hBond position using difference between hBondVector[2] and zSearched[zSearched.length - 1] to start
        // relative search
        for (Atom s2Atom : s2Atoms) {
            s2Atom.move(new double[]{0, 0, hBondVector[2] - zSearched[zSearched.length - 1]});
        }
        //logger.info("Z space: " + Arrays.toString(zSearched));
        //logger.info("Coarse potential surface: " + Arrays.toString(coarsePotentialSurface));

        // Refine the minimum with binary search
        double[] refinedVector = new double[3]; // new c position relative to previous c
        double a;
        double aPotential;
        if (minIndex == 0) {
            a = zSearched[minIndex + 1];
            aPotential = coarsePotentialSurface[minIndex + 1];
        } else if (minIndex == coarsePotentialSurface.length - 1) {
            a = zSearched[minIndex - 1];
            aPotential = coarsePotentialSurface[minIndex - 1];
        } else {
            a = coarsePotentialSurface[minIndex - 1] < coarsePotentialSurface[minIndex + 1] ? // relative to hbond
                    zSearched[minIndex - 1] : zSearched[minIndex + 1];
            aPotential = Math.min(coarsePotentialSurface[minIndex - 1], coarsePotentialSurface[minIndex + 1]);
        }
        double b = zSearched[minIndex]; // relative to hbond
        double bPotential = coarsePotentialSurface[minIndex];
        double c = 0; // Store position of previous c
        double convergence = 1e-5;
        while(Math.abs(aPotential - bPotential) > convergence){
            refinedVector[2] = (a + b) / 2 - c;
            // Move mol2 to new position between a and b
            for(Atom a1: s2Atoms){
                a1.move(refinedVector);
            }
            // Calculate the energy of the system
            forceFieldEnergy.getCoordinates(x);
            double e = forceFieldEnergy.energy(x, false);
            // Set up next step
            minE = Math.min(e, minE);
            if(aPotential > bPotential && e < aPotential){
                a = (a+b)/2;
                aPotential = e;
                c = a;
            } else if(e < bPotential){
                b = (a+b)/2;
                bPotential = e;
                c = b;
            } else{
                // Move to a or b
                c = (a + b) / 2;
                if(aPotential < bPotential){
                    refinedVector[2] = a - c;
                    for(Atom a1: s2Atoms){
                        a1.move(refinedVector);
                    }
                } else{
                    refinedVector[2] = b - c;
                    for(Atom a1: s2Atoms){
                        a1.move(refinedVector);
                    }
                }
                break; // e > bPotential and e > aPotential --> numerical stability err or P.E.S. is not parabolic
            }
        }
        refinedVector[2] = aPotential < bPotential ? a : b;
        //logger.info("Hbond vector: " + Arrays.toString(hBondVector));
        //logger.info("Refined vector (relative to hbond): " + Arrays.toString(refinedVector));
        return refinedVector; // relative to hbond
    }

    // This used to be for minimization but caused issues
    private void systemEnergies(){
        // Monomer one energy
        for(Atom a: s2Atoms){ a.setUse(false); }
        logger.info("\n --------- System 1 Energy Breakdown --------- ");
        double monomerEnergy = forceFieldEnergy.energy(x, true);
        for(Atom a: s2Atoms){ a.setUse(true); }

        // Monomer two energy
        for(Atom a: s1Atoms){ a.setUse(false); }
        logger.info("\n --------- System 2 Energy Breakdown --------- ");
        double monomerEnergy2 = forceFieldEnergy.energy(x, true);
        for(Atom a: s1Atoms){ a.setUse(true); }

        // Log potentials
        logger.info(format("\n %-29s%12.7f kcal/mol", "System energy 1:", monomerEnergy));
        logger.info(format(" %-29s%12.7f kcal/mol", "System energy 2:", monomerEnergy2));

        m1MinEnergy = monomerEnergy;
        m2MinEnergy = monomerEnergy2;
        totalMonomerMinimizedEnergy = monomerEnergy + monomerEnergy2;
    }

    /**
     * This function minimizes a single molecule using static tscan and a minimization engine
     * @param m
     * @return
     */
    private int minimizeMolecule(Molecule m){
        if (this.tScan){
            staticScanMolecule(m);
        }
        Minimize minimize = new Minimize(this.mola, this.forceFieldEnergy, this.algorithmListener);
        minimize.minimize(this.eps, this.maxIter).getCoordinates(this.x);
        return minimize.getStatus();
    }

    private void staticScanMolecule(Molecule m){
        TorsionSearch torsionSearch = new TorsionSearch(this.mola, m, 32, 1);
        torsionSearch.staticAnalysis(0, 100);
        if(!torsionSearch.getStates().isEmpty()){
            AssemblyState minState = torsionSearch.getStates().get(0);
            minState.revertState();
        }
    }
    private int minimizeSystem(Atom a, Atom b) throws Exception{
        if(tScan){
            // Molecules feel each other
            logger.info("\n --------- System 1 Static Torsion Scan --------- ");
            forceFieldEnergy.getCoordinates(x);
            double tscanE = forceFieldEnergy.energy(x, false);
            for(Molecule m: s1) {
                staticScanMolecule(m);
            }
            forceFieldEnergy.getCoordinates(x);
            double tscanEAfter = forceFieldEnergy.energy(x, false);
            logger.info("\n Energy before static torsion scan of system 1: " + tscanE);
            logger.info(" Energy after static torsion scan of system 1: " + tscanEAfter);

            logger.info("\n --------- System 2 Static Torsion Scan --------- ");
            forceFieldEnergy.getCoordinates(x);
            tscanE = forceFieldEnergy.energy(x, false);
            for(Molecule m: s2) {
                staticScanMolecule(m);
            }
            forceFieldEnergy.getCoordinates(x);
            tscanEAfter = forceFieldEnergy.energy(x, false);
            logger.info("\n Energy before static torsion scan of system 2: " + tscanE);
            logger.info(" Energy after static torsion scan of system 2: " + tscanEAfter);
        }
        forceFieldEnergy.getCoordinates(x);
        double e = forceFieldEnergy.energy(x, true);
        if (e > 1000000){
            throw new Exception(" Energy too high to minimize.");
        }
        // Set up restraintBond
        BondType restraint = new BondType(new int[]{a.getAtomicNumber(), b.getAtomicNumber()},
                100.0,
                this.hBondDist,
                BondType.BondFunction.FLAT_BOTTOM_QUARTIC,
                this.flatBottomRadius);
        RestraintBond restraintBond = new RestraintBond(a, b,
                null,
                false,
                0.0, 0.0,
                null);
        restraintBond.setBondType(restraint);
        Minimize minEngine = new Minimize(mola, forceFieldEnergy, algorithmListener);
        try {
            minEngine.minimize(this.eps, this.maxIter);
        } catch (Exception ex){
            // Delete restraintBond no matter what
            a.getBonds().remove(restraintBond);
            b.getBonds().remove(restraintBond);
            a.update();
            b.update();
            mola.getBondList().remove(restraintBond);
            mola.update();
            return -1;
        }
        // Delete restraintBond
        a.getBonds().remove(restraintBond);
        b.getBonds().remove(restraintBond);
        a.update();
        b.update();
        mola.getBondList().remove(restraintBond);
        mola.update();
        return minEngine.getStatus();
    }
    private void setTargetAtoms(Atom[] atoms){
        for(Atom a: atoms){
            if(a.getAtomType().atomicNumber == 7|| a.getAtomType().atomicNumber == 8
                    || a.getAtomType().atomicNumber == 9 || a.getAtomType().atomicNumber == 15
                    || a.getAtomType().atomicNumber == 16 || a.getAtomType().atomicNumber == 17){ // N,O,F,P,S,Cl
                if(a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                    s1TargetAtoms.add(a);
                } else{
                    s2TargetAtoms.add(a);
                }

                // Searching for bonded h's only if we are excluding H's that aren't bonded to electronegative atoms
                if(excludeH) {
                    for (Bond b : a.getBonds()) {
                        int num = b.get1_2(a).getAtomType().atomicNumber;
                        if (num == 1) {
                            if (a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                                s1TargetAtoms.add(a);
                            } else {
                                s2TargetAtoms.add(a);
                            }
                        }
                    }
                }
            }
            else if(a.getAtomType().atomicNumber == 1 && !excludeH){ // Add all H's
                if(a.getMoleculeNumber() == mola.getMoleculeNumbers()[0]) {
                    s1TargetAtoms.add(a);
                } else{
                    s2TargetAtoms.add(a);
                }
            }
        }

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

    private static void alignSystemCOMtoAtomVecWithAxis(Atom a, double[] axis, Atom[] mAtoms){
        // Get center of mass of moleculeOneAtoms
        double[] moleculeOneCOM = getCOM(mAtoms);
        for (Atom mAtom : mAtoms) {
            mAtom.move(new double[]{-moleculeOneCOM[0], -moleculeOneCOM[1], -moleculeOneCOM[2]});
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
     * Gets the center of mass of a set of atoms
     * @param atoms
     * @return x,y,z coordinates of center of mass
     */
    private static double[] getCOM(Atom[] atoms){
        // Get center of mass of moleculeOneAtoms
        double[] COM = new double[3];
        double totalMass = 0.0;
        for(Atom s: atoms){
            double[] pos = s.getXYZ().get();
            COM[0] += pos[0] * s.getMass();
            COM[1] += pos[1] * s.getMass();
            COM[2] += pos[2] * s.getMass();
            totalMass += s.getMass();
        }
        totalMass = 1 / totalMass;
        COM[0] *= totalMass;
        COM[1] *= totalMass;
        COM[2] *= totalMass;

        return COM;
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
        // Define quaternion from axis-angle
        double[] quaternion = new double[4];
        quaternion[0] = cos(theta/2);
        quaternion[1] = crossProduct[0] * sin(theta/2);
        quaternion[2] = crossProduct[1] * sin(theta/2);
        quaternion[3] = crossProduct[2] * sin(theta/2);
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

    private static Atom[] getAtomsFromMoleculeArray(Molecule[] system){
        ArrayList<Atom> atoms = new ArrayList<>();
        for(Molecule m: system){
            atoms.addAll(m.getAtomList());
        }
        return atoms.toArray(new Atom[0]);
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
