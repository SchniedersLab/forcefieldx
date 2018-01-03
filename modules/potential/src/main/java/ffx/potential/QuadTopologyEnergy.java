/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.utils.EnergyException;

/**
 * Implements an error-canceling quad topology, where two large dual-topology 
 * simulation legs are run simultaneously to arrive at a small sum. This 
 * implementation permits sharing of coordinates between the dual topology; 
 * results in energy function of E(A1, A2, B1, B2) = (1-l)A1 + l*A2 + l*B1 +
 * (1-l)B2. When coordinates are shared, this can entail atoms feeling 
 * approximately twice the force as an ordinary atom, possibly requiring a
 * reduced inner timestep.
 * 
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class QuadTopologyEnergy implements CrystalPotential, LambdaInterface {
    private final static Logger logger = Logger.getLogger(QuadTopologyEnergy.class.getName());
    private final DualTopologyEnergy dualTopA;
    private final DualTopologyEnergy dualTopB;
    private final LambdaInterface linterA;
    private final LambdaInterface linterB;
    
    private final int nVarA;
    private final int nVarB;
    private final int nShared;
    private final int uniqueA;
    private final int uniqueB;
    private final int nVarTot;
    
    /**
     * Following arrays keep track of which index in the child DualTopology is 
     * linked to which index in the QuadTopology variable array.
     */
    private final int[] indexAToGlobal;
    private final int[] indexBToGlobal;
    private final int[] indexGlobalToA;
    private final int[] indexGlobalToB;
    
    private final double[] mass;
    private final double[] xA;
    private final double[] xB;
    private final double[] gA;
    private final double[] gB;
    /**
     * tempA and B arrays are used to hold the result of get methods applied
     * on dual topologies A and B; this saves the time of initializing new arrays
     * to hold the data for a very short time. Not thread-safe; if multiple
     * threads are to access this data, will need to refactor the class to 
     * synchronize on these arrays so only one method can use them at a time.
     */
    private final double[] tempA;
    private final double[] tempB;
    
    private STATE state = STATE.BOTH;
    private double lambda;
    private double totalEnergy;
    private double energyA;
    private double energyB;
    private double dEdL, dEdL_A, dEdL_B;
    private double d2EdL2, d2EdL2_A, d2EdL2_B;
    
    private boolean inParallel = false;
    private ParallelTeam team;
    private final EnergyRegion region;
    
    /**
     * Scaling and de-scaling will be applied inside DualTopologyEnergy.
     */
    private double[] scaling;
    private VARIABLE_TYPE[] types = null;
    
     /**
     * General structure: first layer will be the "A/B" layer, consisting of the
     * two dual topologies. The second layer will be the "1/2" layer, consisting
     * of the assemblies/ForceFieldEnergies associated with each dual topology.
     * Arrays will be addressed as [A/B=0/1][1/2=0/1].
     * 
     * For example, if running a dual-force-field calculation, A might be the
     * original molecule, going from AMOEBA to a fixed-charge force field, while
     * B is a modified molecule going from fixed-charge to AMOEBA. Within those,
     * A1 would be AMOEBA-original, A2 would be AMBER-original, B1 would be 
     * AMBER-modified, B2 would be AMOEBA-modified.
     * 
     * This constructor assumes there are no unique atoms in either topology, 
     * and pins everything.
     * 
     * @param dualTopologyA A DualTopologyEnergy
     * @param dualTopologyB A DualTopologyEnergy
     */      
    public QuadTopologyEnergy(DualTopologyEnergy dualTopologyA, DualTopologyEnergy dualTopologyB) {
        this(dualTopologyA, dualTopologyB, null, null);
    }
    
    /**
     * General structure: first layer will be the "A/B" layer, consisting of the
     * two dual topologies. The second layer will be the "1/2" layer, consisting
     * of the assemblies/ForceFieldEnergies associated with each dual topology.
     * Arrays will be addressed as [A/B=0/1][1/2=0/1].
     * 
     * For example, if running a dual-force-field calculation, A might be the
     * original molecule, going from AMOEBA to a fixed-charge force field, while
     * B is a modified molecule going from fixed-charge to AMOEBA. Within those,
     * A1 would be AMOEBA-original, A2 would be AMBER-original, B1 would be 
     * AMBER-modified, B2 would be AMOEBA-modified.
     * 
     * @param dualTopologyA A DualTopologyEnergy
     * @param dualTopologyB A DualTopologyEnergy
     * @param uniqueAList Variables unique to A
     * @param uniqueBList Variables unique to B
     */        
    public QuadTopologyEnergy(DualTopologyEnergy dualTopologyA, DualTopologyEnergy dualTopologyB, List<Integer> uniqueAList, List<Integer> uniqueBList) {
        this.dualTopA = dualTopologyA;
        this.dualTopB = dualTopologyB;
        dualTopB.setCrystal(dualTopA.getCrystal());
        dualTopB.reloadCommonMasses(true);
        linterA = (LambdaInterface) dualTopologyA;
        linterB = (LambdaInterface) dualTopologyB;
        nVarA = dualTopA.getNumberOfVariables();
        nVarB = dualTopB.getNumberOfVariables();
        
        /**
         * Following logic is to A, deal with Integer to int array
         * conversion issues, and B, ensure the set of unique variables is a 
         * unique set. Unique elements may come from either the provided list, or
         * from the non-shared elements of the dual topologies.
         */
        Set<Integer> uniqueASet = new LinkedHashSet<>();
        if (uniqueAList != null) {
            uniqueASet.addAll(uniqueAList);
        }
        int nCommon = dualTopA.getNumSharedVariables();
        for (int i = nCommon; i < nVarA; i++) {
            uniqueASet.add(i);
        }

        uniqueAList = new ArrayList<>(uniqueASet);
        uniqueA = uniqueAList.size();
        int[] uniquesA = new int[uniqueA];
        for (int i = 0; i < uniqueA; i++) {
            uniquesA[i] = uniqueAList.get(i);
        }
        // Following can replace above 5 lines without constructing an intermediate List.
        //uniquesA = uniqueASet.stream().mapToInt(Integer::intValue).toArray();
        //uniqueA = uniquesA.length;
        
        /**
         * Following logic is to A, deal with Integer to int array
         * conversion issues, and B, ensure the set of unique variables is a 
         * unique set. Unique elements may come from either the provided list, or
         * from the non-shared elements of the dual topologies.
         */
        Set<Integer> uniqueBSet = new LinkedHashSet<>();
        if (uniqueBList != null) {
            uniqueBSet.addAll(uniqueBList);
        }
        nCommon = dualTopB.getNumSharedVariables();
        for (int i = nCommon; i < nVarB; i++) {
            uniqueBSet.add(i);
        }
        uniqueBList = new ArrayList<>(uniqueBSet);
        uniqueB = uniqueBList.size();
        int[] uniquesB = new int[uniqueB];
        for (int i = 0; i < uniqueB; i++) {
            uniquesB[i] = uniqueBList.get(i);
        }
        
        nShared = nVarA - uniqueA;
        assert (nShared == nVarB - uniqueB);
        nVarTot = nShared + uniqueA + uniqueB;
        
        indexAToGlobal = new int[nVarA];
        indexBToGlobal = new int[nVarB];
        indexGlobalToA = new int[nVarTot];
        indexGlobalToB = new int[nVarTot];
        // -1 indicates this index in the global array does not point to one in 
        // the dual topology.
        Arrays.fill(indexGlobalToA, -1);
        Arrays.fill(indexGlobalToB, -1);
        
        if (uniqueA > 0) {
            int commonIndex = 0;
            int uniqueIndex = 0;
            for (int i = 0; i < nVarA; i++) {
                if (uniqueIndex < uniqueA && i == uniquesA[uniqueIndex]) {
                    int destIndex = nShared + uniqueIndex;
                    indexAToGlobal[i] = destIndex;
                    indexGlobalToA[destIndex] = i;
                    ++uniqueIndex;
                } else {
                    indexAToGlobal[i] = commonIndex;
                    indexGlobalToA[commonIndex++] = i;
                }
            }
        } else {
            for (int i = 0; i < nVarA; i++) {
                indexAToGlobal[i] = i;
                indexGlobalToA[i] = i;
            }
        }
        
        if (uniqueB > 0) {
            int commonIndex = 0;
            int uniqueIndex = 0;
            for (int i = 0; i < nVarB; i++) {
                if (uniqueIndex < uniqueB && i == uniquesB[uniqueIndex]) {
                    int destIndex = nVarA + uniqueIndex;
                    indexBToGlobal[i] = destIndex;
                    indexGlobalToB[destIndex] = i;
                    ++uniqueIndex;
                } else {
                    indexBToGlobal[i] = commonIndex;
                    indexGlobalToB[commonIndex++] = i;
                }
            }
        } else {
            for (int i = 0; i < nVarB; i++) {
                indexBToGlobal[i] = i;
                indexGlobalToB[i] = i;
            }
        }
        
        xA = new double[nVarA];
        xB = new double[nVarB];
        gA = new double[nVarA];
        gB = new double[nVarB];
        tempA = new double[nVarA];
        tempB = new double[nVarB];
        mass = new double[nVarTot];
        doublesFrom(mass, dualTopA.getMass(), dualTopB.getMass());
        
        region = new EnergyRegion();
        team = new ParallelTeam(1);
    }
    
    /**
     * Copies from an object array of length nVarTot to two object arrays of 
     * length nVarA and nVarB.
     * 
     * @param <T> Type of object
     * @param from Copy from
     * @param toA Copy shared and A-specific to
     * @param toB Copy shared and B-specific to
     */
    private <T> void copyTo(T[] from, T[] toA, T[] toB) {
        if (toA == null) {
            toA = Arrays.copyOf(from, nVarA);
        }
        if (toB == null) {
            toB = Arrays.copyOf(from, nVarB);
        }
        for (int i = 0; i < nVarTot; i++) {
            int index = indexGlobalToA[i];
            if (index >= 0) {
                toA[index] = from[i];
            }
            index = indexGlobalToB[i];
            if (index >= 0) {
                toB[index] = from[i];
            }
        }
    }

    /**
     * Copies from object arrays of length nVarA and nVarB to an object array of
     * length nVarTot; asserts objects in common indices are equal. Should not
     * be used when the result of the common indices should be f(A,B)
     * 
     * @param <T> Type of object
     * @param to Copy to
     * @param fromA Copy shared and A-specific from
     * @param fromB Copy B-specific from
     */
    private <T> void copyFrom(T[] to, T[] fromA, T[] fromB) {
        if (to == null) {
            to = Arrays.copyOf(fromA, nVarTot);
        }
        for (int i = 0; i < nVarA; i++) {
            int index = indexAToGlobal[i];
            to[index] = fromA[i];
        }
        for (int i = 0; i < nVarB; i++) {
            int index = indexBToGlobal[i];
            // Assert either a unique variable, or equals what was added from A.
            assert (index >= nShared || to[index].equals(fromB[i]));
            to[index] = fromB[i];
        }
    }
    
    /**
     * Copies from an double array of length nVarTot to two double arrays of 
     * length nVarA and nVarB.
     *
     * @param from Copy from
     * @param toA Copy shared and A-specific to
     * @param toB Copy shared and B-specific to
     */
    private void doublesTo(double[] from, double[] toA, double[] toB) {
        toA = (toA == null) ? new double[nVarA] : toA;
        toB = (toB == null) ? new double[nVarB] : toB;
        for (int i = 0; i < nVarTot; i++) {
            int index = indexGlobalToA[i];
            if (index >= 0) {
                toA[index] = from[i];
            }
            index = indexGlobalToB[i];
            if (index >= 0) {
                toB[index] = from[i];
            }
        }
    }

    /**
     * Copies from double arrays of length nVarA and nVarB to an object array of
     * length nVarTot; asserts common indices are equal. Should not be used when 
     * the result of the common indices should be f(A,B)
     * 
     * @param to Copy to
     * @param fromA Copy shared and A-specific from
     * @param fromB Copy B-specific from
     */
    private void doublesFrom(double[] to, double[] fromA, double[] fromB) {
        to = (to == null) ? new double[nVarTot] : to;
        for (int i = 0; i < nVarA; i++) {
            to[indexAToGlobal[i]] = fromA[i];
        }
        for (int i = 0; i < nVarB; i++) {
            int index = indexBToGlobal[i];
            // Assert this is either a unique from B or it's equal to what came from A.
            assert (index >= nShared || to[index] == fromB[i]);
            to[index] = fromB[i];
        }
    }

    /**
     * Assigns common indices of to to be sum of fromA and fromB, assigns unique
     * elements to the non-unique indices thereof.
     * 
     * @param to Sum to
     * @param fromA Add shared from and copy A-specific from.
     * @param fromB Add shared from and copy B-specific from.
     */
    private void addDoublesFrom(double[] to, double[] fromA, double[] fromB) {
        to = (to == null) ? new double[nVarTot] : to;
        Arrays.fill(to, 0.0);
        for (int i = 0; i < nVarA; i++) {
            to[indexAToGlobal[i]] = fromA[i];
        }
        for (int i = 0; i < nVarB; i++) {
            to[indexBToGlobal[i]] += fromB[i];
        }
    }
    
    @Override
    public double energy(double[] x) {
        return energy(x, false);
    }

    @Override
    public double energy(double[] x, boolean verbose) {
        region.setX(x);
        region.setVerbose(verbose);
        try {
            team.execute(region);
        } catch (Exception ex) {
            throw new EnergyException(String.format(" Exception in calculating quad-topology energy: %s", ex.toString()), false);
        }

        if (verbose) {
            logger.info(String.format(" Total quad-topology energy: %12.4f", totalEnergy));
        }
        return totalEnergy;
    }
    
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        return energyAndGradient(x, g, false);
    }
    
    @Override
    public double energyAndGradient(double[] x, double[] g, boolean verbose) {
        region.setX(x);
        region.setG(g);
        region.setVerbose(verbose);
        try {
            team.execute(region);
        } catch (Exception ex) {
            throw new EnergyException(String.format(" Exception in calculating quad-topology energy: %s", ex.toString()), false);
        }

        if (verbose) {
            logger.info(String.format(" Total quad-topology energy: %12.4f", totalEnergy));
        }
        return totalEnergy;
    }
    
    @Override
    public void setScaling(double[] scaling) {
        this.scaling = scaling;
        if (scaling != null) {
            double[] scaleA = new double[nVarA];
            double[] scaleB = new double[nVarB];
            doublesTo(scaling, scaleA, scaleB);
            dualTopA.setScaling(scaleA);
            dualTopB.setScaling(scaleB);
        }
    }

    @Override
    public double[] getScaling() {
        return scaling;
    }

    @Override
    public double[] getCoordinates(double[] x) {
        dualTopA.getCoordinates(xA);
        dualTopB.getCoordinates(xB);
        doublesFrom(x, xA, xB);
        return x;
    }

    @Override
    public double[] getMass() {
        return mass;
    }

    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    @Override
    public int getNumberOfVariables() {
        return nVarTot;
    }
    
    /**
     * Returns number of shared variables.
     * @return Shared variables
     */
    public int getNumSharedVariables() {
        return nShared;
    }

    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        if (types == null) {
            VARIABLE_TYPE[] typesA = dualTopA.getVariableTypes();
            VARIABLE_TYPE[] typesB = dualTopB.getVariableTypes();
            if (typesA != null && typesB != null) {
                types = new VARIABLE_TYPE[nVarTot];
                copyFrom(types, dualTopA.getVariableTypes(), dualTopB.getVariableTypes());
            } else {
                logger.fine(" Variable types array remaining null due to null "
                        + "variable types in either A or B dual topology");
            }
        }
        return types;
    }

    @Override
    public void setVelocity(double[] velocity) {
        doublesTo(velocity, tempA, tempB);
        dualTopA.setVelocity(tempA);
        dualTopB.setVelocity(tempB);
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        doublesTo(acceleration, tempA, tempB);
        dualTopA.setVelocity(tempA);
        dualTopB.setVelocity(tempB);
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        doublesTo(previousAcceleration, tempA, tempB);
        dualTopA.setPreviousAcceleration(tempA);
        dualTopB.setPreviousAcceleration(tempB);
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        doublesFrom(velocity, dualTopA.getVelocity(tempA), dualTopB.getVelocity(tempB));
        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        doublesFrom(acceleration, dualTopA.getAcceleration(tempA), dualTopB.getAcceleration(tempB));
        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        doublesFrom(previousAcceleration, dualTopA.getPreviousAcceleration(tempA), dualTopB.getPreviousAcceleration(tempB));
        return previousAcceleration;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        dualTopA.setEnergyTermState(state);
        dualTopB.setEnergyTermState(state);
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setLambda(double lambda) {
        if (!Double.isFinite(lambda) || lambda > 1.0 || lambda < 0.0) {
            throw new ArithmeticException(String.format(" Attempted to set invalid lambda value of %10.6g", lambda));
        }
        this.lambda = lambda;
        dualTopA.setLambda(lambda);
        dualTopB.setLambda(lambda);
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        return dEdL;
    }

    /**
     * Returns true if both dual topologies are zero at the ends.
     *
     * @return If dEdL guaranteed zero at ends.
     */
    @Override
    public boolean dEdLZeroAtEnds() {
        if (!dualTopA.dEdLZeroAtEnds() || !dualTopB.dEdLZeroAtEnds()) {
            return false;
        }
        return true;
    }

    @Override
    public double getd2EdL2() {
        return d2EdL2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        dualTopA.getdEdXdL(tempA);
        dualTopB.getdEdXdL(tempB);
        addDoublesFrom(g, tempA, tempB);
    }
    
    public void setParallel(boolean parallel) {
        this.inParallel = parallel;
        if (team != null) {
            try {
                team.shutdown();
            } catch (Exception e) {
                logger.severe(String.format(" Exception in shutting down old ParallelTeam for DualTopologyEnergy: %s", e.toString()));
            }
        }
        team = parallel ? new ParallelTeam(2) : new ParallelTeam(1);
    }

    @Override
    public Crystal getCrystal() {
        return dualTopA.getCrystal();
    }

    @Override
    public void setCrystal(Crystal crystal) {
        dualTopA.setCrystal(crystal);
        dualTopB.setCrystal(crystal);
    }
    
    private class EnergyRegion extends ParallelRegion {
        
        private double[] x;
        private double[] g;
        private boolean gradient = false;
        
        private final EnergyASection sectA;
        private final EnergyBSection sectB;
        
        public EnergyRegion() {
            sectA = new EnergyASection();
            sectB = new EnergyBSection();
        }
        
        public void setX(double[] x) {
            this.x = x;
        }
        
        public void setG(double[] g) {
            this.g = g;
            setGradient(true);
        }
        
        public void setVerbose(boolean verbose) {
            sectA.setVerbose(verbose);
            sectB.setVerbose(verbose);
        }
        
        public void setGradient(boolean gradient) {
            this.gradient = gradient;
            sectA.setGradient(gradient);
            sectB.setGradient(gradient);
        }

        @Override
        public void start() throws Exception {
            doublesTo(x, xA, xB);
        }
        
        @Override
        public void run() throws Exception {
            execute(sectA, sectB);
        }
        
        @Override
        public void finish() throws Exception {
            totalEnergy = energyA + energyB;
            if (gradient) {
                addDoublesFrom(g, gA, gB);
                dEdL = dEdL_A + dEdL_B;
                d2EdL2 = d2EdL2_A + d2EdL2_B;
            }
            gradient = false;
        }
    }
    
    private class EnergyASection extends ParallelSection {

        private boolean verbose = false;
        private boolean gradient = false;
        
        @Override
        public void run() throws Exception {
            if (gradient) {
                energyA = dualTopA.energyAndGradient(xA, gA, verbose);
                dEdL_A = linterA.getdEdL();
                d2EdL2_A = linterA.getd2EdL2();
            } else {
                energyA = dualTopA.energy(xA, verbose);
            }
            this.verbose = false;
            this.gradient = false;
        }
        
        public void setVerbose(boolean verbose) {
            this.verbose = verbose;
        }
        
        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }
    }
    
    private class EnergyBSection extends ParallelSection {

        private boolean verbose = false;
        private boolean gradient = false;

        @Override
        public void run() throws Exception {
            if (gradient) {
                energyB = dualTopB.energyAndGradient(xB, gB, verbose);
                dEdL_B = linterB.getdEdL();
                d2EdL2_B = linterB.getd2EdL2();
            } else {
                energyB = dualTopB.energy(xB, verbose);
            }
            this.verbose = false;
            this.gradient = false;
        }
        
        public void setVerbose(boolean verbose) {
            this.verbose = verbose;
        }
        
        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }
    }
}

