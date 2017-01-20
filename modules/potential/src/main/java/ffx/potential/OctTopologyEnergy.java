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

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.Potential;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.utils.EnergyException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 *
 * @author jmlitman
 */
public class OctTopologyEnergy implements Potential, LambdaInterface {
    private final static Logger logger = Logger.getLogger(OctTopologyEnergy.class.getName());
    private final QuadTopologyEnergy quadTopGamma;
    private final QuadTopologyEnergy quadTopDelta;
    private final LambdaInterface linterGamma;
    private final LambdaInterface linterDelta;
    
    private final int nVarG;
    private final int nVarD;
    private final int nShared;
    private final int uniqueG;
    private final int uniqueD;
    private final int nVarTot;
    
    /**
     * Following arrays keep track of which index in the child QuadTopology is 
     * linked to which index in the OctTopology variable array.
     */
    private final int[] indexGToGlobal;
    private final int[] indexDToGlobal;
    private final int[] indexGlobalToG;
    private final int[] indexGlobalToD;
    
    private final double[] mass;
    private final double[] xG;
    private final double[] xD; // Laughing smiley not intended.
    private final double[] gG;
    private final double[] gD;
    /**
     * tempA and B arrays are used to hold the result of get methods applied
     * on dual topologies A and B; this saves the time of initializing new arrays
     * to hold the data for a very short time. Not thread-safe; if multiple
     * threads are to access this data, will need to refactor the class to 
     * synchronize on these arrays so only one method can use them at a time.
     */
    private final double[] tempG;
    private final double[] tempD;
    
    private Potential.STATE state = Potential.STATE.BOTH;
    private double lambda;
    private double totalEnergy;
    private double energyG;
    private double energyD;
    private double dEdL, dEdL_G, dEdL_D;
    private double d2EdL2, d2EdL2_G, d2EdL2_D;
    
    private boolean inParallel = false;
    private ParallelTeam team;
    private final EnergyRegion region;
    
    /**
     * Scaling and de-scaling will be applied inside QuadTopologyEnergy.
     */
    private double[] scaling;
    private Potential.VARIABLE_TYPE[] types = null;
    
    /**
     * General structure: this layer has gamma and delta quad topologies (often
     * shortened to G/D), with dual topologies A/B, with force field energies 1/2.
     * 
     * This constructor assumes there are no shared variables in either quad 
     * topology; this behavior is opposite QuadTopologyEnergy behavior.
     * 
     * @param quadTopG A QuadTopologyEnergy
     * @param quadTopD A QuadTopologyEnergy
     */      
    public OctTopologyEnergy(QuadTopologyEnergy quadTopG, QuadTopologyEnergy quadTopD) {
        this(quadTopG, quadTopD, getVars(quadTopG, true), getVars(quadTopD, true));
    }
    
    /**
     * General structure: this layer has gamma and delta quad topologies (often
     * shortened to G/D), with dual topologies A/B, with force field energies 1/2.
     * 
     * This constructor permits QuadTopologyEnergy-like pin-all behavior with
     * uniqueAll true.
     * 
     * @param quadTopG A QuadTopologyEnergy
     * @param quadTopD A QuadTopologyEnergy
     * @param uniqueAll All unique (else all shared)
     */      
    public OctTopologyEnergy(QuadTopologyEnergy quadTopG, QuadTopologyEnergy quadTopD, boolean uniqueAll) {
        this(quadTopG, quadTopD, getVars(quadTopG, uniqueAll), getVars(quadTopD, uniqueAll));
    }  
    
    /**
     * Utility method to return either all variable indices in the quad topology,
     * or null (for sharing none or all variables, respectively).
     * @param quadTop
     * @param getAll
     * @return Unique variable indices
     */
    private static List<Integer> getVars(QuadTopologyEnergy quadTop, boolean getAll) {
        if (getAll) {
            int nVars = quadTop.getNumberOfVariables();
            return IntStream.range(0, nVars).boxed().collect(Collectors.toList());
        } else {
            return null;
        }
    }
    
    /**
     * General structure: this layer has gamma and delta quad topologies (often
     * shortened to G/D), with dual topologies A/B, with force field energies 1/2.
     * 
     * @param quadTopG A QuadTopologyEnergy
     * @param quadTopD A QuadTopologyEnergy
     * @param uniqueGList Variables unique to gamma
     * @param uniqueDList Variables unique to delta
     */        
    public OctTopologyEnergy(QuadTopologyEnergy quadTopG, QuadTopologyEnergy quadTopD, List<Integer> uniqueGList, List<Integer> uniqueDList) {
        this.quadTopGamma = quadTopG;
        this.quadTopDelta = quadTopD;
        //quadTopDelta.reloadCommonMasses(true);
        linterGamma = (LambdaInterface) quadTopG;
        linterDelta = (LambdaInterface) quadTopD;
        nVarG = quadTopGamma.getNumberOfVariables();
        nVarD = quadTopDelta.getNumberOfVariables();
        
        /**
         * Following logic is to A, deal with Integer to int array
         * conversion issues, and B, ensure the set of unique variables is a 
         * unique set. Unique elements may come from either the provided list, or
         * from the non-shared elements of the dual topologies.
         */
        Set<Integer> uniqueGSet = new LinkedHashSet<>();
        if (uniqueGList != null) {
            uniqueGSet.addAll(uniqueGList);
        }
        int nCommon = quadTopGamma.getNumSharedVariables();
        for (int i = nCommon; i < nVarG; i++) {
            uniqueGSet.add(i);
        }
        uniqueGList = new ArrayList(uniqueGSet);
        uniqueG = uniqueGList.size();
        int[] uniquesG = new int[uniqueG];
        for (int i = 0; i < uniqueG; i++) {
            uniquesG[i] = uniqueGList.get(i);
        }
        /*IntStream.range(0, uniqueA).parallel().forEach((i) -> {
            uniquesA[i] = uniqueAList.get(i);
        });*/
        
        /**
         * Following logic is to A, deal with Integer to int array
         * conversion issues, and B, ensure the set of unique variables is a 
         * unique set. Unique elements may come from either the provided list, or
         * from the non-shared elements of the dual topologies.
         */
        Set<Integer> uniqueDSet = new LinkedHashSet<>();
        if (uniqueDList != null) {
            uniqueDSet.addAll(uniqueDList);
        }
        nCommon = quadTopDelta.getNumSharedVariables();
        for (int i = nCommon; i < nVarD; i++) {
            uniqueDSet.add(i);
        }
        uniqueDList = new ArrayList(uniqueDSet);
        uniqueD = uniqueDList.size();
        int[] uniquesD = new int[uniqueD];
        for (int i = 0; i < uniqueD; i++) {
            uniquesD[i] = uniqueDList.get(i);
        }
        /*IntStream.range(0, uniqueA).parallel().forEach((i) -> {
            uniquesA[i] = uniqueAList.get(i);
        });*/
        
        nShared = nVarG - uniqueG;
        assert (nShared == nVarD - uniqueD);
        nVarTot = nShared + uniqueG + uniqueD;
        
        indexGToGlobal = new int[nVarG];
        indexDToGlobal = new int[nVarD];
        indexGlobalToG = new int[nVarTot];
        indexGlobalToD = new int[nVarTot];
        // -1 indicates this index in the global array does not point to one in 
        // the dual topology.
        Arrays.fill(indexGlobalToG, -1);
        Arrays.fill(indexGlobalToD, -1);
        
        if (uniqueG > 0) {
            int commonIndex = 0;
            int uniqueIndex = 0;
            for (int i = 0; i < nVarG; i++) {
                if (uniqueIndex < uniqueG && i == uniquesG[uniqueIndex]) {
                    int destIndex = nShared + uniqueIndex;
                    indexGToGlobal[i] = destIndex;
                    indexGlobalToG[destIndex] = i;
                    //logger.info(String.format("Unique A: i %d destIndex %d uniqueIndex %d", i, destIndex, uniqueIndex));
                    ++uniqueIndex;
                } else {
                    indexGToGlobal[i] = commonIndex;
                    indexGlobalToG[commonIndex++] = i;
                    //logger.info(String.format("Common A: i %d commonIndex %d", i, commonIndex-1));
                }
            }
        } else {
            for (int i = 0; i < nVarG; i++) {
                indexGToGlobal[i] = i;
                indexGlobalToG[i] = i;
            }
            /*IntStream.range(0, nVarA).parallel().forEach((i) -> { 
                indexAToGlobal[i] = i;
                indexGlobalToA[i] = i; 
            });*/
        }
        
        if (uniquesD.length > 0) {
            int commonIndex = 0;
            int uniqueIndex = 0;
            for (int i = 0; i < nVarD; i++) {
                if (uniqueIndex < uniqueD && i == uniquesD[uniqueIndex]) {
                    int destIndex = nVarG + uniqueIndex;
                    indexDToGlobal[i] = destIndex;
                    indexGlobalToD[destIndex] = i;
                    ++uniqueIndex;
                } else {
                    indexDToGlobal[i] = commonIndex;
                    indexGlobalToD[commonIndex++] = i;
                }
            }
        } else {
            /*System.out.println(String.format(" Lengths: Total: %5d Shared: %5d "
                    + "Total A: %5d Unique A: %5d Total B: %5d "
                    + "Unique B: %5d", nVarTot, nShared, nVarA, uniqueA, nVarB, uniqueB));
            System.out.println(String.format(" Arrays:  GlobalToA: %5d "
                    + "AToGlobal: %5d GlobalToB: %5d BToGlobal: %5d", 
                    indexGlobalToA.length, indexAToGlobal.length, indexGlobalToB.length, indexBToGlobal.length));*/
            for (int i = 0; i < nVarD; i++) {
                indexDToGlobal[i] = i;
                indexGlobalToD[i] = i;
            }
            /*IntStream.range(0, nVarB).parallel().forEach((i) -> {
                indexBToGlobal[i] = i;
                indexGlobalToB[i] = i;
            });*/
        }
        
        xG = new double[nVarG];
        xD = new double[nVarD];
        gG = new double[nVarG];
        gD = new double[nVarD];
        tempG = new double[nVarG];
        tempD = new double[nVarD];
        mass = new double[nVarTot];
        doublesFrom(mass, quadTopGamma.getMass(), quadTopDelta.getMass());
        
        region = new EnergyRegion();
        team = new ParallelTeam(1);
    }
    
    /**
     * Copies from an object array of length nVarTot to two object arrays of 
     * length nVarG and nVarD.
     * 
     * @param <T> Type of object
     * @param from Copy from
     * @param toG Copy shared and gamma-specific to
     * @param toD Copy shared and delta-specific to
     */
    private <T> void copyTo(T[] from, T[] toG, T[] toD) {
        if (toG == null) {
            toG = Arrays.copyOf(from, nVarG);
        }
        if (toD == null) {
            toD = Arrays.copyOf(from, nVarD);
        }
        for (int i = 0; i < nVarTot; i++) {
            int index = indexGlobalToG[i];
            if (index >= 0) {
                toG[index] = from[i];
            }
            index = indexGlobalToD[i];
            if (index >= 0) {
                toD[index] = from[i];
            }
        }
    }

    /**
     * Copies from object arrays of length nVarG and nVarD to an object array of
     * length nVarTot; asserts objects in common indices are equal. Should not
     * be used when the result of the common indices should be f(G,D)
     * 
     * @param <T> Type of object
     * @param to Copy to
     * @param fromG Copy shared and gamma-specific from
     * @param fromD Copy delta-specific from
     */
    private <T> void copyFrom(T[] to, T[] fromG, T[] fromD) {
        if (to == null) {
            to = Arrays.copyOf(fromG, nVarTot);
        }
        for (int i = 0; i < nVarG; i++) {
            int index = indexGToGlobal[i];
            to[index] = fromG[i];
        }
        for (int i = 0; i < nVarD; i++) {
            int index = indexDToGlobal[i];
            // Assert either a unique variable, or equals what was added from A.
            assert (index >= nShared || to[index].equals(fromD[i]));
            to[index] = fromD[i];
        }
    }
    
    /**
     * Copies from an double array of length nVarTot to two double arrays of 
     * length nVarG and nVarD.
     * 
     * @param <T> Type of object
     * @param from Copy from
     * @param toG Copy shared and gamma-specific to
     * @param toD Copy shared and delta-specific to
     */
    private void doublesTo(double[] from, double[] toG, double[] toD) {
        toG = (toG == null) ? new double[nVarG] : toG;
        toD = (toD == null) ? new double[nVarD] : toD;
        for (int i = 0; i < nVarTot; i++) {
            int index = indexGlobalToG[i];
            if (index >= 0) {
                toG[index] = from[i];
            }
            index = indexGlobalToD[i];
            if (index >= 0) {
                toD[index] = from[i];
            }
        }
    }

    /**
     * Copies from double arrays of length nVarG and nVarD to an object array of
     * length nVarTot; asserts common indices are equal. Should not be used when 
     * the result of the common indices should be f(G,D)
     * 
     * @param to Copy to
     * @param fromG Copy shared and gamma-specific from
     * @param fromD Copy delta-specific from
     */
    private void doublesFrom(double[] to, double[] fromG, double[] fromD) {
        to = (to == null) ? new double[nVarTot] : to;
        for (int i = 0; i < nVarG; i++) {
            to[indexGToGlobal[i]] = fromG[i];
        }
        for (int i = 0; i < nVarD; i++) {
            int index = indexDToGlobal[i];
            // Assert this is either a unique from B or it's equal to what came from A.
            assert (index >= nShared || to[index] == fromD[i]);
            to[index] = fromD[i];
        }
    }
    
    private void addDoublesFrom(double[] to, double[] fromG, double[] fromD) {
        addDoublesFrom(to, fromG, fromD, 1.0);
    }

    /**
     * Assigns common indices of to to be sum of fromG and fromD, assigns unique
     * elements to the non-unique indices thereof. Additionally scales the results
     * by some scalar multiplier.
     * 
     * @param to Sum to
     * @param fromG Add shared from and copy gamma-specific from.
     * @param fromD Add shared from and copy delta-specific from.
     * @param scalingFactor Scale values by this factor.
     */
    private void addDoublesFrom(double[] to, double[] fromG, double[] fromD, double scalingFactor) {
        to = (to == null) ? new double[nVarTot] : to;
        Arrays.fill(to, 0.0);
        for (int i = 0; i < nVarG; i++) {
            to[indexGToGlobal[i]] = fromG[i] * scalingFactor;
        }
        for (int i = 0; i < nVarD; i++) {
            to[indexDToGlobal[i]] += (fromD[i] * scalingFactor);
        }
    }
    
    /**
     * Assigns common indices of to to be difference of fromG and fromD, assigns 
     * unique elements to the non-unique indices thereof (multiplied by -1 for
     * delta).
     * 
     * @param to Sum to
     * @param fromG Add shared from and copy gamma-specific from.
     * @param fromD Subtract shared from and copy -1 * delta-specific from.
     */
    private void subtractDoublesFrom(double[] to, double[] fromG, double[] fromD) {
        subtractDoublesFrom(to, fromG, fromD, 1.0);
    }
    
    /**
     * Assigns common indices of to to be difference of fromG and fromD, assigns 
     * unique elements to the non-unique indices thereof (multiplied by -1 for
     * delta).
     * 
     * @param to Sum to
     * @param fromG Add shared from and copy gamma-specific from.
     * @param fromD Subtract shared from and copy -1 * delta-specific from.
     * @param scalingFactor Scale values by this factor
     */
    private void subtractDoublesFrom(double[] to, double[] fromG, double[] fromD, double scalingFactor) {
        to = (to == null) ? new double[nVarTot] : to;
        Arrays.fill(to, 0.0);
        for (int i = 0; i < nVarG; i++) {
            to[indexGToGlobal[i]] = fromG[i] * scalingFactor;
        }
        for (int i = 0; i < nVarD; i++) {
            to[indexDToGlobal[i]] -= (fromD[i] * scalingFactor);
        }
    }
    
    @Override
    public double energy(double[] x) {
        return energy(x, false);
    }

    @Override
    public double energy(double[] x, boolean verbose) {
        //if (inParallel) {
        region.setX(x);
        region.setVerbose(verbose);
        try {
            team.execute(region);
        } catch (Exception ex) {
            throw new EnergyException(String.format(" Exception in calculating quad-topology energy: %s", ex.toString()), false);
        }
        /*} else {
            doublesTo(x, xA, xB);
            energyA = dualTopA.energy(xA, verbose);
            energyB = dualTopB.energy(xB, verbose);
            totalEnergy = energyA + energyB;
        }*/
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
        //if (inParallel) {
        region.setX(x);
        region.setG(g);
        region.setVerbose(verbose);
        try {
            team.execute(region);
        } catch (Exception ex) {
            throw new EnergyException(String.format(" Exception in calculating quad-topology energy: %s", ex.toString()), false);
        }
        /*} else {
            doublesTo(x, xA, xB);

            energyA = dualTopA.energyAndGradient(xA, gA, verbose);
            dEdL_A = linterA.getdEdL();
            d2EdL2_A = linterA.getd2EdL2();

            energyB = dualTopB.energyAndGradient(xB, gB, verbose);
            dEdL_B = linterB.getdEdL();
            d2EdL2_B = linterB.getd2EdL2();

            subtractDoublesFrom(g, gA, gB);

            dEdL = dEdL_A - dEdL_B;
            d2EdL2 = d2EdL2_A - d2EdL2_B;
            totalEnergy = energyA - energyB;
        }*/
        if (verbose) {
            logger.info(String.format(" Total quad-topology energy: %12.4f", totalEnergy));
        }
        return totalEnergy;
    }
    
    @Override
    public void setScaling(double[] scaling) {
        this.scaling = scaling;
        if (scaling != null) {
            double[] scaleG = new double[nVarG];
            double[] scaleD = new double[nVarD];
            doublesTo(scaling, scaleG, scaleD);
            quadTopGamma.setScaling(scaleG);
            quadTopDelta.setScaling(scaleD);
        }
    }

    @Override
    public double[] getScaling() {
        return scaling;
    }

    @Override
    public double[] getCoordinates(double[] x) {
        quadTopGamma.getCoordinates(xG);
        quadTopDelta.getCoordinates(xD);
        doublesFrom(x, xG, xD);
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

    @Override
    public Potential.VARIABLE_TYPE[] getVariableTypes() {
        if (types == null) {
            Potential.VARIABLE_TYPE[] typesG = quadTopGamma.getVariableTypes();
            Potential.VARIABLE_TYPE[] typesD = quadTopDelta.getVariableTypes();
            if (typesG != null && typesD != null) {
                types = new Potential.VARIABLE_TYPE[nVarTot];
                copyFrom(types, quadTopGamma.getVariableTypes(), quadTopDelta.getVariableTypes());
            } else {
                logger.fine(" Variable types array remaining null due to null "
                        + "variable types in either A or B dual topology");
            }
        }
        return types;
    }

    @Override
    public void setVelocity(double[] velocity) {
        doublesTo(velocity, tempG, tempD);
        quadTopGamma.setVelocity(tempG);
        quadTopDelta.setVelocity(tempD);
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        doublesTo(acceleration, tempG, tempD);
        quadTopGamma.setVelocity(tempG);
        quadTopDelta.setVelocity(tempD);
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        doublesTo(previousAcceleration, tempG, tempD);
        quadTopGamma.setPreviousAcceleration(tempG);
        quadTopDelta.setPreviousAcceleration(tempD);
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        doublesFrom(velocity, quadTopGamma.getVelocity(tempG), quadTopDelta.getVelocity(tempD));
        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        doublesFrom(acceleration, quadTopGamma.getAcceleration(tempG), quadTopDelta.getAcceleration(tempD));
        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        doublesFrom(previousAcceleration, quadTopGamma.getPreviousAcceleration(tempG), quadTopDelta.getPreviousAcceleration(tempD));
        return previousAcceleration;
    }

    @Override
    public void setEnergyTermState(Potential.STATE state) {
        this.state = state;
        quadTopGamma.setEnergyTermState(state);
        quadTopDelta.setEnergyTermState(state);
    }

    @Override
    public Potential.STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setLambda(double lambda) {
        this.lambda = lambda;
        double lambdaG = 0.5 * lambda;
        quadTopGamma.setLambda(lambdaG);
        quadTopDelta.setLambda((1.0 - lambdaG));
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        return dEdL;
    }

    @Override
    public double getd2EdL2() {
        return d2EdL2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        quadTopGamma.getdEdXdL(tempG);
        quadTopDelta.getdEdXdL(tempD);
        addDoublesFrom(g, tempG, tempD, 0.5);
        //subtractDoublesFrom(g, tempG, tempD, 0.5);
    }
    
    public void setParallel(boolean parallel) {
        this.inParallel = parallel;
        if (team != null) {
            try {
                team.shutdown();
            } catch (Exception e) {
                logger.severe(String.format(" Exception in shutting down old ParallelTeam for QuadTopologyEnergy: %s", e.toString()));
            }
        }
        team = parallel ? new ParallelTeam(2) : new ParallelTeam(1);
    }
    
    private class EnergyRegion extends ParallelRegion {
        
        private double[] x;
        private double[] g;
        private boolean gradient = false;
        
        private final EnergyGSection sectG;
        private final EnergyDSection sectD;
        
        public EnergyRegion() {
            sectG = new EnergyGSection();
            sectD = new EnergyDSection();
        }
        
        public void setX(double[] x) {
            this.x = x;
        }
        
        public void setG(double[] g) {
            this.g = g;
            setGradient(true);
        }
        
        public void setVerbose(boolean verbose) {
            sectG.setVerbose(verbose);
            sectD.setVerbose(verbose);
        }
        
        public void setGradient(boolean gradient) {
            this.gradient = gradient;
            sectG.setGradient(gradient);
            sectD.setGradient(gradient);
        }

        @Override
        public void start() throws Exception {
            doublesTo(x, xG, xD);
        }
        
        @Override
        public void run() throws Exception {
            execute(sectG, sectD);
        }
        
        @Override
        public void finish() throws Exception {
            totalEnergy = energyG - energyD;
            if (gradient) {
                subtractDoublesFrom(g, gG, gD);
                dEdL = 0.5 * (dEdL_G + dEdL_D);
                d2EdL2 = 0.25 * (d2EdL2_G - d2EdL2_D);
            }
            gradient = false;
        }
    }
    
    private class EnergyGSection extends ParallelSection {

        private boolean verbose = false;
        private boolean gradient = false;
        
        @Override
        public void run() throws Exception {
            if (gradient) {
                energyG = quadTopGamma.energyAndGradient(xG, gG, verbose);
                dEdL_G = linterGamma.getdEdL();
                d2EdL2_G = linterGamma.getd2EdL2();
            } else {
                energyG = quadTopGamma.energy(xG, verbose);
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
    
    private class EnergyDSection extends ParallelSection {

        private boolean verbose = false;
        private boolean gradient = false;

        @Override
        public void run() throws Exception {
            if (gradient) {
                energyD = quadTopDelta.energyAndGradient(xD, gD, verbose);
                dEdL_D = linterDelta.getdEdL();
                d2EdL2_D = linterDelta.getd2EdL2();
            } else {
                energyD = quadTopDelta.energy(xD, verbose);
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
