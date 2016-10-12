/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import ffx.numerics.Potential;
import ffx.potential.bonded.LambdaInterface;
import java.util.logging.Logger;

/**
 * Implements an error-canceling quad topology, where two large dual-topology 
 * simulation legs are run simultaneously to arrive at a small sum. This is a
 * simplistic implementation running both dual topologies linked only by a common
 * lambda; will be likely be refactored into one implementation of an interface
 * or abstract class when more sophisticated implementations are added.
 * 
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class QuadTopologyEnergy implements Potential, LambdaInterface {
    private final static Logger logger = Logger.getLogger(QuadTopologyEnergy.class.getName());
    private final DualTopologyEnergy dualTopA;
    private final DualTopologyEnergy dualTopB;
    private final LambdaInterface linterA;
    private final LambdaInterface linterB;
    private final int nVarA;
    private final int nVarB;
    private final int nVarTot;
    private final double[] mass;
    private STATE state = STATE.BOTH;
    private double lambda;
    private double totalEnergy;
    private double energyA;
    private double energyB;
    private double dEdL, dEdL_A, dEdL_B;
    private double d2EdL2, d2EdL2_A, d2EdL2_B;
    private double[] xA;
    private double[] xB;
    private double[] gA;
    private double[] gB;
    /**
     * tempA and B arrays are used to hold the result of get methods applied
     * on dual topologies A and B; this saves the time of initializing new arrays
     * to hold the data for a very short time. Not thread-safe; if multiple
     * threads are to access this data, will need to refactor the class to 
     * synchronize on these arrays so only one method can use them at a time.
     */
    private final double[] tempA;
    private final double[] tempB;
    /**
     * Scaling and de-scaling will be applied inside DualTopologyEnergy. This
     * array will just be a concatenation of the two scaling arrays.
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
     * @param dualTopologyA A DualTopologyEnergy
     * @param dualTopologyB A DualTopologyEnergy
     */    
    public QuadTopologyEnergy(DualTopologyEnergy dualTopologyA, DualTopologyEnergy dualTopologyB) {
        this.dualTopA = dualTopologyA;
        this.dualTopB = dualTopologyB;
        linterA = (LambdaInterface) dualTopologyA;
        linterB = (LambdaInterface) dualTopologyB;
        nVarA = dualTopA.getNumberOfVariables();
        nVarB = dualTopB.getNumberOfVariables();
        nVarTot = nVarA + nVarB;
        mass = new double[nVarTot];
        System.arraycopy(dualTopA.getMass(), 0, mass, 0, nVarA);
        System.arraycopy(dualTopB.getMass(), 0, mass, nVarA, nVarB);
        xA = new double[nVarA];
        xB = new double[nVarB];
        gA = new double[nVarA];
        gB = new double[nVarB];
        tempA = new double[nVarA];
        tempB = new double[nVarB];
        
    }

    @Override
    public double energy(double[] x) {
        splitCoordinates(x);
        energyA = dualTopA.energy(xA);
        energyB = dualTopB.energy(xB);
        totalEnergy = energyA + energyB;
        unsplitCoordinates(x);
        return totalEnergy;
    }

    @Override
    public double energyAndGradient(double[] x, double[] g) {
        splitCoordinates(x);
        
        energyA = dualTopA.energyAndGradient(xA, gA);
        dEdL_A = linterA.getdEdL();
        d2EdL2_A = linterA.getd2EdL2();

        energyB = dualTopB.energyAndGradient(xB, gB);
        dEdL_B = linterB.getdEdL();
        d2EdL2_B = linterB.getd2EdL2();
        
        unsplitGradient(x, g);
        dEdL = dEdL_A + dEdL_B;
        d2EdL2 = d2EdL2_A + d2EdL2_B;
        totalEnergy = energyA + energyB;
        return totalEnergy;
    }

    @Override
    public void setScaling(double[] scaling) {
        this.scaling = scaling;
        if (scaling != null) {
            double[] scaleA = new double[nVarA];
            System.arraycopy(scaling, 0, scaleA, 0, nVarA);
            dualTopA.setScaling(scaleA);

            double[] scaleB = new double[nVarB];
            System.arraycopy(scaling, nVarA, scaleB, 0, nVarB);
            dualTopB.setScaling(scaleB);
        }
    }

    @Override
    public double[] getScaling() {
        return scaling;
    }

    /**
     * Reduce xA and xB into x.
     * @param x 
     */
    private void unsplitCoordinates(double x[]) {
        System.arraycopy(xA, 0, x, 0, nVarA);
        System.arraycopy(xB, 0, x, nVarA, nVarB);
    }

    /**
     * Reduce xA and xB into x, gA and gB into g.
     * @param x
     * @param g 
     */
    private void unsplitGradient(double x[], double g[]) {
        unsplitCoordinates(x);
        if (g == null) {
            g = new double[nVarTot];
        }
        System.arraycopy(gA, 0, g, 0, nVarA);
        System.arraycopy(gB, 0, g, nVarA, nVarB);
    }

    /**
     * Split x into xA and xB.
     * @param x 
     */
    private void splitCoordinates(double x[]) {
        System.arraycopy(x, 0, xA, 0, nVarA);
        System.arraycopy(x, nVarA, xB, 0, nVarB);
    }

    @Override
    public double[] getCoordinates(double[] x) {
        if (x == null) {
            x = new double[nVarTot];
        }
        dualTopA.getCoordinates(xA);
        dualTopB.getCoordinates(xB);
        System.arraycopy(xA, 0, x, 0, nVarA);
        System.arraycopy(xB, 0, x, nVarA, nVarB);
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
    public VARIABLE_TYPE[] getVariableTypes() {
        if (types == null) {
            types = new VARIABLE_TYPE[nVarTot];
            System.arraycopy(dualTopA.getVariableTypes(), 0, types, 0, nVarA);
            System.arraycopy(dualTopB.getVariableTypes(), 0, types, nVarA, nVarB);
        }
        return types;
    }

    @Override
    public void setVelocity(double[] velocity) {
        double[] velSlice = new double[nVarA];
        System.arraycopy(velocity, 0, velSlice, 0, nVarA);
        dualTopA.setVelocity(velSlice);
        
        velSlice = new double[nVarB];
        System.arraycopy(velocity, nVarA, velSlice, 0, nVarB);
        dualTopB.setVelocity(velSlice);
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        double[] accSlice = new double[nVarA];
        System.arraycopy(acceleration, 0, accSlice, 0, nVarA);
        dualTopA.setAcceleration(accSlice);
        
        accSlice = new double[nVarB];
        System.arraycopy(acceleration, nVarA, accSlice, 0, nVarB);
        dualTopB.setAcceleration(accSlice);
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        double[] priorAccSlice = new double[nVarA];
        System.arraycopy(previousAcceleration, 0, priorAccSlice, 0, nVarA);
        dualTopA.setPreviousAcceleration(priorAccSlice);
        
        priorAccSlice = new double[nVarB];
        System.arraycopy(previousAcceleration, nVarA, priorAccSlice, 0, nVarB);
        dualTopB.setPreviousAcceleration(priorAccSlice);
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        if (velocity == null) {
            velocity = new double[nVarTot];
        }
        System.arraycopy(dualTopA.getVelocity(tempA), 0, velocity, 0, nVarA);
        System.arraycopy(dualTopB.getVelocity(tempB), 0, velocity, nVarA, nVarB);
        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        if (acceleration == null) {
            acceleration = new double[nVarTot];
        }
        System.arraycopy(dualTopA.getAcceleration(tempA), 0, acceleration, 0, nVarA);
        System.arraycopy(dualTopB.getAcceleration(tempB), 0, acceleration, nVarA, nVarB);
        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        if (previousAcceleration == null) {
            previousAcceleration = new double[nVarTot];
        }
        System.arraycopy(dualTopA.getPreviousAcceleration(tempA), 0, previousAcceleration, 0, nVarA);
        System.arraycopy(dualTopB.getPreviousAcceleration(tempB), 0, previousAcceleration, nVarA, nVarB);
        return previousAcceleration;
    }

    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setLambda(double lambda) {
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

    @Override
    public double getd2EdL2() {
        return d2EdL2;
    }

    @Override
    public void getdEdXdL(double[] g) {
        dualTopA.getdEdXdL(tempA);
        dualTopB.getdEdXdL(tempB);
        if (g == null) {
            g = new double[nVarTot];
        }
        System.arraycopy(tempA, 0, g, 0, nVarA);
        System.arraycopy(tempB, 0, g, nVarA, nVarB);
    }
}
