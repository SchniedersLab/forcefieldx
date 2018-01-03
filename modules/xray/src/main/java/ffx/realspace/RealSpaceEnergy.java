/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.realspace;

import java.util.logging.Logger;

import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.xray.RefinementMinimize.RefinementMode;
import ffx.xray.RefinementModel;

/**
 * Combine the Real Space target and chemical potential energy.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 *
 * @since 1.0
 */
public class RealSpaceEnergy implements LambdaInterface, Potential {

    private static final Logger logger = Logger.getLogger(RealSpaceEnergy.class.getName());

    /**
     * The Real Space data to refine against.
     */
    private final RealSpaceData realSpaceData;
    /**
     * The Refinement Model that contains info on mapping between alternate
     * conformers.
     */
    private final RefinementModel refinementModel;
    /**
     * The number of parameters that are being refined.
     */
    private int nXYZ;
    /**
     * The refinement mode.
     */
    private RefinementMode refinementMode;
    /**
     * If true, the XYZ coordinates will be refined.
     */
    private boolean refineXYZ = false;
    /**
     * Optimization scaling used to improve convergence.
     */
    protected double[] optimizationScaling = null;
    /**
     * Value of the lambda state variable.
     */
    protected double lambda = 1.0;
    /**
     * Total energy of the refinement.
     */
    private double totalEnergy;
    /**
     * The RealSpaceEnergy is updated with FAST varying energy terms.
     */
    private STATE state = STATE.BOTH;

    /**
     * Diffraction data energy target
     *
     * @param realSpaceData {@link RealSpaceData} object to associate with the
     * target
     * @param nxyz number of xyz parameters
     * @param nb number of b factor parameters
     * @param nocc number of occupancy parameters
     * @param refinementMode the {@link RefinementMinimize.RefinementMode} type
     * of refinement requested
     */
    public RealSpaceEnergy(RealSpaceData realSpaceData, int nxyz, int nb, int nocc,
            RefinementMode refinementMode) {
        this.realSpaceData = realSpaceData;
        this.refinementModel = realSpaceData.getRefinementModel();
        this.refinementMode = refinementMode;
        this.nXYZ = nxyz;
        setRefinementBooleans();
    }

    /**
     * <p>
     * Getter for the field <code>refinementMode</code>.</p>
     *
     * @return a {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public RefinementMode getRefinementMode() {
        return refinementMode;
    }

    /**
     * <p>
     * Setter for the field <code>refinementMode</code>.</p>
     *
     * @param refinementMode a
     * {@link ffx.xray.RefinementMinimize.RefinementMode} object.
     */
    public void setRefinementMode(RefinementMode refinementMode) {
        this.refinementMode = refinementMode;
        setRefinementBooleans();
    }

    /**
     * If the refinement mode has changed, this should be called to update which
     * parameters are being fit.
     */
    private void setRefinementBooleans() {
        // reset, if previously set
        refineXYZ = false;

        if (refinementMode == RefinementMode.COORDINATES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS
                || refinementMode == RefinementMode.COORDINATES_AND_OCCUPANCIES
                || refinementMode == RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES) {
            refineXYZ = true;
        }
    }

    /**
     * The parameters passed in are only for "active" atoms.
     *
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {

        double e = 0.0;
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        if (refineXYZ) {
            int index = 0;
            double xyz[] = new double[3];
            for (Atom a : refinementModel.getTotalAtomArray()) {
                if (a.isActive()) {
                    int i = index * 3;
                    xyz[0] = x[i];
                    xyz[1] = x[i + 1];
                    xyz[2] = x[i + 2];
                    a.setXYZ(xyz);
                    a.setXYZGradient(0.0, 0.0, 0.0);
                    a.setLambdaXYZGradient(0.0, 0.0, 0.0);
                    index++;
                }
            }
        }

        // Target function for real space refinement.
        if (refineXYZ) {
            e = realSpaceData.computeRealSpaceTarget();
        }

        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }

        totalEnergy = e;
        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        double e = 0.0;
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }

        if (refineXYZ) {
            int index = 0;
            double xyz[] = new double[3];
            for (Atom a : refinementModel.getTotalAtomArray()) {
                if (a.isActive()) {
                    int i = index * 3;
                    xyz[0] = x[i];
                    xyz[1] = x[i + 1];
                    xyz[2] = x[i + 2];
                    a.setXYZ(xyz);
                    a.setXYZGradient(0.0, 0.0, 0.0);
                    a.setLambdaXYZGradient(0.0, 0.0, 0.0);
                    index++;
                }
            }
        }

        /**
         * Target function for real space refinement
         */
        if (refineXYZ) {
            e = realSpaceData.computeRealSpaceTarget();
            /**
             * Pack gradients into gradient array
             */
            getXYZGradients(g);
        }

        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }
        totalEnergy = e;
        return e;
    }

    /**
     * Get the number of xyz parameters being fit.
     *
     * @return the number of xyz parameters
     */
    public int getNXYZ() {
        return nXYZ;
    }

    /**
     * Set the number of xyz parameters.
     *
     * @param nxyz requested number of xyz parameters
     */
    public void setNXYZ(int nxyz) {
        this.nXYZ = nxyz;
    }

    /**
     * Fill gradient array with xyz gradient.
     *
     * @param g array to add gradient to
     */
    public void getXYZGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                a.getXYZGradient(grad);
                g[index++] = grad[0];
                g[index++] = grad[1];
                g[index++] = grad[2];
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double x[]) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        int index = 0;
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                x[index++] = a.getX();
                x[index++] = a.getY();
                x[index++] = a.getZ();
            }
        }
        return x;
    }

    /**
     * Set atomic coordinates positions.
     *
     * @param x an array of coordinates for active atoms.
     */
    public void setCoordinates(double x[]) {
        assert (x != null);
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                xyz[0] = x[index++];
                xyz[1] = x[index++];
                xyz[2] = x[index++];
                a.moveTo(xyz);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling
    ) {
        optimizationScaling = scaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        double mass[] = new double[nXYZ];
        int i = 0;
        if (refineXYZ) {
            for (Atom a : refinementModel.getTotalAtomArray()) {
                if (a.isActive()) {
                    double m = a.getMass();
                    mass[i++] = m;
                    mass[i++] = m;
                    mass[i++] = m;
                }
            }
        }
        return mass;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return nXYZ;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda
    ) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            realSpaceData.setLambda(lambda);
        } else {
            String message = String.format(" Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        return realSpaceData.getdEdL();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        return 0.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradient
    ) {
        realSpaceData.getdEdXdL(gradient);
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {

        int nActive = 0;
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                nActive++;
            }
        }

        VARIABLE_TYPE type[] = new VARIABLE_TYPE[nActive * 3];

        int index = 0;
        for (int i = 0; i < nActive; i++) {
            type[index++] = VARIABLE_TYPE.X;
            type[index++] = VARIABLE_TYPE.Y;
            type[index++] = VARIABLE_TYPE.Z;
        }
        return type;
    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setEnergyTermState(STATE state
    ) {
        this.state = state;
    }

    @Override
    public void setVelocity(double[] velocity
    ) {
        if (velocity == null) {
            return;
        }
        int index = 0;
        double vel[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                vel[0] = velocity[index++];
                vel[1] = velocity[index++];
                vel[2] = velocity[index++];
                a.setVelocity(vel);
            }
        }
    }

    @Override
    public void setAcceleration(double[] acceleration
    ) {
        if (acceleration == null) {
            return;
        }
        int index = 0;
        double accel[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                accel[0] = acceleration[index++];
                accel[1] = acceleration[index++];
                accel[2] = acceleration[index++];
                a.setAcceleration(accel);
            }
        }
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration
    ) {
        if (previousAcceleration == null) {
            return;
        }
        int index = 0;
        double prev[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                prev[0] = previousAcceleration[index++];
                prev[1] = previousAcceleration[index++];
                prev[2] = previousAcceleration[index++];
                a.setPreviousAcceleration(prev);
            }
        }
    }

    @Override
    public double[] getVelocity(double[] velocity
    ) {
        int n = getNumberOfVariables();
        if (velocity == null || velocity.length < n) {
            velocity = new double[n];
        }
        int index = 0;
        double v[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                a.getVelocity(v);
                velocity[index++] = v[0];
                velocity[index++] = v[1];
                velocity[index++] = v[2];
            }
        }
        return velocity;
    }

    @Override
    public double[] getAcceleration(double[] acceleration
    ) {
        int n = getNumberOfVariables();
        if (acceleration == null || acceleration.length < n) {
            acceleration = new double[n];
        }
        int index = 0;
        double acc[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                a.getAcceleration(acc);
                acceleration[index++] = acc[0];
                acceleration[index++] = acc[1];
                acceleration[index++] = acc[2];
            }
        }
        return acceleration;
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration
    ) {
        int n = getNumberOfVariables();
        if (previousAcceleration == null || previousAcceleration.length < n) {
            previousAcceleration = new double[n];
        }
        int index = 0;
        double prev[] = new double[3];
        for (Atom a : refinementModel.getTotalAtomArray()) {
            if (a.isActive()) {
                a.getPreviousAcceleration(prev);
                previousAcceleration[index++] = prev[0];
                previousAcceleration[index++] = prev[1];
                previousAcceleration[index++] = prev[2];
            }
        }
        return previousAcceleration;
    }
}
