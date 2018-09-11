/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms;

import java.io.File;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;

/**
 * The AlgorithmFunctions interface specifies default methods for L-BFGS
 * minimization and molecular dynamics, similar to pre-existing Groovy method
 * closures. Enables FFX classes and scripts to perform these functions using
 * either default FFX methods, a local implementation, or another implementation
 * replacing User Interfaces.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface AlgorithmFunctions extends PotentialsFunctions {

    /**
     * Runs molecular dynamics on an assembly.
     *
     * @param assembly       a {@link ffx.potential.MolecularAssembly} object.
     * @param nStep          Timesteps
     * @param timeStep       Time per step
     * @param printInterval  a double.
     * @param saveInterval   a double.
     * @param temperature    a double.
     * @param initVelocities Initialize velocities from Maxwell-Boltzmann distribution
     * @param dyn            Dynamics file
     */
    abstract public void md(MolecularAssembly assembly, int nStep, double timeStep,
                            double printInterval, double saveInterval, double temperature,
                            boolean initVelocities, File dyn);

    /**
     * Relax the coordinates of a MolecularAssembly and minimize its potential energy
     *
     * @param assembly a {@link ffx.potential.MolecularAssembly} object.
     * @param eps      RMS gradient convergence criteria
     * @return A <code>Potential</code>
     */
    abstract public Potential minimize(MolecularAssembly assembly, double eps);

    /**
     * Returns a default Listener if available (null by default).
     *
     * @return An AlgorithmListener or null.
     */
    default public AlgorithmListener getDefaultListener() {
        return null;
    }
}
