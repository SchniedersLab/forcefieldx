/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.algorithms;

import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parsers.PotentialsFunctions;
import java.io.File;

/**
 * The AlgorithmFunctions interface specifies default methods for LBFGS minimization
 * and molecular dynamics, similar to pre-existing Groovy method closures. Enables
 * FFX classes and scripts to perform these functions using either default FFX 
 * methods, a local implementation, or another implementation replacing User Interfaces.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface AlgorithmFunctions extends PotentialsFunctions {
    public void md(MolecularAssembly assembly, int nStep, double timeStep, 
            double printInterval, double saveInterval, double temperature, 
            boolean initVelocities, File dyn);
    public Potential minimize(MolecularAssembly assembly, double eps);
}
