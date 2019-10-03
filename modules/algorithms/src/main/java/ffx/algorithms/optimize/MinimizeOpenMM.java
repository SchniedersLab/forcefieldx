//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.algorithms.optimize;

import java.util.logging.Logger;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

import com.sun.jna.ptr.PointerByReference;

import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.uiowa.jopenmm.OpenMM_Vec3;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalEnergyMinimizer_minimize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_get;

import ffx.algorithms.AlgorithmListener;
import ffx.numerics.Potential;
import ffx.numerics.optimization.LineSearch;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

/**
 * OpenMM accelerated L-BFGS minimization.
 *
 * @author Hernan Bernabe
 */
public class MinimizeOpenMM extends Minimize {

    private static final Logger logger = Logger.getLogger(MinimizeOpenMM.class.getName());

    public MinimizeOpenMM(MolecularAssembly molecularAssembly) {
        super(molecularAssembly, molecularAssembly.getPotentialEnergy(), null);
    }

    public MinimizeOpenMM(MolecularAssembly molecularAssembly, ForceFieldEnergyOpenMM forceFieldEnergyOpenMM) {
        super(molecularAssembly, forceFieldEnergyOpenMM, null);
    }

    public MinimizeOpenMM(MolecularAssembly molecularAssembly, ForceFieldEnergyOpenMM forceFieldEnergyOpenMM,
                          AlgorithmListener algorithmListener) {
        super(molecularAssembly, forceFieldEnergyOpenMM, algorithmListener);
    }


    /**
     * Note the OpenMM L-BFGS minimizer does not accept the parameter "m"
     * for the number of previous steps used to estimate the Hessian.
     *
     * @param m             The number of previous steps used to estimate the Hessian (ignored).
     * @param eps           The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     * @return The potential.
     */
    @Override
    public Potential minimize(int m, double eps, int maxIterations) {
        return minimize(eps, maxIterations);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps           The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     * @return a {@link ffx.numerics.Potential} object.
     */
    @Override
    public Potential minimize(double eps, int maxIterations) {

        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();

        if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
            time = -System.nanoTime();
            ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy;

            // Get the OpenMM Context.
            PointerByReference context = forceFieldEnergyOpenMM.getContext();

            // Respect the use flag, and lambda state.
            Atom[] atoms = molecularAssembly.getAtomArray();
            forceFieldEnergyOpenMM.updateParameters(atoms);

            // Respect (in)active atoms.
            forceFieldEnergyOpenMM.setActiveAtoms();

            double[] x2 = x;
            int nAtoms = atoms.length;
            if (x.length != nAtoms * 3) {
                // In FFX, the x coordinate array only includes active atoms, but OpenMM operates on all atoms.
                x2 = new double[nAtoms * 3];
                int index = 0;
                for (Atom atom : atoms) {
                    x2[index++] = atom.getX();
                    x2[index++] = atom.getY();
                    x2[index++] = atom.getZ();
                }
            } else {
                forceFieldEnergyOpenMM.getCoordinates(x);
            }

            forceFieldEnergyOpenMM.setOpenMMPositions(x2, x2.length);

            // Calculate the starting energy before optimization.
            double e = forceFieldEnergyOpenMM.energy(x2);
            logger.info(format(" Initial energy:                 %12.6f (kcal/mol)", e));

            // Run the minimization.
            OpenMM_LocalEnergyMinimizer_minimize(context, eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ), maxIterations);

            // Get the minimized coordinates, forces and potential energy back from OpenMM.
            int infoMask = OpenMM_State_Positions + OpenMM_State_Energy + OpenMM_State_Forces;
            PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);
            energy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
            PointerByReference positions = OpenMM_State_getPositions(state);
            PointerByReference forces = OpenMM_State_getForces(state);

            // Load updated coordinate position.
            forceFieldEnergyOpenMM.getOpenMMPositions(positions, nAtoms, x2);

            // Compute the RMS gradient.
            int index = 0;
            double totalForce = 0;
            for (int i = 0; i < nAtoms; i++) {
                OpenMM_Vec3 forceOpenMM = OpenMM_Vec3Array_get(forces, i);
                double fx = forceOpenMM.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double fy = forceOpenMM.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double fz = forceOpenMM.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                if (atoms[i].isActive()) {
                    grad[index++] = -fx;
                    grad[index++] = -fy;
                    grad[index++] = -fz;
                    totalForce += fx * fx + fy * fy + fz * fz;
                    if (isNaN(totalForce) || isInfinite(totalForce)) {
                        String message = format(" The gradient of variable %d is %8.3f.", i, totalForce);
                        logger.warning(message);
                    }
                }
            }
            rmsGradient = sqrt(totalForce / n);

            // Clean up.
            OpenMM_State_destroy(state);

            double[] ffxGrad = new double[n];
            forceFieldEnergy.getCoordinates(x);
            double ffxEnergy = forceFieldEnergyOpenMM.energyAndGradientFFX(x, ffxGrad);
            double grmsFFX = 0.0;
            for (int i = 0; i < n; i++) {
                double gi = ffxGrad[i];
                if (isNaN(gi) || isInfinite(gi)) {
                    String message = format(" The gradient of variable %d is %8.3f.", i, gi);
                    logger.warning(message);
                }
                grmsFFX += gi * gi;
            }
            grmsFFX = sqrt(grmsFFX / n);

            time += System.nanoTime();
            logger.info(format(" Final energy for OpenMM         %12.6f vs. FFX %12.6f in %8.3f (sec).", energy, ffxEnergy, time * 1.0e-9));
            logger.info(format(" Convergence criteria for OpenMM %12.6f vs. FFX %12.6f (kcal/mol/A).", rmsGradient, grmsFFX));

        }

        if (algorithmListener != null) {
            algorithmListener.algorithmUpdate(molecularAssembly);
        }

        return forceFieldEnergy;
    }

    /**
     * MinimizeOpenMM does not currently support the OptimizationListener interface.
     *
     * @since 1.0
     */
    @Override
    public boolean optimizationUpdate(int iteration, int functionEvaluations, double rmsGradient,
                                      double rmsCoordinateChange, double energy, double energyChange,
                                      double angle, LineSearch.LineSearchResult lineSearchResult) {
        logger.warning(" MinimizeOpenMM does not support updates at each optimization step.");
        return false;
    }


}
