/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import java.util.logging.Logger;
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

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;

/**
 * OpenMM accelerated L-BFGS minimization.
 *
 * @author Hernan Bernabe
 */
public class MinimizeOpenMM {

    public static final Logger logger = Logger.getLogger(MinimizeOpenMM.class.getName());

    MolecularAssembly molecularAssembly;

    public MinimizeOpenMM(MolecularAssembly molecularAssembly) {
        this.molecularAssembly = molecularAssembly;
    }

    /**
     * <p>
     * minimize</p>
     *
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize() {
        return minimize(1.0, Integer.MAX_VALUE);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps The convergence criteria.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps) {
        return minimize(eps, Integer.MAX_VALUE);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps           The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps, int maxIterations) {
        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();

        if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
            long time = -System.nanoTime();
            ForceFieldEnergyOpenMM forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy;
            PointerByReference context = forceFieldEnergyOpenMM.getContext();

            // Run the OpenMM minimization.
            OpenMM_LocalEnergyMinimizer_minimize(context, eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ), maxIterations);

            // Get the minimized coordinates, forces and potential energy back from OpenMM.
            int infoMask = OpenMM_State_Positions + OpenMM_State_Energy + OpenMM_State_Forces;
            PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);
            double energyOpenMM = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
            PointerByReference positions = OpenMM_State_getPositions(state);
            PointerByReference forces = OpenMM_State_getForces(state);

            // Load updated coordinate position.
            int numParticles = forceFieldEnergyOpenMM.getNumParticles();
            double x[] = new double[numParticles * 3];
            forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles, x);

            // Compute the RMS gradient.
            double totalForce = 0;
            for (int i = 0; i < numParticles; i++) {
                OpenMM_Vec3 forceOpenMM = OpenMM_Vec3Array_get(forces, i);
                double fx = forceOpenMM.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double fy = forceOpenMM.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double fz = forceOpenMM.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                totalForce += fx * fx + fy * fy + fz * fz;
            }
            double grmsOpenMM = sqrt(totalForce / (3 * numParticles));

            // Clean up.
            OpenMM_State_destroy(state);

            double grad[] = new double[numParticles * 3];
            double energy = forceFieldEnergy.energyAndGradient(x, grad);
            double grms = 0.0;
            for (int i = 0; i < grad.length; i++) {
                double gi = grad[i];
                if (gi == Double.NaN
                        || gi == Double.NEGATIVE_INFINITY
                        || gi == Double.POSITIVE_INFINITY) {
                    String message = format(" The gradient of variable %d is %8.3f.", i, gi);
                    logger.warning(message);
                }
                grms += gi * gi;
            }
            grms = sqrt(grms / (3 * numParticles));

            time += System.nanoTime();
            logger.info(format(" Convergence criteria for OpenMM %12.6f vs. FFX %12.6f (kcal/mol/A).",
                    grmsOpenMM, grms));
            logger.info(format(" Final energy for         OpenMM %12.6f vs. FFX %12.6f (kcal/mol) in %8.3f (sec).",
                    energyOpenMM, energy, time * 1.0e-9));
        }

        return forceFieldEnergy;
    }

}
