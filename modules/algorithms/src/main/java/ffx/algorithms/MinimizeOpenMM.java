/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import com.sun.jna.ptr.PointerByReference;
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
import edu.uiowa.jopenmm.OpenMM_Vec3;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import static java.lang.String.format;
import java.util.logging.Level;
import java.util.logging.Logger;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
     * @param eps The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps, int maxIterations) {
        // Run the OpenMM minimization.
        long time = -System.nanoTime();

        ForceFieldEnergyOpenMM forceFieldEnergyOpenMM;

        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        int n = forceFieldEnergy.getNumberOfVariables();
        double coords[] = new double[n];
        double grad[] = new double[n];

        forceFieldEnergy.getCoordinates(coords);

        double energy = forceFieldEnergy.energyAndGradient(coords, grad);

        logger.info(format(" Minimize initial energy: %16.8f", energy));

        if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
            forceFieldEnergyOpenMM = (ForceFieldEnergyOpenMM) forceFieldEnergy;
            PointerByReference context = forceFieldEnergyOpenMM.getContext();
            OpenMM_LocalEnergyMinimizer_minimize(context, eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ), maxIterations);

            // Get the minimized coordinates and potential energy back from OpenMM. 
            int infoMask = OpenMM_State_Positions + OpenMM_State_Energy + OpenMM_State_Forces;
            PointerByReference state = OpenMM_Context_getState(context, infoMask, forceFieldEnergyOpenMM.enforcePBC);
            double currentPotentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
            PointerByReference positions = OpenMM_State_getPositions(state);
            PointerByReference forces = OpenMM_State_getForces(state);
            int numParticles = forceFieldEnergyOpenMM.getNumParticles();
            double x[] = new double[numParticles * 3];
            double forcesArr[] = new double[numParticles * 3];
            for (int i = 0; i < numParticles; i++) {
                int offset = i * 3;
                OpenMM_Vec3 forceOpenMM = OpenMM_Vec3Array_get(forces, i);
                forcesArr[offset] = forceOpenMM.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                forcesArr[offset + 1] = forceOpenMM.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                forcesArr[offset + 2] = forceOpenMM.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
            }
            
            double totalForce = 0;
            
            for (int i = 0; i < forcesArr.length; i++){
                totalForce += forcesArr[i] * forcesArr[i];
            }
            
            double grmsOpenMM = sqrt(totalForce/numParticles);
            
            forceFieldEnergyOpenMM.getOpenMMPositions(positions, numParticles, x);
            
            
            energy = forceFieldEnergy.energyAndGradient(x, grad);
            
            double rms = sqrt(numParticles);
            double grms = 0.0;
            

            for (int i = 0; i < grad.length; i++) {
                double gi = grad[i];
                if (gi == Double.NaN
                        || gi == Double.NEGATIVE_INFINITY
                        || gi == Double.POSITIVE_INFINITY) {
                    String message = format("The gradient of variable %d is %8.3f.", i, gi);
                    logger.warning(message);
                }
                double gis = gi;

                grms += gis * gis;
            }
            grms = sqrt(grms) / rms;

            time += System.nanoTime();
            logger.info(String.format(
                    " OpenMM minimization finished with energy %16.8f with GRMS %16.8f in %8.3f sec.",
                    currentPotentialEnergy, grmsOpenMM, time * 1.0e-9));
            logger.fine(String.format(" FFX calculated energy %16.8f with GRMS %16.8f", energy, grms));

            // Clean up.
            OpenMM_State_destroy(state);
        }

        return forceFieldEnergy;
    }

}
