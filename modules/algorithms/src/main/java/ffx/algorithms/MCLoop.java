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
package ffx.algorithms;

import java.util.List;
import java.util.Random;
import java.util.logging.Logger;
import java.util.concurrent.ThreadLocalRandom;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Loop;
import java.util.ArrayList;

/**
 * @author Armin Avdic
 */
public class MCLoop implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(MCLoop.class.getName());
    /**
     * The MoleularAssembly.
     */
    private final MolecularAssembly molAss;
    /**
     * The MolecularDynamics object controlling the simulation.
     */
    private MolecularDynamics molDyn;
    /**
     * The MD thermostat.
     */
    private final Thermostat thermostat;
    /**
     * Boltzmann's constant is kcal/mol/Kelvin.
     */
    private static final double boltzmann = 0.0019872041;
    /**
     * Energy of the system at initialization.
     */
    private final double systemReferenceEnergy;
    /**
     * The current MD step.
     */
    private int stepCount = 0;
    /**
     * Number of simulation steps between MC move attempts.
     */
    private int mcStepFrequency;
    /**
     * Number of accepted MD moves.
     */
    private int numMovesAccepted;
    /**
     * Storage for coordinates before MC move.
     */
    private double[] oldCoords;
    /**
     * KIC generation of loop solutions. See doi:10.1002/jcc.10416
     */    
    Loop loop;
    /**
     * Number of KIC iterations per MC move.
     */ 
    private int iterations;
    /**
     * List of active atoms.
     */ 
    private Atom[] atoms;
     /**
     * First and last residue numbers of loop.
     */ 
    private final int firstResidue;
    private final int endResidue;
    /**
     * Everyone's favorite.
     */
    private final Random rng = new Random();
    /**
     * The ForceFieldEnergy object being used by MD.
     */
    private final ForceFieldEnergy forceFieldEnergy;
     /**
     * The LambdaInterface object being used by OSRW.
     */
    private LambdaInterface lambdaInterface;
    private boolean skipAlgorithm = false;

    /**
     * Construct a Monte-Carlo loop switching mechanism.
     *
     * @param molAss the molecular assembly
     * @param mcStepFrequency number of MD steps between switch attempts
     * @param thermostat the MD thermostat
     */
    MCLoop(MolecularAssembly molAss, int mcStepFrequency, Thermostat thermostat, int firstResidue, int endResidue) {
        numMovesAccepted = 0;

        this.molAss = molAss;
        this.atoms = molAss.getAtomArray();
        this.forceFieldEnergy = molAss.getPotentialEnergy();
        this.mcStepFrequency = (mcStepFrequency == 0) ? Integer.MAX_VALUE : mcStepFrequency;
        this.thermostat = thermostat;
        systemReferenceEnergy = molAss.getPotentialEnergy().energy(false, true);
        this.firstResidue = firstResidue;
        this.endResidue = endResidue; 
        this.iterations = 1;

        if ((endResidue - firstResidue) < 3){
            logger.info("MCLoop requires at least 3 residues. First and last residues are anchors.");
            skipAlgorithm = true;
        }
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Running MCLoop:\n"));
        sb.append(String.format("     mcStepFrequency: %4d\n", mcStepFrequency));
        sb.append(String.format("     referenceEnergy: %7.2f\n", systemReferenceEnergy));
        logger.info(sb.toString());

        loop = new Loop(molAss);
    }

    /**
     * Identify Random Portion of 3 residues from Loop and run KIC.
     */
   // private  List <double[]> generateLoops() {
        //TODO: When general KIC is implemented
        //int subLength = rng.nextInt ( end_res - stt_res - 1);
        //int subLength = 3; 
        //int loopLength = end_res - stt_res;
        //int offset = rng.nextInt( loopLength );
      //  return loop.getSolutions().size();
   //     return loop.generateLoops(stt_res, end_res);
  //  }

    /**
     * The primary driver. Called by the MD engine at each dynamics step.
     * @param molAss
     */
    @Override
    public boolean mcUpdate(MolecularAssembly molAss) {

        stepCount++;
        if (skipAlgorithm == true){
            return false;
        }
        // Decide on the type of step to be taken.
        if ((stepCount % mcStepFrequency != 0 )) {
            // Not yet time for an MC step, return to MD.
            return false;
        }
        
        if(lambdaInterface != null){
            if (lambdaInterface.getLambda() > 0.1){
                logger.info(String.format(" KIC procedure skipped (Lambda > 0.1)."));
                return false;
            }
        }
        
        atoms = molAss.getAtomArray();

        // Randomly choose a target sub portion of loop to KIC.
        int midResidue;
        midResidue = ThreadLocalRandom.current().nextInt(firstResidue + 1, endResidue);
   
        List <double[]> loopSolutions;     
        loopSolutions = loop.generateLoops(midResidue - 1, midResidue + 1);

        for(int i = 1; i < iterations; i++){
            //pick random subloop
            midResidue = ThreadLocalRandom.current().nextInt(firstResidue + 1, endResidue);
            //pick random solution
            if (loopSolutions.size() > 0){
                List <double[]> tempLoops = loop.generateLoops(midResidue - 1, midResidue + 1,loopSolutions.get(rng.nextInt(loopSolutions.size())));
                for (double[] tempLoop : tempLoops) {
                    loopSolutions.add(tempLoop);
                }
            } else {
                loopSolutions = loop.generateLoops(midResidue - 1, midResidue + 1);
            }
                    
        }
        int numLoopsFound = loopSolutions.size();
        // Check whether KIC found alternative loops
        if (numLoopsFound <= 1) {
            return false;
        }     
        
        // Perform the MC move.
        boolean accepted = tryLoopStep(loopSolutions);
        return accepted;                 
    }

    /**
     * Perform a loop MC move.
     * @param loopSolutions
     * @return accept/reject
     */
    private boolean tryLoopStep(List <double[]> loopSolutions) {

        // Choose from the list of available loops and save current coordinates
        double[] newCoords = loopSolutions.get(rng.nextInt(loopSolutions.size()));       
        oldCoords = storeActiveCoordinates();

        // Optimize the system.
        Minimize minimize1 = new Minimize(null, forceFieldEnergy, null);
        minimize1.minimize();
        double originalLoopEnergy = forceFieldEnergy.energy(false,true);
        
        // Perform move and analysis of chosen loop.
        performLoopMove(newCoords);
        Minimize minimize2 = new Minimize(null, forceFieldEnergy, null);
        minimize2.minimize();
        double newLoopEnergy = forceFieldEnergy.energy(false,true);
        
        // Remove the scaling of coordinates & gradient set by the minimizer.
        forceFieldEnergy.setScaling(null);
                
        double temperature = thermostat.getCurrentTemperature();
        double kT = boltzmann * temperature;
        // Test the MC criterion for a loop move.
        double dE = Math.abs(originalLoopEnergy - newLoopEnergy);
        if (newLoopEnergy < originalLoopEnergy){
            dE = -dE;
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format(" Assessing possible MC loop move step:\n"));
        sb.append(String.format("     original loop: %16.8f\n", originalLoopEnergy));
        sb.append(String.format("     possible loop: %16.8f\n", newLoopEnergy));
        sb.append(String.format("     -----\n"));
        
        // Test Monte-Carlo criterion.
        if (dE < 0 ) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        double criterion = exp(-dE / kT);

        double metropolis = random();
        sb.append(String.format("     criterion:  %9.4f\n", criterion));
        sb.append(String.format("     rng:        %9.4f\n", metropolis));
        if ((metropolis < criterion )) {
            sb.append(String.format("     Accepted!"));
            logger.info(sb.toString());
            numMovesAccepted++;
            return true;
        }
        sb.append(String.format("     Denied."));
        logger.info(sb.toString());

        // Move was rejected, undo the loop move
        performLoopMove(oldCoords);
        return false;
    }

    /**
     * Perform the requested coordinate move
     *
     * @param newCoordinates
     */
    private void performLoopMove(double[] newCoordinates) {
        int index = 0;
        for (Atom a : atoms) {
            double x = newCoordinates[index++];
            double y = newCoordinates[index++];
            double z = newCoordinates[index++];
            a.moveTo(x, y, z);
        }
    }

    /**
     * Get the current MC acceptance rate.
     *
     * @return the acceptance rate.
     */
    public double getAcceptanceRate() {
        // Intentional integer division.
        int numTries = stepCount / mcStepFrequency;
        return (double) numMovesAccepted / numTries;
    }

    public void setIterations(int iterations){
        this.iterations = iterations;
    }
    public void addMolDyn(MolecularDynamics molDyn) {
        this.molDyn = molDyn;
    }
    
    public void addLambdaInterface(LambdaInterface lambdaInterface) {
        this.lambdaInterface = lambdaInterface;
    }
    
    private double[] storeActiveCoordinates() {
        double [] coords = new double[atoms.length*3];
        int index = 0;
        for (Atom a : atoms) {
            coords[index++] = a.getX();
            coords[index++] = a.getY();
            coords[index++] = a.getZ();
        }
        return coords;
    }

    private enum MCOverride {
        ACCEPT, REJECT, NONE;
    }
}
