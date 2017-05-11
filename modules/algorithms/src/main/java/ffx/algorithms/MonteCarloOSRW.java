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

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.mc.BoltzmannMC;
import ffx.algorithms.mc.LambdaMove;
import ffx.algorithms.mc.MDMove;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 *
 * @author mrtollefson
 */
public class MonteCarloOSRW extends BoltzmannMC {
    
    private final MolecularDynamics molecularDynamics;
    private final Potential potential;
    private final OSRW osrw;
    private final MolecularAssembly molecularAssembly;
    private final CompositeConfiguration properties;
    private final AlgorithmListener listener;
    private final Thermostats requestedThermostat;
    private final Integrators requestedIntegrator;
    private LambdaInterface linter;
    private double lambda = 0.0; 
    private final int mdSteps = 100;
    private final double timeStep = 1.0;
    private final double printInterval = 0.01;
    private final double temperature = 310.0;
    
    /**
     * 
     * @param molecularDynamics
     * @param potentialEnergy
     * @param osrw 
     * @param molecularAssembly 
     * @param properties 
     * @param listener 
     * @param requestedThermostat 
     * @param requestedIntegrator 
     */
    public MonteCarloOSRW(MolecularDynamics molecularDynamics, Potential potentialEnergy,
            OSRW osrw, MolecularAssembly molecularAssembly, CompositeConfiguration properties,
            AlgorithmListener listener, Thermostats requestedThermostat, Integrators requestedIntegrator) {
        this.molecularDynamics = molecularDynamics;
        this.potential = potentialEnergy;
        this.osrw = osrw;
        this.molecularAssembly = molecularAssembly;
        this.properties = properties;
        this.listener = listener;
        this.requestedThermostat = requestedThermostat;
        this.requestedIntegrator = requestedIntegrator;
    }

    public void setLambda(double lambda) {
        this.lambda = lambda;
        linter.setLambda(lambda);
    }
    
    public double getLambda() {
        return lambda;
    }
    
    /**
     * The goal is to sample coordinates (X) and converge "dU/dL" for every state (lambda)
     * along the thermodynamic path.
     * 
     * 1.) At a fixed lambda, run a defined length MD trajectory to "move" coordinates and dU/dL.
     * 2.) Accept / Reject the MD move using the OSRW energy.
     * 
     * 3.) Randomly change the value of Lambda.
     * 4.) Accept / Reject the Lambda move using the OSRW energy.
     * 
     * 5.) Add 
     */
   public void sample(){
       double [] coordinates = null; 
       double currentEnergy = potential.energy(potential.getCoordinates(coordinates));
       
      MDMove mdMove = new MDMove(molecularAssembly, potential, properties,listener, requestedThermostat, requestedIntegrator);
      mdMove.move(); //runs 100 steps of molecular dynamics to move coordinates
      
      double proposedEnergy = potential.energy(potential.getCoordinates(coordinates));
      
      if(evaluateMove(currentEnergy, proposedEnergy)){
          logger.info(String.format(" Monte Carlo step accepted: e1 -> e2 %10.6f -> %10.6f", currentEnergy, proposedEnergy));
          currentEnergy = proposedEnergy;
          //Accept MD move using OSRW energy
      }
      else{
          mdMove.revertMove();
      }
      
       
       LambdaMove lambdaMove = new LambdaMove(lambda, linter); 
       lambdaMove.move();
       double proposedLambda = linter.getLambda();
       molecularDynamics.dynamic(mdSteps, timeStep, printInterval, 10.0, temperature, true, null); //Run some MD at new lambda
       proposedEnergy = molecularDynamics.getTotalEnergy(); //Get proposed energy
       
       int[][] recursionKernel = osrw.getRecursionKernel();
       
       if(evaluateMove(currentEnergy, proposedEnergy)){
           logger.info(String.format(" Monte Carlo step accepted: e1 -> e2 %10.6f -> %10.6f", currentEnergy, proposedEnergy));
           osrw.setLambda(proposedLambda);
           //recursionKernel = recursionKernel + 1; Add one to recursionKernel matrix at current lambda bin
       }
       else{
           lambdaMove.revertMove();
           //recursionKernel = recursionKernel + 1; Add one to recursionKernel matrix at current lambda bin
       }
   }

    @Override
    protected double currentEnergy() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    protected void storeState() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void revertStep() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}

