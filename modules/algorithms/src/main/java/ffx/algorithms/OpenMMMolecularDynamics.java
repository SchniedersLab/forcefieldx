/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import com.sun.jna.ptr.PointerByReference;
import ffx.potential.MolecularAssembly;
import ffx.crystal.Crystal;
import ffx.potential.OpenMMForceFieldEnergy;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import java.util.Random;
import ffx.potential.parsers.DYNFilter;


//import static simtk.openmm.OpenMMAmoebaLibrary.*;
import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;


/**
 *
 * @author hbernabe
 */
public class OpenMMMolecularDynamics {
    private OpenMMForceFieldEnergy openMMForceFieldEnergy;
    private static final Logger logger = Logger.getLogger(OpenMMMolecularDynamics.class.getName());
    private double temperature;
    private double saveInterval;
    private double restartFrequency;
    private final List<AssemblyInfo> assemblies;
    private final MolecularAssembly molecularAssembly;
    private boolean saveSnapshotAsPDB = true;
    private long time;
    private int initialStep = 0;
    private int finalSteps;
    private String integrator;
    private double timeStep;
    private double frictionCoeff;
    private double lambda;
    private PointerByReference context;
    private Random random;
    private int randomNumber;
    private File dyn;
    private File restartFile = null;
    private boolean loadRestart = false;
    private boolean initialized = false;
    private int numberOfVariables;
    private boolean done = true;
    private boolean initVelocities;
    private double[] x;
    private double[] v;
    private double[] a;
    private double[] aPrevious;
    private double[] grad;
    private double[] mass;
    private DYNFilter dynFilter = null;
    private PointerByReference positions;
    private PointerByReference velocities;
    private PointerByReference forces;
    private String positionsString;
    
    
    public OpenMMMolecularDynamics(MolecularAssembly assembly, OpenMMForceFieldEnergy openMMForceFieldEnergy, CompositeConfiguration properties, double temperature, 
            double saveInterval, double restartFrequency, String integrator, double timeStep, double frictionCoeff, File dyn, boolean initVelocities){
        
        this.openMMForceFieldEnergy = openMMForceFieldEnergy;
        this.temperature = temperature;
        this.saveInterval = saveInterval;
        this.restartFrequency = restartFrequency;
        this.molecularAssembly = assembly;
        this.integrator = integrator;
        this.timeStep = timeStep;
        this.frictionCoeff = frictionCoeff;
        this.dyn = dyn;
        this.initVelocities = initVelocities;
        assemblies = new ArrayList<>();
        assemblies.add(new AssemblyInfo(assembly));
        
        assemblies.get(0).props = properties;
        random = new Random();
        randomNumber = random.nextInt();
        mass = openMMForceFieldEnergy.getMass();
        numberOfVariables = openMMForceFieldEnergy.getNumberOfVariables();
        x = new double [numberOfVariables];
        v = new double [numberOfVariables];
        a = new double [numberOfVariables];
        aPrevious = new double [numberOfVariables];
        grad = new double [numberOfVariables];
    }
    
    /**
     * takeSteps moves the simulation forward in time a user defined number of steps
     * and integrates the equations of motion for each step. This method ensures that
     * the algorithm reports back only when the time interval (steps) specified by the
     * user is completed.
     * @param openMMForceFieldEnergy
     * @param numSteps 
     */
    
    public void takeSteps(OpenMMForceFieldEnergy openMMForceFieldEnergy, int numSteps){
        PointerByReference integrator = openMMForceFieldEnergy.getIntegrator();
        double time = System.nanoTime();
        OpenMM_Integrator_step(integrator, numSteps);
        double took = (double) ((System.nanoTime() - time)* 1.0e-9);
        
    }
    
    /**
     * openMM_Update obtains the state of the simulation from OpenMM, getting
     * positions and velocities back from the OpenMM data structure.
     * 
     * @param openMMForceFieldEnergy
     * @param i
     */
    public void openMM_Update(OpenMMForceFieldEnergy openMMForceFieldEnergy, int i){
        int infoMask;
        //double aMass;
        //double positionConvert;
        //double velocityConvert;
        //double forceConvert;
        double totalEnergy;
        double potentialEnergy;
        double kineticEnergy;
        
        
        context = openMMForceFieldEnergy.getContext();
        infoMask = OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Forces + OpenMM_State_Energy;
        
        PointerByReference state = OpenMM_Context_getState(context, infoMask, 0);
        
        potentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        kineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        
        // assert energies = stuff we recalculated.
        totalEnergy = potentialEnergy + kineticEnergy;
        if (initialStep == 0){
            logger.log(Level.INFO, " Initial energy levels");
            logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
            logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
            logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", kineticEnergy, potentialEnergy, totalEnergy, temperature));
            initialStep = 1;
        }
        else if (i % 100 == 0){
            time = System.nanoTime() - time;
            logger.info(String.format(" %7d%15.4f%13.4f%13.4f%9.2f%9.3f", i, kineticEnergy, potentialEnergy,
                        totalEnergy, temperature, time * 1.0e-9));
            time = System.nanoTime();
        }
        
        if(i != 0){
            openMMForceFieldEnergy.updateOpenMMPositions(state); // maybe also velocities, accelerations?
        }
        
        if (saveInterval > 0 && i % saveInterval == 0 && i != 0){
            for (AssemblyInfo ai : assemblies) {
                    if (ai.archiveFile != null && !saveSnapshotAsPDB) {
                        if (ai.xyzFilter.writeFile(ai.archiveFile, true)) {
                            logger.info(String.format(" Appended snap shot to %s", ai.archiveFile.getName()));
                        } else {
                            logger.warning(String.format(" Appending snap shot to %s failed", ai.archiveFile.getName()));
                        }
                    } else if (saveSnapshotAsPDB) {
                        if (ai.pdbFilter.writeFile(ai.pdbFile, false)) {
                            logger.info(String.format(" Wrote PDB file to %s", ai.pdbFile.getName()));
                        } else {
                            logger.warning(String.format(" Writing PDB file to %s failed.", ai.pdbFile.getName()));
                        }
                    }
                }
        }
        
        if (restartFrequency > 0 && i % restartFrequency == 0 && i != 0){
            /*
            positions = OpenMM_State_getPositions(state);
        
            x = positions.getPointer().getDoubleArray(0, numberOfVariables);
        
            velocities = OpenMM_State_getVelocities(state);
            v = velocities.getPointer().getDoubleArray(0, numberOfVariables);
        
            forces = OpenMM_State_getForces(state);
            a = forces.getPointer().getDoubleArray(0, numberOfVariables);
        
            for (int j = 0; j < numberOfVariables; j ++){
                a[j] = a[j]/mass[j];
                aPrevious[j] = a[j];
            }
            */
            positions = OpenMM_State_getPositions(state);
            x = openMMForceFieldEnergy.getOpenMMPositions(positions, numberOfVariables);
            
            velocities = OpenMM_State_getVelocities(state);
            v = openMMForceFieldEnergy.getOpenMMVelocities(velocities, numberOfVariables);
            
            forces = OpenMM_State_getForces(state);
            a = openMMForceFieldEnergy.getOpenMMAccelerations(forces, numberOfVariables, mass);
            aPrevious = openMMForceFieldEnergy.getOpenMMAccelerations(forces, numberOfVariables, mass);
        
        
            /*
            positions = OpenMM_State_getPositions(state);
            velocities = OpenMM_State_getVelocities(state);
            forces = OpenMM_State_getForces(state);
            for (int j = 0; j < numberOfVariables; j++){
            x[j] = positions.getPointer().getDouble(j);
            }
        
            for (int j = 0; j < numberOfVariables; j++){
                v[j] = velocities.getPointer().getDouble(j);
            }
        
            for (int j = 0; j < numberOfVariables; j++){
                double temp = forces.getPointer().getDouble(j);
                a[j] = temp/mass[j];
                aPrevious[j] = a[j];
            }
            */
        
            }
        
          /**
             * Write out restart files every saveRestartFileFrequency steps.
             */
            if (restartFrequency > 0 && i % restartFrequency == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }
            
    }
    
    /**
     * start sets up context, write out file name, restart file name, sets the integrator and determines
     * whether the simulation is starting out from a previous molecular dynamics run (.dyn) or if the
     * initial velocities are determined by a Maxwell Boltzmann distribution. This method then calls methods
     * openMMUPdate and takeSteps to run the molecular dynamics simulation.
     * @param numSteps
     * @param intervalSteps 
     */
    
    public void start(int numSteps, int intervalSteps){
        int i = 0;
        done = false;
        assemblies.stream().parallel().forEach((ainfo) -> {
            MolecularAssembly mola = ainfo.getAssembly();
            CompositeConfiguration aprops = ainfo.props;
            File file = mola.getFile();
            String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
            File archFile = ainfo.archiveFile;
            if (archFile == null) {
                archFile = new File(filename + ".arc");
                ainfo.archiveFile = XYZFilter.version(archFile);
            }
            if (ainfo.pdbFile == null) {
                String extName = FilenameUtils.getExtension(file.getName());
                if (extName.toLowerCase().startsWith("pdb")) {
                    ainfo.pdbFile = file;
                } else {
                    ainfo.pdbFile = new File(filename + ".pdb");
                }
            }
            if (ainfo.xyzFilter == null) {
                ainfo.xyzFilter = new XYZFilter(file, mola, mola.getForceField(), aprops);
            }
            if (ainfo.pdbFilter == null) {
                ainfo.pdbFilter = new PDBFilter(ainfo.pdbFile, mola, mola.getForceField(), aprops);
            }
        });
        //openMMForceFieldEnergy.setIntegrator(integrator, timeStep, frictionCoeff, temperature);
        
        context = openMMForceFieldEnergy.getContext();
        
        
        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());
        
        if (dyn == null){
            this.restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }
        
        
        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }
        
        if (!initialized){
            
            if (loadRestart){
                Crystal crystal = molecularAssembly.getCrystal();
                // possibly add check to see if OpenMM supports this space group.
                if (!dynFilter.readDYN(restartFile, crystal, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    return;
                } else {
                    molecularAssembly.getPotentialEnergy().setCrystal(crystal);
                    openMMForceFieldEnergy.setOpenMMPositions(x, numberOfVariables);
                    openMMForceFieldEnergy.setOpenMMVelocities(v, numberOfVariables);
                }
            }
            
            
            else {
                openMMForceFieldEnergy.loadOpenMMPositions();
                //openMMForceFieldEnergy.setOpenMMPositions();
                
                if (initVelocities){
                    
                    OpenMM_Context_setVelocitiesToTemperature(context, temperature, randomNumber);
                }
            }
            
        }
        logger.info(String.format("Out of position and velocity assignment"));

        time = System.nanoTime();
        while (i < numSteps){
            openMM_Update(openMMForceFieldEnergy, i);
            takeSteps(openMMForceFieldEnergy, intervalSteps);
            i = i + intervalSteps;
            finalSteps = i;
        }
        
        openMM_Update(openMMForceFieldEnergy, finalSteps);     
    }
    
    /**
     * A simple container class to hold all the infrastructure associated with
     * a MolecularAssembly for MolecularDynamics; assembly, properties, archive
     * and PDB files, PDB and XYZ filters. Direct access to package-private
     * members breaks encapsulation a bit, but the private inner class shouldn't
     * be used externally anyways.
     */
    private class AssemblyInfo {
        private final MolecularAssembly mola;
        CompositeConfiguration props = null;
        File archiveFile = null;
        File pdbFile = null;
        PDBFilter pdbFilter = null;
        XYZFilter xyzFilter = null;

        public AssemblyInfo(MolecularAssembly assembly) {
            this.mola = assembly;
        }

        public MolecularAssembly getAssembly() {
            return mola;
        }
    }
    
    public void setLambda(double newLambda){
        lambda = newLambda;
    }
    
    
    
}
