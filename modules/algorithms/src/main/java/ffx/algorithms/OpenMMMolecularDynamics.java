/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import com.sun.jna.ptr.PointerByReference;
import ffx.potential.MolecularAssembly;
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
    private File restartFile = null;
    private int numberOfVariables;
    private double[] x;
    private double[] v;
    private double[] a;
    private double[] aPrevious;
    private double[] grad;
    private double[] mass;
    private DYNFilter dynFilter = null;
    
    
    public OpenMMMolecularDynamics(MolecularAssembly assembly, OpenMMForceFieldEnergy openMMForceFieldEnergy, CompositeConfiguration properties, double temperature, 
            double saveInterval, double restartFrequency, String integrator, double timeStep, double frictionCoeff){
        
        this.openMMForceFieldEnergy = openMMForceFieldEnergy;
        this.temperature = temperature;
        this.saveInterval = saveInterval;
        this.restartFrequency = restartFrequency;
        this.molecularAssembly = assembly;
        this.integrator = integrator;
        this.timeStep = timeStep;
        this.frictionCoeff = frictionCoeff;
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
    
    public void takeSteps(OpenMMForceFieldEnergy openMMForceFieldEnergy, int numSteps){
        PointerByReference integrator = openMMForceFieldEnergy.getIntegrator();
        double time = System.nanoTime();
        OpenMM_Integrator_step(integrator, numSteps);
        double took = (double) ((System.nanoTime() - time)* 1.0e-9);
        //logger.log(Level.INFO, "Total time elapsed : {0} (sec) per 100 steps", took);
        //logger.info(String.format(" Total time elapsed : %1.3f (sec) per 100 steps", took));
        
    }
    
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
        
        //positionConvert = 1.0/OpenMM_NmPerAngstrom;
        //velocityConvert = 1.0/OpenMM_NmPerAngstrom;
        //forceConvert = 10.0;
        
        //PointerByReference positionArray = OpenMM_Vec3Array_create();
        
        potentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        kineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        totalEnergy = potentialEnergy + kineticEnergy;
        
        if (initialStep == 0){
            logger.log(Level.INFO, " Initial energy levels");
            logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
            logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
            logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", kineticEnergy, potentialEnergy, totalEnergy, temperature));
            initialStep = 1;
            //logger.log(Level.INFO, "At step {0} energy levels are:", i);
            //logger.log(Level.INFO, "OpenMM Potential Energy: {0}", potentialEnergy);
            //logger.log(Level.INFO, "OpenMM Kinetic Energy: {0}", kineticEnergy);
            //logger.log(Level.INFO, "OpenMM Total Energy: {0}", totalEnergy);
        }
        else if (i % 100 == 0){
            time = System.nanoTime() - time;
            logger.info(String.format(" %7d%15.4f%13.4f%13.4f%9.2f%9.3f", i, kineticEnergy, potentialEnergy,
                        totalEnergy, temperature, time * 1.0e-9));
            time = System.nanoTime();
        }
        
        if(i != 0){
            openMMForceFieldEnergy.updateOpenMMPositions(state);
        }
        /*
        openMMForceFieldEnergy.getCoordinates(x);
        openMMForceFieldEnergy.setVelocity(v);
        openMMForceFieldEnergy.setAcceleration(a);
        openMMForceFieldEnergy.setPreviousAcceleration(aPrevious);
        */
        if (saveInterval > 0 && i % saveInterval == 0 && i != 0){
            //logger.info(String.format(" Wrote PDB file to %s", openMMForceFieldEnergy.molecularAssembly.pdbFile.getName()));
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
        
          /**
             * Write out restart files every saveRestartFileFrequency steps.
             */
            /*if (restartFrequency > 0 && i % restartFrequency == 0 && i != 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }*/
    }
    
    public void start(int numSteps, int intervalSteps){
        int i = 0;
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
        
        OpenMM_Context_setVelocitiesToTemperature(context, temperature, randomNumber);
        
        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());
        
        this.restartFile = new File(firstFileName + ".dyn");

        //openMM_Update(openMMForceFieldEnergy, i);
        //int finalSteps = 0;
        //long time = System.nanoTime();
        time = System.nanoTime();
        while (i < numSteps){
            openMM_Update(openMMForceFieldEnergy, i);
            takeSteps(openMMForceFieldEnergy, intervalSteps);
            i = i + intervalSteps;
            finalSteps = i;
        }
        
        openMM_Update(openMMForceFieldEnergy, finalSteps);
        //long took = (long) ((System.nanoTime() - time) * 1e-9);
        //logger.log(Level.INFO, "Total time elapsed : {0} (sec)", took);
        //logger.log(Level.INFO, "Final energy at step {0}:", finalSteps);
        
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
