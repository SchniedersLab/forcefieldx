/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import static ffx.algorithms.Thermostat.kB;
import ffx.potential.MolecularAssembly;
import ffx.crystal.Crystal;
import ffx.potential.OpenMMForceFieldEnergy;
import ffx.potential.extended.ExtendedSystem;
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
import static java.lang.String.format;
import static simtk.openmm.OpenMMLibrary.*;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static simtk.openmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;

/**
 * Runs Molecular Dynamics using OpenMM implementation
 *
 * @author Hernan V. Bernabe
 */
public class OpenMMMolecularDynamics extends MolecularDynamics{

    private OpenMMForceFieldEnergy openMMForceFieldEnergy;
    private static final Logger logger = Logger.getLogger(OpenMMMolecularDynamics.class.getName());
    private static final int DEFAULT_INTERVAL_STEPS = 100;
    private double temperature;
    private double currentTemperature;
    private double saveInterval;
    private double restartFrequency;
    private int intervalSteps;
    private double printInterval;
    private AlgorithmListener listener;
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
    private PointerByReference integratorReference;
    private Pointer andersenTemp;
    private Random random;
    private int randomNumber;
    private File dyn;
    private File restartFile = null;
    private String fileType = "XYZ";
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
    private double collisionFreq;

    /**
     * Constructs an OpenMMMolecularDynamics object, to perform molecular dynamics using native OpenMM routines, avoiding
     * the cost of communicating coordinates, gradients, and energies back and forth across the PCI bus.
     * @param assembly MolecularAssembly to operate on
     * @param openMMForceFieldEnergy OpenMMForceFieldEnergy Potential. Cannot be any other type of Potential.
     * @param properties Associated properties
     * @param listener
     * @param thermostat May have to be slightly modified for native OpenMM routines
     * @param integratorMD May have to be slightly modified for native OpenMM routines
     */
    public OpenMMMolecularDynamics(MolecularAssembly assembly, OpenMMForceFieldEnergy openMMForceFieldEnergy, CompositeConfiguration properties, AlgorithmListener listener, Thermostats thermostat, Integrators integratorMD) {

        super(assembly, openMMForceFieldEnergy, properties, listener, thermostat, integratorMD);
        this.openMMForceFieldEnergy = openMMForceFieldEnergy;
        openMMForceFieldEnergy.addCOMMRemover(false);
        this.molecularAssembly = assembly;
        this.listener = listener;
        assemblies = new ArrayList<>();
        assemblies.add(new AssemblyInfo(assembly));
        assemblies.get(0).props = properties;
        random = new Random();
        randomNumber = random.nextInt();
        mass = openMMForceFieldEnergy.getMass();
        numberOfVariables = openMMForceFieldEnergy.getNumberOfVariables();
        x = new double[numberOfVariables];
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        grad = new double[numberOfVariables];
    }

    /**
     * UNSUPPORTED: OpenMMMolecularDynamics is not presently capable of handling extended system variables. Will
     * throw an UnsupportedOperationException.
     * @param system
     * @param printFrequency
     */
    @Override
    public void attachExtendedSystem(ExtendedSystem system, int printFrequency) {
        throw new UnsupportedOperationException(" OpenMMMolecularDynamics does not support extended system variables!");
    }

    /**
     * UNSUPPORTED: OpenMMMolecularDynamics is not presently capable of handling extended system variables. Will
     * throw an UnsupportedOperationException.
     */
    @Override
    public void detachExtendedSystem() {
        throw new UnsupportedOperationException(" OpenMMMolecularDynamics does not support extended system variables!");
    }

    /**
     * takeSteps moves the simulation forward in time a user defined number of
     * steps and integrates the equations of motion for each step. This method
     * ensures that the algorithm reports back only when the time interval
     * (steps) specified by the user is completed.
     *
     * @param openMMForceFieldEnergy
     * @param intervalSteps
     */
    private void takeSteps(int intervalSteps) {
        double time = System.nanoTime();
        OpenMM_Integrator_step(integratorReference, intervalSteps);
        double took = (double) ((System.nanoTime() - time) * 1.0e-9);

    }

    /**
     * openMM_Update obtains the state of the simulation from OpenMM, getting
     * positions and velocities back from the OpenMM data structure.
     *
     * @param i
     */
    private void openMM_Update(int i) {
        int infoMask;
        double totalEnergy;
        double potentialEnergy;
        double kineticEnergy;

        context = openMMForceFieldEnergy.getContext();
        infoMask = OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Forces + OpenMM_State_Energy;

        PointerByReference state = OpenMM_Context_getState(context, infoMask, 0);

        potentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        kineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;

        if (i == 0) {
            velocities = OpenMM_State_getVelocities(state);
            v = openMMForceFieldEnergy.getOpenMMVelocities(velocities, numberOfVariables);

            double e = 0.0;
            for (int ii = 0; ii < v.length; ii++) {
                double velocity = v[ii];
                double v2 = velocity * velocity;
                e += mass[ii] * v2;
            }

            int dof = v.length;
            currentTemperature = e / (kB * dof);
        }

        if (i != 0) {

            positions = OpenMM_State_getPositions(state);
            x = openMMForceFieldEnergy.getOpenMMPositions(positions, numberOfVariables);

            velocities = OpenMM_State_getVelocities(state);
            v = openMMForceFieldEnergy.getOpenMMVelocities(velocities, numberOfVariables);

            double e = 0.0;
            for (int ii = 0; ii < v.length; ii++) {
                double velocity = v[ii];
                double v2 = velocity * velocity;
                e += mass[ii] * v2;
            }

            int dof = v.length;
            currentTemperature = e / (kB * dof);

            forces = OpenMM_State_getForces(state);
            a = openMMForceFieldEnergy.getOpenMMAccelerations(forces, numberOfVariables, mass);
            aPrevious = openMMForceFieldEnergy.getOpenMMAccelerations(forces, numberOfVariables, mass);

        }

        // assert energies = stuff we recalculated.
        totalEnergy = potentialEnergy + kineticEnergy;
        if (initialStep == 0) {
            logger.log(Level.INFO, " Initial energy levels");
            logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
            logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
            logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", kineticEnergy, potentialEnergy, totalEnergy, currentTemperature));
            initialStep = 1;
        } else if (i % printInterval == 0) {
            time = System.nanoTime() - time;
            logger.info(String.format(" %7d%15.4f%13.4f%13.4f%9.2f%9.3f", i, kineticEnergy, potentialEnergy,
                    totalEnergy, currentTemperature, time * 1.0e-9));
            time = System.nanoTime();
        }

        if (saveInterval > 0 && i % (saveInterval * 1000) == 0 && i != 0) {
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
        if (restartFrequency > 0 && i % (restartFrequency * 1000) == 0 && i != 0) {
            if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
            } else {
                logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
            }
        }

        OpenMM_State_destroy(state);
    }

    @Override
    public void init(int numSteps, double timeStep, double printInterval, double saveInterval, String fileType, double restartFrequency, double temperature, boolean initVelocities, File dyn){
        this.temperature = temperature;
        this.saveInterval = saveInterval;
        this.timeStep = timeStep;
        this.printInterval = (int) printInterval;
        this.dyn = dyn;
        this.initVelocities = initVelocities;
        currentTemperature = temperature;

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
        openMMForceFieldEnergy.setIntegrator(integrator, timeStep, frictionCoeff, temperature, collisionFreq);

        integratorReference = openMMForceFieldEnergy.getIntegrator();

        context = openMMForceFieldEnergy.getContext();

        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());

        if (dyn == null) {
            this.restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        if (!initialized) {

            if (loadRestart) {
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
            } else {
                openMMForceFieldEnergy.loadOpenMMPositions();
                //openMMForceFieldEnergy.setOpenMMPositions();

                if (initVelocities) {

                    OpenMM_Context_setVelocitiesToTemperature(context, temperature, randomNumber);
                }
            }

        }
    }
    
    /**
     * start sets up context, write out file name, restart file name, sets the
     * integrator and determines whether the simulation is starting out from a
     * previous molecular dynamics run (.dyn) or if the initial velocities are
     * determined by a Maxwell Boltzmann distribution. This method then calls
     * methods openMMUPdate and takeSteps to run the molecular dynamics
     * simulation.
     *
     * @param numSteps
     * @param intervalSteps
     */
    @Override
    public void dynamic(int numSteps, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {
        int i = 0;
        
        init(numSteps, timeStep, printInterval, saveInterval, fileType, temperature, restartFrequency, initVelocities, dyn);
        /*
        this.temperature = temperature;
        this.saveInterval = saveInterval;
        this.timeStep = timeStep;
        this.dyn = dyn;
        this.initVelocities = initVelocities;
        currentTemperature = temperature;
        this.printInterval = (int) printInterval;

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
        openMMForceFieldEnergy.setIntegrator(integrator, timeStep, frictionCoeff, temperature, collisionFreq);

        integratorReference = openMMForceFieldEnergy.getIntegrator();

        context = openMMForceFieldEnergy.getContext();

        String firstFileName = FilenameUtils.removeExtension(molecularAssembly.getFile().getAbsolutePath());

        if (dyn == null) {
            this.restartFile = new File(firstFileName + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        if (!initialized) {

            if (loadRestart) {
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
            } else {
                openMMForceFieldEnergy.loadOpenMMPositions();
                //openMMForceFieldEnergy.setOpenMMPositions();

                if (initVelocities) {

                    OpenMM_Context_setVelocitiesToTemperature(context, temperature, randomNumber);
                }
            }

        }
        */
        time = System.nanoTime();
        while (i < numSteps) {
            openMM_Update(i);
            takeSteps(intervalSteps);
            i = i + intervalSteps;
            finalSteps = i;
        }

        openMM_Update(finalSteps);
    }

    /**
     * A simple container class to hold all the infrastructure associated with a
     * MolecularAssembly for MolecularDynamics; assembly, properties, archive
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

    public void setLambda(double newLambda) {
        lambda = newLambda;
    }

    /**
     * Method to set file type from groovy scripts.
     *
     * @param fileType the type of snapshot files to write.
     */
    public void setFileType(String fileType) {
        this.fileType = fileType;
    }

    /**
     * Method to set the Restart Frequency.
     *
     * @param restartFrequency the time between writing restart files.
     */
    public void setRestartFrequency(double restartFrequency) {
        this.restartFrequency = restartFrequency;
    }

    /**
     * Method to set the Integrator string
     * 
     * @param integrator string name of the integrator to be used during the simulation
     */
    public void setIntegratorString(String integrator){
        this.integrator = integrator;
    }
    
    /**
     * Method to set the coefficient of friction (for Langevin integrator)
     * 
     * @param frictionCoeff coefficient of friction that acts on the atoms during simulation
     */
    public void setFrictionCoefficient(double frictionCoeff){
        this.frictionCoeff = frictionCoeff;
    }
    
    /**
     * Method to set the collision frequency (for Andersen thermostat)
     * 
     * @param collisionFreq rate at which atoms collide in the Andersen thermostat
     */
    public void setCollisionFrequency(double collisionFreq){
        this.collisionFreq = collisionFreq;
    }
    
    public void setIntervalSteps(int intervalSteps){
        this.intervalSteps = intervalSteps;
    }
}
