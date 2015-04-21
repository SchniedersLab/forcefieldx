/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

import java.io.File;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.fill;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.parsers.DYNFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.XYZFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Run NVE or NVT molecular dynamics.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class MolecularDynamics implements Runnable, Terminatable {

    private static final Logger logger = Logger.getLogger(MolecularDynamics.class.getName());
    private final MolecularAssembly molecularAssembly;
    private final Potential potential;
    private final CompositeConfiguration properties;
    private AlgorithmListener algorithmListener;
    private MonteCarloListener monteCarloListener;
    private Thermostat thermostat;
    private Integrator integrator;
    private File archiveFile = null;
    private File restartFile = null;
    private File pdbFile = null;
    private XYZFilter xyzFilter = null;
    private DYNFilter dynFilter = null;
    private PDBFilter pdbFilter = null;
    private int printFrequency = 100;
    private int saveSnapshotFrequency = 1000;
    private int removeCOMMotionFrequency = 100;
    private boolean initVelocities = true;
    private boolean loadRestart = false;
    private boolean initialized = false;
    private boolean done = true;
    private boolean terminate = false;
//    private final int numberOfVariables;
//    private final double[] x;
//    private final double[] v;
//    private final double[] a;
//    private final double[] aPrevious;
//    private final double[] grad;
//    private final double[] mass;
    private int numberOfVariables;
    private double[] x;
    private double[] v;
    private double[] a;
    private double[] aPrevious;
    private double[] grad;
    private double[] mass;
    private int nSteps = 1000;
    private double targetTemperature = 298.15;
    private double dt = 1.0;
    private double currentTemperature;
    private double currentKineticEnergy;
    private double currentPotentialEnergy;
    private double currentTotalEnergy;
    private boolean saveSnapshotAsPDB = true;
    private int saveRestartFileFrequency = 1000;
    private String fileType = "PDB";
    private double restartFrequency = 0.1;

    /**
     * <p>
     * Constructor for MolecularDynamics.</p>
     *
     * @param assembly a {@link ffx.potential.MolecularAssembly} object.
     * @param potentialEnergy a {@link ffx.numerics.Potential} object.
     * @param properties a
     * {@link org.apache.commons.configuration.CompositeConfiguration} object.
     * @param listener a {@link ffx.algorithms.AlgorithmListener} object.
     * @param requestedThermostat a
     * {@link ffx.algorithms.Thermostat.Thermostats} object.
     * @param requestedIntegrator a
     * {@link ffx.algorithms.Integrator.Integrators} object.
     */
    public MolecularDynamics(MolecularAssembly assembly,
            Potential potentialEnergy,
            CompositeConfiguration properties,
            AlgorithmListener listener,
            Thermostats requestedThermostat,
            Integrators requestedIntegrator) {
        this.molecularAssembly = assembly;
        this.algorithmListener = listener;
        this.potential = potentialEnergy;

        this.properties = properties;
        mass = potentialEnergy.getMass();
        numberOfVariables = potentialEnergy.getNumberOfVariables();
        x = new double[numberOfVariables];
        v = new double[numberOfVariables];
        a = new double[numberOfVariables];
        aPrevious = new double[numberOfVariables];
        grad = new double[numberOfVariables];

        /**
         * If an Integrator wasn't passed to the MD constructor, check for one
         * specified as a property.
         */
        if (requestedIntegrator == null) {
            String integrate = properties.getString("integrate", "beeman").trim();
            try {
                requestedIntegrator = Integrators.valueOf(integrate);
            } catch (Exception e) {
                requestedIntegrator = Integrators.BEEMAN;
            }

        }
        switch (requestedIntegrator) {
            case RESPA:
                integrator = new Respa(numberOfVariables, x, v, a, aPrevious, mass);
                break;
            case STOCHASTIC:
                double friction = properties.getDouble("friction", 91.0);
                Stochastic stochastic = new Stochastic(friction, numberOfVariables, x, v, a, mass);
                if (properties.containsKey("randomseed")) {
                    stochastic.setRandomSeed(properties.getInt("randomseed", 0));
                }
                integrator = stochastic;
                /**
                 * The stochastic dynamics integration procedure will thermostat
                 * the system. The ADIABTIC thermostat just serves to report the
                 * temperature and initialize velocities if necessary.
                 */
                requestedThermostat = Thermostats.ADIABATIC;
                break;
            case VELOCITYVERLET:
                integrator = new VelocityVerlet(numberOfVariables, x, v, a, mass);
                break;
            case BEEMAN:
            default:
                integrator = new BetterBeeman(numberOfVariables, x, v, a, aPrevious, mass);
        }

        /**
         * If a Thermostat wasn't passed to the MD constructor, check for one
         * specified as a property.
         */
        if (requestedThermostat == null) {
            String thermo = properties.getString("thermostat", "Berendsen").trim();
            try {
                requestedThermostat = Thermostats.valueOf(thermo);
            } catch (Exception e) {
                requestedThermostat = Thermostats.BERENDSEN;
            }
        }

        switch (requestedThermostat) {
            case BERENDSEN:
                double tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Berendsen(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
                break;
            case BUSSI:
                tau = properties.getDouble("tau-temperature", 0.2);
                thermostat = new Bussi(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes(), 300.0, tau);
                break;
            case ADIABATIC:
            default:
                thermostat = new Adiabatic(numberOfVariables, x, v, mass, potentialEnergy.getVariableTypes());
        }

        if (properties.containsKey("randomseed")) {
            thermostat.setRandomSeed(properties.getInt("randomseed", 0));
        }

        /**
         * For StochasticDynamics, center of mass motion will not be removed.
         */
        if (integrator instanceof Stochastic) {
            thermostat.removingCenterOfMassMotion(false);
        }

        done = true;
    }
    
    /**
     * Reinitialize the MD engine after a chemical change by specifying the inserted/removed atoms (fast).
     * Currently only works for a SINGLE-ATOM insertion OR deletion.
     * @param insertedAtoms
     * @param removedAtoms 
     */
    public void reInit(List<Atom> insertedAtoms, List<Atom> removedAtoms) {
        int debugLogLevel = 0;
        String debugLogString = System.getProperty("debug");
        if (debugLogString != null) {
            debugLogLevel = Integer.parseInt(debugLogString);
        }
        long startTime = System.nanoTime();
        
        done = false;
        
        mass = potential.getMass();
        int oldVars = numberOfVariables;
        int newVars = potential.getNumberOfVariables();
        int numOldAtoms = oldVars / 3;
        int numNewAtoms = newVars / 3;
        
        // Initialize all arrays with new numberOfVariables.
        double xNew[] = new double[newVars];
        double vNew[] = new double[newVars];
        double aNew[] = new double[newVars];
        double aPreviousNew[] = new double[newVars];
        double gradNew[] = new double[newVars];

        Atom newAtomArray[] = molecularAssembly.getAtomArray();
        if (numNewAtoms != newAtomArray.length) {
            logger.warning(String.format("numNewAtoms: %d\tnewAtomArray.length: %d", numNewAtoms, newAtomArray.length));
        }
        
        if (insertedAtoms.isEmpty() && !removedAtoms.isEmpty()) {
            if (removedAtoms.size() > 1) {
                logger.severe("MD reInit(): option not yet implemented.");
            }
            // We lost one atom from the system.
            // Locate its position and shift all x,v,etc. array elements up from that point on.
            Atom removed = removedAtoms.get(0);
            int remPosition = (removed.getXYZIndex() - 1) * 3;
            for (int i = 0; i < newVars; i++) {
                if (i < remPosition) {
                    xNew[i] = x[i];
                    vNew[i] = v[i];
                    aNew[i] = a[i];
                    aPreviousNew[i] = aPrevious[i];
                    gradNew[i] = grad[i];
                } else {
                    xNew[i] = x[i+3];
                    vNew[i] = v[i+3];
                    aNew[i] = a[i+3];
                    aPreviousNew[i] = aPrevious[i+3];
                    gradNew[i] = grad[i+3];
                }
            }            
        } else if (!insertedAtoms.isEmpty() && removedAtoms.isEmpty()) {
            if (insertedAtoms.size() > 1) {
                logger.severe("MD reInit(): option not yet implemented.");
            }
            // We gained one new atom in the system.
            // Locate its position and insert its x,v,etc. array elements there, shifting all others' down.
            Atom inserted = insertedAtoms.get(0);
            int insPosition = (inserted.getXYZIndex() - 1) * 3;
            for (int i = 0; i < oldVars; i++) {
                if (i < insPosition) {
                    xNew[i] = x[i];
                    vNew[i] = v[i];
                    aNew[i] = a[i];
                    aPreviousNew[i] = aPrevious[i];
                    gradNew[i] = grad[i];
                } else {
                    xNew[i+3] = x[i];
                    vNew[i+3] = v[i];
                    aNew[i+3] = a[i];
                    aPreviousNew[i+3] = aPrevious[i];
                    gradNew[i+3] = grad[i];
                }
            }
        } else {
            // Not yet implemented.
            logger.severe("MD reInit() can't yet handle adding and removing atoms simultaneously.");
        }
        
        // Initialize new-atom velocities from a Maxwell distribution.
        for (Atom inserted : insertedAtoms) {
            int insPosition = (inserted.getXYZIndex()-1) * 3;
            double maxwell[] = thermostat.maxwellIndividual(inserted.getMass());
            vNew[insPosition] = maxwell[0];
            vNew[insPosition+1] = maxwell[1];
            vNew[insPosition+2] = maxwell[2];
            logger.info(" Assigned maxwell velocity for new atom " + inserted);
        }
        
        numberOfVariables = newVars;
        x = xNew;
        v = vNew;
        a = aNew;
        aPrevious = aPreviousNew;
        grad = gradNew;
        
        potential.getCoordinates(x);
        thermostat.setNumberOfVariables(newVars, x, v, mass);
        ((BetterBeeman) integrator).setNumberOfVariables(newVars, x, v, a, aPrevious, mass);

        long took = (long) ((System.nanoTime() - startTime) * 1e-6);
        // logger.info(String.format(" MD reInit(): %d ms", took));
        
        done = true;
        return;
    }
    
    /**
     * Reinitialize the MD engine after an arbitrary chemical change.
     * WARNING: this works but it's god-awful slow (like two-minutes-per slow).
     * @param oldAtomArray a deep copy of the all-atom array from before you made a change
     */
    public void reInit(Atom oldAtomArray[]) {
        long startTime = System.nanoTime();
        int debugLogLevel = 0;
        String debugLogString = System.getProperty("debug");
        if (debugLogString != null) {
            debugLogLevel = Integer.parseInt(debugLogString);
        }
        if (debugLogLevel > 0) {
            logger.info(String.format(" ---------\n REINIT MD\n ---------"));
        }
        
        done = false;
        if (oldAtomArray == null) {
            logger.severe("Old atom array was null.");
        }
        
        mass = potential.getMass();
        int oldVars = numberOfVariables;
        int newVars = potential.getNumberOfVariables();
        int numOldAtoms = oldVars / 3;
        int numNewAtoms = newVars / 3;
        
        // Initialize all arrays with new numberOfVariables.
        double xNew[] = new double[newVars];
        double vNew[] = new double[newVars];
        double aNew[] = new double[newVars];
        double aPreviousNew[] = new double[newVars];
        double gradNew[] = new double[newVars];
        
        int smallVars = (oldVars < newVars ? oldVars : newVars);
        for (int i = 0; i < smallVars; i++) {
            xNew[i] = x[i];
            vNew[i] = v[i];
            aNew[i] = a[i];
            aPreviousNew[i] = aPrevious[i];
            gradNew[i] = grad[i];
        }

        Atom newAtomArray[] = molecularAssembly.getAtomArray();
        if (numNewAtoms != newAtomArray.length) {
            logger.warning(String.format("numNewAtoms: %d\tnewAtomArray.length: %d", numNewAtoms, newAtomArray.length));
        }
        
        // for each new atom: for each old atom: if they match, copy the x v a aPrev grad arrays from old to new;
        //      flag those new atoms which don't have an old component as 'halted' (to be given maxwell velocity);
        //      keep track of oldIndex values which have already been match for quicker looping
        List<Atom> matchedAtoms = new ArrayList<>();
        List<Atom> haltedAtoms = new ArrayList<>();
        int oldIndex, newIndex;
        for (newIndex = 0; newIndex < numNewAtoms; newIndex++) {
            Atom newAtom = newAtomArray[newIndex];
            String newString = newAtom.toNameNumberString();
            boolean found = false;
            for (oldIndex = 0; oldIndex < numOldAtoms; oldIndex++) {
                Atom oldAtom = oldAtomArray[oldIndex];
                String oldString = oldAtom.toNameNumberString();
                if (newAtom == oldAtom || newString.equals(oldString)) {
                    found = true;
                    int newPosition = newIndex*3, oldPosition = oldIndex*3;
                    if (newPosition == oldPosition) {
                        break;
                    }
                    /*
                    StringBuilder sb = new StringBuilder();
                    sb.append(String.format(" old:%d  new:%d  atom:%s\n", oldIndex, newIndex, oldAtom));
                    sb.append(String.format("   xva: [%.2g %.2g %.2g] [%.2g %.2g %.2g] [%.2g %.2g %.2g]",
                            x[oldPosition], x[oldPosition+1], x[oldPosition+2],
                            v[oldPosition], v[oldPosition+1], v[oldPosition+2],
                            a[oldPosition], a[oldPosition+1], a[oldPosition+2]));
                    if (newIndex <= 20) {
                        logger.info(sb.toString());
                    }*/                    
                    xNew[newPosition] = x[oldPosition];
                    xNew[newPosition+1] = x[oldPosition+1];
                    xNew[newPosition+2] = x[oldPosition+2];
                    vNew[newPosition] = v[oldPosition];
                    vNew[newPosition+1] = v[oldPosition+1];
                    vNew[newPosition+2] = v[oldPosition+2];
                    aNew[newPosition] = a[oldPosition];
                    aNew[newPosition+1] = a[oldPosition+1];
                    aNew[newPosition+2] = a[oldPosition+2];
                    aPreviousNew[newPosition] = aPrevious[oldPosition];
                    aPreviousNew[newPosition+1] = aPrevious[oldPosition+1];
                    aPreviousNew[newPosition+2] = aPrevious[oldPosition+2];
                    gradNew[newPosition] = grad[oldPosition];
                    gradNew[newPosition+1] = grad[oldPosition+1];
                    gradNew[newPosition+2] = grad[oldPosition+2];
                    break;
                }
            }
            if (!found) {
                logger.info(String.format(" No previous record of atom: %s.", newAtom));
                haltedAtoms.add(newAtom);
            }
        }
        
        for (Atom halted : haltedAtoms) {
            // Using a maxwell distribution.
            int haltedPosition = (halted.getXYZIndex() - 1) * 3;
            double maxwell[] = thermostat.maxwellIndividual(halted.getMass());
            vNew[haltedPosition] = maxwell[0];
            vNew[haltedPosition+1] = maxwell[1];
            vNew[haltedPosition+2] = maxwell[2];
            logger.info(" Assigned maxwell velocity for new atom " + halted);
            
            /* Using a nearby bonded atom (not stat-mech rigorous).
            try {
                Bond bond = (Bond) halted.getBondList().get(0);
                Atom mate = bond.get1_2(halted);
                int matePosition = (mate.getXYZIndex() - 1) * 3;
                int haltedPosition = (halted.getXYZIndex() - 1) * 3;
                vNew[haltedPosition] = vNew[matePosition];
                vNew[haltedPosition+1] = vNew[matePosition+1];
                vNew[haltedPosition+2] = vNew[matePosition+2];
                aNew[haltedPosition] = aNew[matePosition];
                aNew[haltedPosition+1] = aNew[matePosition+1];
                aNew[haltedPosition+2] = aNew[matePosition+2];
                logger.info(" Assigned velocity for new atom " + halted + " using bond with " + mate);
            } catch (IndexOutOfBoundsException ex) {
                logger.info(" No bond available to init velocity for atom " + halted);
            }
            */
        }
        
        numberOfVariables = newVars;
        x = xNew;
        v = vNew;
        a = aNew;
        aPrevious = aPreviousNew;
        grad = gradNew;
        
        potential.getCoordinates(x);
        thermostat.setNumberOfVariables(newVars, x, v, mass);
        ((BetterBeeman) integrator).setNumberOfVariables(newVars, x, v, a, aPrevious, mass);

        long took = (long) ((System.nanoTime() - startTime) * 1e-6);
        logger.info(String.format(" MD reInit(): %d ms", took));
        
        done = true;
        return;
    }

    /**
     * <p>
     * Setter for the field <code>thermostat</code>.</p>
     *
     * @param thermostat a {@link ffx.algorithms.Thermostat} object.
     */
    public void setThermostat(Thermostat thermostat) {
        this.thermostat = thermostat;
    }

    /**
     * <p>
     * Getter for the field <code>thermostat</code>.</p>
     *
     * @return a {@link ffx.algorithms.Thermostat} object.
     */
    public Thermostat getThermostat() {
        return thermostat;
    }

    /**
     * <p>
     * Setter for the field <code>x</code>.</p>
     *
     * @param x a double array to set the current parameters to.
     */
    public void setParameters(double x[]) {
        System.arraycopy(x, 0, this.x, 0, numberOfVariables);
    }

    /**
     * <p>
     * Getter for the field <code>x</code>.</p>
     *
     * @return a double array with the current parameters
     */
    public double[] getParameters() {
        return x;
    }

    /**
     * <p>
     * Setter for the field <code>archiveFile</code>.</p>
     *
     * @param archive a {@link java.io.File} object.
     */
    public void setArchiveFile(File archive) {
        this.archiveFile = archive;
    }

    /**
     * <p>
     * Getter for the field <code>archiveFile</code>.</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getArchiveFile() {
        return archiveFile;
    }
    
    public void addMCListener(MonteCarloListener monteCarloListener) {
        this.monteCarloListener = monteCarloListener;
    }

    /**
     * <p>
     * init</p>
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param fileType a String.
     * @param restartFrequency the number of steps between writing restart
     * files.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final String fileType, final double restartFrequency,
            final double temperature, final boolean initVelocities, final File dyn) {

        /**
         * Return if already running.
         */
        if (!done) {
            logger.warning(" Programming error - attempt to modify parameters of a running MolecularDynamics instance.");
            return;
        }

        this.nSteps = nSteps;
        /**
         * Convert the time step from femtoseconds to picoseconds.
         */
        dt = timeStep * 1.0e-3;

        /**
         * Convert the print interval to a print frequency.
         */
        printFrequency = 100;
        if (printInterval >= this.dt) {
            printFrequency = (int) (printInterval / this.dt);
        }

        /**
         * Convert the save interval to a save frequency.
         */
        saveSnapshotFrequency = 1000;
        if (saveInterval >= this.dt) {
            saveSnapshotFrequency = (int) (saveInterval / this.dt);
        }

        /**
         * Set snapshot file type.
         */
        saveSnapshotAsPDB = true;
        if (fileType.equals("XYZ")) {
            saveSnapshotAsPDB = false;
        } else if (!fileType.equals("PDB")) {
            logger.warning("Snapshot file type unrecognized; saving snaphshots as PDB.\n");
        }

        /**
         * Convert restart frequency to steps.
         */
        saveRestartFileFrequency = 1000;
        if (restartFrequency >= this.dt) {
            saveRestartFileFrequency = (int) (restartFrequency / this.dt);
        }

        File file = molecularAssembly.getFile();
        String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
        if (archiveFile == null) {
            archiveFile = new File(filename + ".arc");
            archiveFile = XYZFilter.version(archiveFile);
        }

        if (dyn == null) {
            this.restartFile = new File(filename + ".dyn");
            loadRestart = false;
        } else {
            this.restartFile = dyn;
            loadRestart = true;
        }

        if (xyzFilter == null) {
            xyzFilter = new XYZFilter(file, molecularAssembly,
                    molecularAssembly.getForceField(), properties);
        }

        if (dynFilter == null) {
            dynFilter = new DYNFilter(molecularAssembly.getName());
        }

        if (pdbFilter == null) {
            pdbFile = new File(filename + "_dyn.pdb");
            pdbFilter = new PDBFilter(new File(filename + "_dyn.pdb"), molecularAssembly, null, null);
        }

        this.targetTemperature = temperature;
        this.initVelocities = initVelocities;
    }

    /**
     * A version of init with the original method header. Redirects to the new
     * method with default values for added parameters. Needed by (at least)
     * ReplicaExchange, which calls this directly.
     *
     * @param nSteps the number of MD steps.
     * @param timeStep the time step.
     * @param printInterval the number of steps between loggging updates.
     * @param saveInterval the number of steps between saving snapshots.
     * @param temperature the target temperature.
     * @param initVelocities true to reset velocities from a Maxwell
     * distribution.
     * @param dyn the Dynamic restart file.
     */
    public void init(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            final File dyn) {
        init(nSteps, timeStep, printInterval, saveInterval, "PDB", 0.1, temperature, initVelocities, dyn);
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is
     * done.
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param fileType a String (either XYZ or PDB).
     * @param restartFrequency a double specifying the restart frequency.
     * @param dyn a {@link java.io.File} object.
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            String fileType, double restartFrequency, final File dyn) {
        this.fileType = fileType;
        this.restartFrequency = restartFrequency;
        dynamic(nSteps, timeStep, printInterval, saveInterval, temperature,
                initVelocities, dyn);
    }

    /**
     * Blocking molecular dynamics. When this method returns, the MD run is
     * done.
     *
     * @param nSteps a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void dynamic(final int nSteps, final double timeStep, final double printInterval,
            final double saveInterval, final double temperature, final boolean initVelocities,
            final File dyn) {

        /**
         * Return if already running; Could happen if two threads call dynamic
         * on the same MolecularDynamics instance.
         */
        if (!done) {
            logger.warning(" Programming error - a thread invoked dynamic when it was already running.");
            return;
        }

        if (integrator instanceof Stochastic) {
            logger.info(format("\n Stochastic dynamics in the NVT ensemble\n"));
        } else if (!(thermostat instanceof Adiabatic)) {
            logger.info(format("\n Molecular dynamics in the NVT ensemble\n"));
        } else {
            logger.info(format("\n Molecular dynamics in the NVE ensemble\n"));
        }

        init(nSteps, timeStep, printInterval, saveInterval, fileType, restartFrequency, temperature, initVelocities, dyn);

        done = false;

        if (dyn != null) {
            logger.info(format(" Continuing from " + dyn.getAbsolutePath()));
        }
        logger.info(String.format(" Number of steps: %8d", nSteps));
        logger.info(String.format(" Time step:       %8.3f (fsec)", timeStep));
        logger.info(String.format(" Print interval:  %8.3f (psec)", printInterval));
        logger.info(String.format(" Save interval:   %8.3f (psec)", saveInterval));
        logger.info(String.format(" Archive file: %s", archiveFile.getName()));
        logger.info(String.format(" Restart file: %s", restartFile.getName()));

        Thread dynamicThread = new Thread(this);
        dynamicThread.start();
        synchronized (this) {
            try {
                while (dynamicThread.isAlive()) {
                    wait(100);
                }
            } catch (InterruptedException e) {
                String message = " Molecular dynamics interrupted.";
                logger.log(Level.WARNING, message, e);
            }
        }
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
     * Set the number of time steps between removal of center of mass kinetic
     * energy.
     *
     * @param removeCOMMotionFrequency Number of time steps between center of
     * mass removal.
     */
    public void setRemoveCOMMotionFrequency(int removeCOMMotionFrequency) {
        if (removeCOMMotionFrequency < 0) {
            removeCOMMotionFrequency = 0;
        }
        this.removeCOMMotionFrequency = removeCOMMotionFrequency;
        if (removeCOMMotionFrequency != 0) {
            thermostat.removingCenterOfMassMotion(true);
        } else {
            thermostat.removingCenterOfMassMotion(false);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        done = false;
        terminate = false;

        /**
         * Set the target temperature.
         */
        thermostat.setTargetTemperature(targetTemperature);
        if (integrator instanceof Stochastic) {
            Stochastic stochastic = (Stochastic) integrator;
            stochastic.setTemperature(targetTemperature);
        }

        /**
         * Set the step size.
         */
        integrator.setTimeStep(dt);

        if (!initialized) {
            /**
             * Initialize from a restart file.
             */
            if (loadRestart) {
                if (!dynFilter.readDYN(restartFile, x, v, a, aPrevious)) {
                    String message = " Could not load the restart file - dynamics terminated.";
                    logger.log(Level.WARNING, message);
                    done = true;
                    return;
                }
            } else {
                /**
                 * Initialize from using current atomic coordinates.
                 */
                potential.getCoordinates(x);
                /**
                 * Initialize atomic velocities from a Maxwell-Boltzmann
                 * distribution or set to 0.
                 */
                if (initVelocities) {
                    thermostat.maxwell(targetTemperature);
                } else {
                    fill(v, 0.0);
                }
            }
        } else {
            /**
             * If MD has already been run (ie. Annealing or RepEx), then
             * initialize velocities if requested.
             */
            if (initVelocities) {
                thermostat.maxwell(targetTemperature);
            }
        }

        /**
         * Compute the current potential energy.
         */
        potential.setScaling(null);
        currentPotentialEnergy = potential.energyAndGradient(x, grad);

        /**
         * Compute the current kinetic energy.
         */
        thermostat.kineticEnergy();
        currentKineticEnergy = thermostat.getKineticEnergy();
        currentTemperature = thermostat.getCurrentTemperature();
        currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

        /**
         * Initialize current and previous accelerations.
         */
        if (!initialized) {
            if (!loadRestart) {
                for (int i = 0; i < numberOfVariables; i++) {
                    a[i] = -Thermostat.convert * grad[i] / mass[i];
                }
                if (aPrevious != null) {
                    System.arraycopy(a, 0, aPrevious, 0, numberOfVariables);
                }
            }
            initialized = true;
        }

        logger.info(String.format("\n      Time      Kinetic    Potential        Total     Temp      CPU"));
        logger.info(String.format("      psec     kcal/mol     kcal/mol     kcal/mol        K      sec\n"));
        logger.info(String.format("          %13.4f%13.4f%13.4f %8.2f ", currentKineticEnergy, currentPotentialEnergy, currentTotalEnergy, currentTemperature));

        /**
         * Integrate Newton's equations of motion for the requested number of
         * steps, unless early termination is requested.
         */
        long time = System.nanoTime();
        for (int step = 1; step <= nSteps; step++) {
            
            /**
             * Do the half-step thermostat operation.
             */
            thermostat.halfStep(dt);

            /**
             * Do the half-step integration operation.
             */
            integrator.preForce(potential);

            /**
             * Compute the potential energy and gradients.
             */
            currentPotentialEnergy = potential.energyAndGradient(x, grad);

            /**
             * Add the potential energy of the slow degrees of freedom.
             */
            if (integrator instanceof Respa) {
                Respa r = (Respa) integrator;
                currentPotentialEnergy += r.getHalfStepEnergy();
            }

            /**
             * Do the full-step integration operation.
             */
            integrator.postForce(grad);

            /**
             * Compute the full-step kinetic energy.
             */
            thermostat.kineticEnergy();

            /**
             * Do the full-step thermostat operation.
             */
            thermostat.fullStep(dt);

            /**
             * Recompute the kinetic energy after the full-step thermostat
             * operation.
             */
            thermostat.kineticEnergy();

            /**
             * Remove center of mass motion ever ~100 steps.
             */
            if (thermostat.removingCOM() && step % removeCOMMotionFrequency == 0) {
                thermostat.centerOfMassMotion(true, false);
            }

            /**
             * Collect current kinetic energy, temperature, and total energy.
             */
            currentKineticEnergy = thermostat.getKineticEnergy();
            currentTemperature = thermostat.getCurrentTemperature();
            currentTotalEnergy = currentKineticEnergy + currentPotentialEnergy;

            /**
             * Log the current state every printFrequency steps.
             */
            if (step % printFrequency == 0) {
                double simTime = step * dt;
                time = System.nanoTime() - time;
                logger.info(String.format(" %7.3e%13.4f%13.4f%13.4f%9.2f%9.3f", simTime, currentKineticEnergy, currentPotentialEnergy,
                        currentTotalEnergy, currentTemperature, time * 1.0e-9));
                time = System.nanoTime();
            }

            /**
             * Write out snapshots in selected format every
             * saveSnapshotFrequency steps.
             */
            if (saveSnapshotFrequency > 0 && step % saveSnapshotFrequency == 0) {
                if (archiveFile != null && saveSnapshotAsPDB == false) {
                    if (xyzFilter.writeFile(archiveFile, true)) {
                        logger.info(String.format(" Appended snap shot to " + archiveFile.getName()));
                    } else {
                        logger.warning(String.format(" Appending snap shot to " + archiveFile.getName() + " failed"));
                    }
                } else if (saveSnapshotAsPDB == true) {
                    if (pdbFilter.writeFile(pdbFile, false)) {
                        logger.info(String.format(" Wrote PDB file to " + pdbFile.getName()));
                    }
                }
            }

            /**
             * Write out restart files every saveRestartFileFrequency steps.
             */
            if (saveRestartFileFrequency > 0 && step % saveRestartFileFrequency == 0) {
                if (dynFilter.writeDYN(restartFile, molecularAssembly.getCrystal(), x, v, a, aPrevious)) {
                    logger.info(String.format(" Wrote dynamics restart file to " + restartFile.getName()));
                } else {
                    logger.info(String.format(" Writing dynamics restart file to " + restartFile.getName() + " failed"));
                }
            }

            /**
             * Notify the algorithmListener.
             */
            if (algorithmListener != null && step % printFrequency == 0) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            if (monteCarloListener != null) {
                long startTime = System.nanoTime();
                monteCarloListener.mcUpdate(molecularAssembly);
                potential.getCoordinates(x);
                long took = (long) ((System.nanoTime() - startTime) * 1e-6);
//                logger.info(String.format(" mcUpdate() took: %d ms", took));
                
                /* TESTING
                for (Atom atom : molecularAssembly.getAtomArray()) {
                    if (atom.toString().contains("HE2 GLH 27")) {
                        int index = (atom.getXYZIndex() - 1)*3;
                        logger.info(String.format(" xyz: %d    %s\n     [%g %g %g] [%g %g %g]\n", 
                                atom.xyzIndex, atom, v[index], v[index+1], v[index+2],
                                a[index], a[index+1], a[index+2]));
                    }
                } */
            }

            /**
             * Check for a termination request.
             */
            if (terminate) {
                logger.info(String.format("\n Terminating after %8d time steps\n", step));
                break;
            }
        }

        /**
         * Log normal completion.
         */
        if (!terminate) {
            logger.info(String.format(" Completed %8d time steps\n", nSteps));
        }

        /**
         * Reset the done and terminate flags.
         */
        done = true;
        terminate = false;
    }

    /**
     * Get the total system energy.
     *
     * @return total energy.
     */
    public double getTotalEnergy() {
        return currentTotalEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating dynamics.\n", e);
                }
            }
        }
    }
}
