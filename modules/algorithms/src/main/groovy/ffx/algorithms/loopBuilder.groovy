/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

// LOOP BUILDER

// Java Imports
import java.util.Scanner;

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Paralle Java Imports
import edu.rit.pj.Comm;
import edu.rit.mp.DoubleBuf;
import edu.rit.mp.IntegerBuf;

// FFX Imports
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Minimize;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.OSRW;
import ffx.algorithms.TransitionTemperedOSRW;
import ffx.algorithms.RotamerOptimization
import ffx.algorithms.RotamerOptimization.Direction;
import ffx.algorithms.SimulatedAnnealing;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.algorithms.MCLoop;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.ResidueEnumerations;
import ffx.potential.bonded.ResidueEnumerations.CommonAminoAcid3;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.ForceFieldEnergy;
import ffx.numerics.Potential;

// Default convergence criteria.
double eps = 0.1;

// Temperture in degrees Kelvin.
double temperature = 298.15;

// Time step in femtoseconds.
double timeStep = 1.0;

// Frequency to log thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to write out coordinates in picoseconds.
double saveInterval = 100;

// Frequency to write out restart information in picoseconds.
double restartInterval = 1.0;

// Number of molecular dynamics steps: default is 100 nanoseconds.
int nSteps = 50000;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = Thermostats.BERENDSEN;

// Integrators [ BEEMAN, RESPA, STOCHASTIC ]
Integrators integrator = Integrators.BEEMAN;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// File type of coordinate snapshots to write out.
String fileType = "PDB";

// Value of Lambda.
double lambda = 0.0;

// Monte-Carlo step frequencies for loop moves.
int mcStepFrequency = 1000;

// Rotamer Optimization
boolean runRotamer = false;

// Simulated Annealing
boolean runSimulatedAnnealing = false;

// OSRW
boolean runOSRW = true;

// Monte Carlo with KIC
boolean runMCLoop = false;

// Transition Tempered OSRW
boolean runTTOSRW = false;

// Local minimization mode
boolean localMin = false;

RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary();
// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc loopBuilder [options] <filename1>');
cli.h(longOpt:'help', 'Print this help message.');
cli.e(longOpt:'eps', args:1, argName:'1.0', 'RMS gradient convergence criteria');
cli.n(longOpt:'steps', args:1, argName:'10000', 'Number of molecular dynamics steps.');
cli.c(longOpt:'chain', args:1, argName:' ', 'Single character chain name to limit optimization to single chain.');
cli.d(longOpt:'dt', args:1, argName:'2.5', 'Time discretization step (fsec).');
cli.i(longOpt:'integrator', args:1, argName:'Beeman', 'Integrator: [Beeman / Respa / Stochastic]');
cli.r(longOpt:'report', args:1, argName:'0.01', 'Interval to report thermodyanamics (psec).');
cli.m(longOpt:'minimize','Local minimization of loop residues (need -s and -f flags).');
cli.w(longOpt:'write', args:1, argName:'100.0', 'Interval to write out coordinates (psec).');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.g(longOpt:'bias', args:1, argName:'0.01', 'Gaussian bias magnitude (kcal/mol).');
cli.osrw(longOpt:'OSRW', 'Run OSRW.');
cli.tt(longOpt:'ttOSRW', 'Run Transition Tempered OSRW');
cli.sa(longOpt:'simulated annealing', 'Run simulated annealing.');
cli.rot(longOpt:'rotamer', 'Run rotamer optimization.');
cli.mc(longOpt:'MC Loop','Run Monte Carlo KIC');
cli.a(longOpt:'all', 'Run optimal pipeline of algorithms.');
cli.s(longOpt:'start', args:1, argName:'1', 'Starting residue of existing loop.');
cli.f(longOpt:'final', args:1, argName:'-1', 'Final residue of an existing loop.');
cli.mcn(longOpt:'mcStepFreq', args:1, argName:'10', 'Number of MD steps between Monte-Carlo protonation changes.')

def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Starting and ending loop atoms.
if (options.s && options.f) {
    loopStart = Integer.parseInt(options.s);
    loopStop = Integer.parseInt(options.f);
} else if (options.s || options.f){
    logger.info("Starting atom and final atom numbers are need to use this option.");
}

// Chain Name.
if (options.c) {
    chain = options.c.toCharacter();
}

// Local minimization mode
if (options.m) {
    localMin = true;
}

// Set Monte Carlo step frequency
if (options.mcn) {
    mcStepFrequency = Integer.parseInt(options.mcn);
}

// Load the time steps in femtoseconds.
if (options.d) {
    timeStep = Double.parseDouble(options.d);
}

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = Integrators.BEEMAN;
    }
}

// Report interval in picoseconds.
if (options.r) {
    printInterval = Double.parseDouble(options.r);
}

// Write interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Load convergence criteria.
if (options.e) {
    eps = Double.parseDouble(options.e);
}

// Run OSRW
if (options.osrw){
    runOSRW = true;
}

// Gaussian bias magnitude (kcal/mol).
if (options.g) {
    biasMag = Double.parseDouble(options.g);
}

// Run Simulated Annealing
if (options.sa){
    runSimulatedAnnealing = true;
}

// Run Rotamer Optimization
if (options.rot){
    runRotamer = true;
}

// Run Transition Tempered OSRW
if (options.tt){
    runTTOSRW = true;
}

// Default
if (!(options.osrw && options.sa)){
    runOSRW = true;
}

// Run MC Loop Optimization
if (options.mc){
    runMCLoop = true;
  //  runOSRW = false;
    MCLoop mcLoop;
}

// Robust Default
if (options.a){
    runOSRW = true;
    runRotamer = true;
}
//build loop with PDBFilter if an existing loop is not provided
if(!(options.s && options.f)){
    System.setProperty("buildLoops", "true");
}
System.setProperty("vdwterm", "false");

List<String> arguments = options.arguments();
String filename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    filename = arguments.get(0);
} else {
    return cli.usage();
}

File structureFile = new File(FilenameUtils.normalize(filename));
structureFile = new File(structureFile.getAbsolutePath());
String baseFilename = FilenameUtils.removeExtension(structureFile.getName());
File histogramRestart = new File(baseFilename + ".his");
File lambdaRestart = new File(baseFilename + ".lam");
File dyn = new File(baseFilename + ".dyn");

Comm world = Comm.world();
int size = world.size();
int rank = 0;
double[] energyArray = new double[world.size()];
for(int i =0; i < world.size(); i++){
    energyArray[i] = Double.MAX_VALUE;
}

// For a multi-process job, try to get the restart files from rank sub-directories.
if (size > 1) {
    rank = world.rank();
    File rankDirectory = new File(structureFile.getParent() + File.separator + Integer.toString(rank));
    if (!rankDirectory.exists()) {
        rankDirectory.mkdir();
    }
    lambdaRestart = new File(rankDirectory.getPath() + File.separator + baseFilename + ".lam");
    dyn = new File(rankDirectory.getPath() + File.separator + baseFilename + ".dyn");
    structureFile = new File(rankDirectory.getPath() + File.separator + structureFile.getName());
}

if (!dyn.exists()) {
    dyn = null;
}

open(filename);

// If this is a multi-process job, set the structure file to come from the subdirectory.
if (size > 1) {
    active.setFile(structureFile);
}

// Get a reference to the first system's ForceFieldEnergy.
ForceFieldEnergy forceFieldEnergy = active.getPotentialEnergy();
forceFieldEnergy.setPrintOnFailure(false, false);
// Set built atoms active/use flags to true (false for other atoms).
Atom[] atoms = active.getAtomArray();

//if existing loop is used, set loop atoms to match atoms built with PDBFilter
if(options.s && options.f){
    for (int i = 0; i < atoms.length; i++){
        Atom ai = atoms[i];
        if(!options.c || chain == ai.getChainID()){
            if(ai.getResidueNumber() >= loopStart && ai.getResidueNumber() <= loopStop){
                ai.setBuilt(true);
            }
        }
    }
} else {
    //create array of built residues
    ArrayList<Residue> loopResidues = new ArrayList<>();
    for (int i = 0; i < active.getChains().size(); i++){
        ArrayList<Residue> allResidues = active.getChains()[i].getResidues();

        for (int j = 0; j < allResidues.size(); j++) {
            Residue temp = allResidues[j];
            if (temp.getBackboneAtoms().get(0).getBuilt()) {
                loopResidues.add(allResidues[j]);
            }
        }
    }
    loopStart = loopResidues.get(0).getResidueNumber();
    loopStop = loopResidues.get(0).getResidueNumber();
    for (int i = 0; i < loopResidues.size(); i++){
        if(loopResidues.get(i).getChainID() == loopResidues.get(0).getChainID()){
            if (loopStop + 1 == loopResidues.get(i).getResidueNumber()){
                loopStop = loopResidues.get(i).getResidueNumber();
            }
        }
    }
}

for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    if (ai.getBuilt()) {
        ai.setActive(true);
        ai.setUse(true);
    } else {
        ai.setActive(false);
        ai.setUse(true);
    }
}

logger.info("\n Running minimize on built atoms of " + active.getName());
logger.info(" RMS gradient convergence criteria: " + eps);

// Minimization without vdW.
if (localMin){
    runOSRW = false;
} else {
    e = minimize(eps);
    energy();
}

boolean loopBuildError = false;

if(runOSRW){
    // Run OSRW.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "true");
    System.setProperty("polarization", "none");
    System.setProperty("intramolecular-softcore", "true");
    System.setProperty("intermolecular-softcore", "true");
    System.setProperty("lambdaterm", "true");
    System.setProperty("lambda-torsions", "true");
    System.setProperty("ligand-vapor-elec","false");
    System.setProperty("vdw-cutoff", "9.0");
    System.setProperty("lambda-bias-cutoff", "3");
    if (options.g) {
        System.setProperty("bias-gaussian-mag",String.format("%f",biasMag));
    } else {
        System.setProperty("bias-gaussian-mag", "0.002");
    }

    // Set the thermostat time constant (in psec) to the the time step (i.e. to give velocity rescaling).
    System.setProperty("tau-temperature", String.format("%f", timeStep * 1.0e-3));

    System.setProperty("lambda-bin-width", "0.01");
    System.setProperty("tau-temperature","0.05");

    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        if (ai.getBuilt()) {
            ai.setApplyLambda(true);
        } else {
            ai.setApplyLambda(false);
        }
    }

    forceFieldEnergy= new ForceFieldEnergy(active);
    forceFieldEnergy.setPrintOnFailure(false, false);
    forceFieldEnergy.setLambda(lambda);
    energy();

    // Turn off checks for overlapping atoms, which is expected for lambda=0.
    forceFieldEnergy.getCrystal().setSpecialPositionCutoff(0.0);

    boolean asynchronous = true;

    Potential osrw;
    if(runTTOSRW){
        osrw = new TransitionTemperedOSRW(forceFieldEnergy, forceFieldEnergy, lambdaRestart, histogramRestart, active.getProperties(),
            (temperature), timeStep, printInterval, saveInterval, asynchronous, sh);
    } else {
        osrw =  new OSRW(forceFieldEnergy, forceFieldEnergy, lambdaRestart, histogramRestart, active.getProperties(),
            (temperature), timeStep, printInterval, saveInterval, asynchronous, sh);
    }

    osrw.setLambda(lambda);
    osrw.setThetaMass(1.0e-19);
    osrw.setOptimization(true, active);
    // Create the MolecularDynamics instance.
    MolecularDynamics molDyn = new MolecularDynamics(active, osrw, active.getProperties(),
        null, thermostat, integrator);


    if(runMCLoop){
        mcLoop = new MCLoop(active, mcStepFrequency, molDyn.getThermostat(),loopStart,loopStop);
        molDyn.addMCListener(mcLoop);
        mcLoop.addMolDyn(molDyn);
        mcLoop.addLambdaInterface(osrw.getLambdaInterface());
        mcLoop.setIterations(20);

    }


    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,
        fileType, restartInterval, dyn);

    logger.info("Obtaining low energy coordinates");
    double[] lowEnergyCoordinates = osrw.getLowEnergyLoop();
    double currentOSRWOptimum = osrw.getOSRWOptimum();
    if (lowEnergyCoordinates != null){
        forceFieldEnergy.setCoordinates(lowEnergyCoordinates);
    } else {
        logger.info("OSRW stage did not succeed in finding a loop.");
        loopBuildError = true;
    }
}
 /*
if (runMCLoop){
    // Monte Carlo with KIC
    System.setProperty("vdwterm", "false");
    System.setProperty("polarization", "none");
    logger.info("\n Running molecular dynamics on " + baseFilename);

    forceFieldEnergy= new ForceFieldEnergy(active);
    // create the MD object
    MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), active.getProperties(), null, thermostat, integrator);

    // create the Monte-Carlo listener and connect it to the MD
    mcLoop = new MCLoop(active, mcStepFrequency, molDyn.getThermostat(),loopStart,loopStop);
    molDyn.addMCListener(mcLoop);
    mcLoop.addMolDyn(molDyn);

    molDyn.dynamic(nSteps, timeStep, printInterval, saveInterval, temperature, initVelocities,fileType,restartInterval,dyn);
}   */


if (runSimulatedAnnealing) {
    // Minimize with vdW.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "false");
    forceFieldEnergy = new ForceFieldEnergy(active);
    forcefieldEnergy.setPrintOnFailure(false, false);
    e = minimize(eps);

    // SA with vdW.
    logger.info("\n Running simulated annealing on " + active.getName());
    double[] heatUpTemperatures = [150,250,400,700,1000];
    // Number of molecular dynamics steps at each temperature.
    int steps = 267; //267 at 3
    // Time step in femtoseconds.
    timeStep = 3.0;
    // Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
    thermostat = Thermostats.BERENDSEN;
    // Integrators [ BEEMAN, RESPA, STOCHASTIC]
    if (!options.i) {
        integrator = Integrators.RESPA;
    }

    SimulatedAnnealing simulatedAnnealing = new SimulatedAnnealing(active, forceFieldEnergy, active.getProperties(), null, thermostat, integrator);
    simulatedAnnealing.annealToTargetValues(heatUpTemperatures, steps, timeStep);

    double[] annealingTargetTemperatures = [1000, 800, 600, 500, 400, 300];
    steps = 800; //800 at 3
    simulatedAnnealing.annealToTargetValues(annealingTargetTemperatures,steps,timeStep);
}


for (int i = 0; i <= atoms.length; i++) {
    Atom ai = atoms[i - 1];
    ai.setUse(true);
    ai.setApplyLambda(false);
}

if(!loopBuildError){
    // Optimize with the full AMOEBA potential energy.
    System.setProperty("vdwterm", "true");
    System.setProperty("mpoleterm", "true");
    System.setProperty("polarization", "direct");
    System.setProperty("intramolecularSoftcore", "false");
    System.setProperty("intermolecularSoftcore", "false");
    System.setProperty("lambdaterm", "false");
    System.setProperty("lambda-torsions", "false");

    forceFieldEnergy = new ForceFieldEnergy(active);
    forcefieldEnergy.setPrintOnFailure(false, false);
    e = minimize(eps);
    energy();
    if (size > 1){
        structureFile = new File("postOSRW."+ String.format("%d",world.rank())+"." + structureFile.getName());
    } else{
        structureFile = new File("postOSRW."+structureFile.getName());
    }
    saveAsPDB(structureFile);
}


if (runOSRW && size > 1){

    DoubleBuf receiveBuffer = DoubleBuf.buffer(energyArray);
    if (!(world.rank() == 0)){
        world.receive(world.rank()-1,receiveBuffer);
        energyArray[world.rank()-1] = receiveBuffer.get(world.rank()-1);
    }
    if(!loopBuildError){
        energyArray[world.rank()] = active.getPotentialEnergy().getTotalEnergy();
    } else {
        energyArray[world.rank()] = Double.MAX_VALUE;
    }
    if (world.rank() < world.size()-1){
        DoubleBuf sendBuffer = DoubleBuf.buffer(energyArray);
        world.send(world.rank()+1,sendBuffer);
    }
    world.barrier();
    if(world.rank() == world.size()-1){
        for(int i = 0; i < world.size(); i++) {
            String resultFileName = "Loop.txt";
            File rankAndEnergyFile = new File(resultFileName);
            BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(rankAndEnergyFile, true));

            String rankString = Integer.toString(i);
            String energyString = Double.toString(energyArray[i]);

            bufferedWriter.write( rankString + ":" + energyString);
            bufferedWriter.newLine();
            bufferedWriter.flush();
        }
    }
    saveAsPDB(structureFile);
}

if (runRotamer){

    for (int i = 0; i <= atoms.length; i++) {
        Atom ai = atoms[i - 1];
        ai.setActive(true);
        ai.setUse(true);
    }

    forceFieldEnergy = new ForceFieldEnergy(active);
    forcefieldEnergy.setPrintOnFailure(false, false);
    Polymer[] polymers = active.getChains();
    ArrayList<Residue> fullResidueList = polymers[0].getResidues();
    ArrayList<Residue> residuesToRO = new ArrayList<>();

    //Rotamer Optimization inclusion list building (grab built residues)
    for (int i = 0; i < fullResidueList.size(); i++) {
        Residue r = fullResidueList[i];
        if (r.getBackboneAtoms().get(0).getBuilt()) {
            residuesToRO.add(fullResidueList[i]);
        }
    }

    int startResID = residuesToRO.get(0).getResidueNumber();
    int finalResID = residuesToRO.get(residuesToRO.size() - 1).getResidueNumber();

    //Find best loop generated by multiple walkers
    if (runOSRW && size > 1){
        world.barrier();
        int bestRank;

        if (world.rank() == 0){
            int[] loopRanks = new int[size];
            double[] loopEnergies = new double[size];

            double lowestEnergy = Double.MAX_VALUE;
            BufferedReader reader = new BufferedReader(new FileReader("Loop.txt"));
            String line = null;
            int i = 0;
            while((line = reader.readLine()) != null){
                String[] lineData;
                lineData = line.split(":");
                //lineData[0] contains rank information (see Loop.txt)
                //lineData[1] contains energy information (see Loop.txt)
                if(Double.parseDouble(lineData[1]) < lowestEnergy){
                    bestRank = Integer.parseInt(lineData[0]);
                    lowestEnergy = Double.parseDouble(lineData[1]);
                }
                i++;
            }
        }
        world.barrier();
        IntegerBuf broadcastBuf = IntegerBuf.buffer(bestRank);
        world.broadcast(0,broadcastBuf);
        bestRank = broadcastBuf.get(0);

        if(world.rank()==bestRank){
            energy();
        }

        active.destroy();
        File bestRankDirectory = new File(structureFile.getParentFile().getParent() + File.separator + Integer.toString(bestRank));
        File bestStructureFile = new File(bestRankDirectory.getPath() + File.separator + structureFile.getName());

        //buildLoops=false needed to avoid error in PDBFilter because saved loops do not have dbref/seqres information
        System.setProperty("buildLoops", "false");
        logger.info(String.format("Path to best loop " + bestRankDirectory.getPath() + File.separator + bestStructureFile.getName()));
        open(bestRankDirectory.getPath() + File.separator + bestStructureFile.getName());
        structureFile = bestStructureFile;
    }

    forceFieldEnergy = new ForceFieldEnergy(active);
    forcefieldEnergy.setPrintOnFailure(false, false);
    boolean threeBodyTerm = false;
    RotamerOptimization rotamerOptimization;

    energy();

    logger.info(String.format("Rotamer Optimization"));
    rotamerOptimization = new RotamerOptimization(active, forceFieldEnergy, null);

    rotamerOptimization.setThreeBodyEnergy(threeBodyTerm);
    rLib.setUseOrigCoordsRotamer(true);
    //rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards);
    rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson);

    //Rotamer Optimization inclusion list building (grab residues within 7A of the built loop)
    boolean expandList = true
    double expansionDistance = 7.0;
    rLib.setUseOrigCoordsRotamer(true);

    if (expandList) {
        // Do a sliding-window rotamer optimization on loop window with a radius-inclusion criterion.
        rotamerOptimization.setForcedResidues(startResID,finalResID);
        rotamerOptimization.setWindowSize(1);
        rotamerOptimization.setDistanceCutoff(expansionDistance);
    }
    rotamerOptimization.setResidues(startResID, finalResID);
    residuesToRO = rotamerOptimization.getResidues();

    RotamerLibrary.measureRotamers(residuesToRO, false);
    rotamerOptimization.optimize(RotamerOptimization.Algorithm.SLIDING_WINDOW);
}
saveAsPDB(structureFile);
