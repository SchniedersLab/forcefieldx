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

// MOLECULAR & STOCHASTIC DYNAMICS

// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.DiscountPh;
import ffx.algorithms.Protonate.DynamicsLauncher;
import ffx.algorithms.Integrator.Integrators;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.extended.TitrationESV;
import ffx.potential.extended.TitrationESV.TitrationUtils;

// Number of molecular dynamics steps
int nSteps = 1000000;

// Time step in femtoseconds.
double dt = 1.0;

// Frequency to print out thermodynamics information in picoseconds.
double printInterval = 0.01;

// Frequency to save out coordinates in picoseconds.
double saveInterval = 0.1;

// Temperature in degrees Kelvin.
double temperature = 298.15;

// Thermostats [ ADIABATIC, BERENDSEN, BUSSI ]
Thermostats thermostat = null;

// Integrators [ BEEMAN, RESPA, STOCHASTIC, VELOCITYVERLET]
Integrators integrator = null;

// Reset velocities (ignored if a restart file is given)
boolean initVelocities = true;

// Interval to write out restart file (psec)
double restartFrequency = 1000;

// File type of snapshots.
String fileType = "PDB";

// Monte-Carlo step frequencies for titration and rotamer moves.
int titrationFrequency = 10;
int titrationDuration = 1000;
int rotamerMoveRatio = 0;

// Simulation pH
double pH = 7.4;

// Titrating residue list.
List<String> rlTokens = new ArrayList<>();

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc discount [options] <filename>');
cli.h(longOpt:'help', 'Print this message.');
//cli.b(longOpt:'thermostat', args:1, argName:'Berendsen', 'Thermostat: [Adiabatic / Berendsen / Bussi]');
cli.d(longOpt:'dt', args:1, argName:'1.0', 'Time discretization (fsec).');
//cli.i(longOpt:'integrate', args:1, argName:'Beeman', 'Integrator: [Beeman / RESPA / Stochastic / VELOCITYVERLET]');
cli.l(longOpt:'log', args:1, argName:'0.01', 'Interval to log thermodyanamics (psec).');
cli.n(longOpt:'steps', args:1, argName:'1000000', 'Number of molecular dynamics steps.');
cli.p(longOpt:'polarization', args:1, argName:'Mutual', 'Polarization: [None / Direct / Mutual]');
cli.t(longOpt:'temperature', args:1, argName:'298.15', 'Temperature in degrees Kelvin.');
cli.w(longOpt:'save', args:1, argName:'0.1', 'Interval to write out coordinates (psec).');
cli.s(longOpt:'restart', args:1, argName:'0.1', 'Interval to write out restart file (psec).');
cli.f(longOpt:'file', args:1, argName:'PDB', 'Choose file type to write to [PDB/XYZ]');
//cli.ra(longOpt:'resAll', 'Titrate all residues.');
cli.rl(longOpt:'resList', args:1, 'Titrate a list of residues (eg A4.A8.B2.B34)');
//cli.rn(longOpt:'resName', args:1, 'Titrate a list of residue names (eg "LYS,TYR,HIS")');
//cli.rw(longOpt:'resWindow', args:1, 'Titrate all residues with intrinsic pKa within [arg] units of simulation pH.');
cli.pH(longOpt:'pH', args:1, argName:'7.4', 'Constant simulation pH.');
cli.mc(longOpt:'titrationFrequency', args:1, argName:'100', 'Number of steps between Monte-Carlo proton attempts.')
cli.mcd(longOpt:'titrationDuration', args:1, argName:'100', 'Number of steps for which to run continuous proton dynamics during MC move.');
cli.mcr(longOpt:'rotamerMoveRatio', args:1, argName:'0', 'Number of steps between Monte-Carlo rotamer attempts.')
//cli.tt(longOpt:'titrateTermini', args:1, argName:'false', 'Titrate amino acid chain ends.');
def options = cli.parse(args);

if (options.h) {
    return cli.usage();
}

// Required Options
if (!options.rl) {
    logger.warning("Must list titrating residues (-rl).");
    return;
} else {
    def tok = (options.rl).tokenize('.');
    for (String t : tok) {
        rlTokens.add(t);
    }
}

// Suggested Options
if (!options.pH) {
    logger.warning("No solution pH specified; 7.4 assumed.");
    pH = 7.4;
} else {
    pH = Double.parseDouble(options.pH);
}

if (options.rw) {
    window = Double.parseDouble(options.rw);
}

if (options.mc) {
    titrationFrequency = Integer.parseInt(options.mc);
}

if (options.mcd) {
    titrationDuration = Integer.parseInt(options.mcd);
}

if (options.mcr) {
    rotamerMoveRatio = Integer.parseInt(options.mcr);
}

// Load the number of molecular dynamics steps.
if (options.n) {
    nSteps = Integer.parseInt(options.n);
}

// Write dyn interval in picoseconds
if (options.s) {
    restartFrequency = Double.parseDouble(options.s);
}

//
if (options.f) {
    fileType = options.f.toUpperCase();
}
// Load the time steps in femtoseconds.
if (options.d) {
    dt = Double.parseDouble(options.d);
}

// Report interval in picoseconds.
if (options.l) {
    printInterval = Double.parseDouble(options.l);
}

// Write snapshot interval in picoseconds.
if (options.w) {
    saveInterval = Double.parseDouble(options.w);
}

// Temperature in degrees Kelvin.
if (options.t) {
    temperature = Double.parseDouble(options.t);
}

if (options.p) {
    System.setProperty("polarization", options.p);
}

// Thermostat.
if (options.b) {
    try {
        thermostat = Thermostats.valueOf(options.b.toUpperCase());
    } catch (Exception e) {
        thermostat = null;
    }
}

if (options.mcmd) {
    mcRunTime = Integer.parseInt(options.mcmd);
}

// Integrator.
if (options.i) {
    try {
        integrator = Integrators.valueOf(options.i.toUpperCase());
    } catch (Exception e) {
        integrator = null;
    }
}

System.setProperty("forcefield","AMOEBA_PROTEIN_2013");
System.setProperty("pme-qi","true");
System.setProperty("polarization", "NONE");
System.setProperty("ligand-vapor-elec","false");
System.setProperty("no-ligand-condensed-scf","false");
System.setProperty("mpoleterm","false");
//System.setProperty("mpoleterm","true");

/*  Available options:  Key                         Values (default first)
    System.getProperty("cphmd-seedMode")            CURRENT_VALUES|RANDOM|HALF_LAMBDA
    System.getProperty("cphmd-override")            NONE|ACCEPT|REJECT
    System.getProperty("cphmd-snapshotsType")       SEPARATE|INTERLEAVED|NONE
    System.getProperty("cphmd-histidineMode")       HIE_ONLY|HID_ONLY|SINGLE|DOUBLE
    System.getProperty("cphmd-debugLog")            int     ()
    System.getProperty("cphmd-referenceOverride")   double  ()
    System.getProperty("cphmd-tempMonitor")         double  (6000.0)
    System.getProperty("cphmd-logTimings")          boolean
    System.getProperty("cphmd-termini")             boolean
    System.getProperty("cphmd-zeroReferences")      boolean
*/

List<String> arguments = options.arguments();
String modelfilename = null;
if (arguments != null && arguments.size() > 0) {
    // Read in command line.
    modelfilename = arguments.get(0);
    open(modelfilename);
} else if (active == null) {
    return cli.usage();
} else {
    modelfilename = active.getFile();
}

logger.info("\n Running molecular dynamics on " + modelfilename);

// Restart File
File dyn = new File(FilenameUtils.removeExtension(modelfilename) + ".dyn");
if (!dyn.exists()) {
    dyn = null;
}

// Create the ExtendedSystem.
MolecularAssembly mola = (MolecularAssembly) active;
ExtendedSystem esvSystem = new ExtendedSystem(mola);

// Find the requested residues and create TitrationESVs from them.
List<ExtendedVariable> esvList = new ArrayList<>();
Polymer[] polymers = active.getChains();
for (String token : rlTokens) {    
    char chainID = token.charAt(0);
    int resNum = Integer.parseInt(token.substring(1));    
    Residue target = null;
    for (Polymer p : polymers) {
        char pid = p.getChainID().charValue();
        if (pid == chainID) {
            for (Residue res : p.getResidues()) {
                if (res.getResidueNumber() == resNum) {
                    target = res;
                    break;
                }
            }
            if (target != null) {
                break;
            }
        }
    }
    if (target == null) {
        logger.severe("Couldn't find target residue " + token);
    }
    
    MultiResidue titrating = TitrationUtils.titrationFactory(mola, target);
    TitrationESV esv = new TitrationESV(titrating, 7.4, biasMag);
    esvSystem.addVariable(esv);
}

// Attach populated esvSystem to the potential.
ForceFieldEnergy ffe = mola.getPotentialEnergy();
ffe.attachExtendedSystem(esvSystem);

// Create the MolecularDynamics.
MolecularDynamics molDyn = new MolecularDynamics(active, active.getPotentialEnergy(), active.getProperties(), sh, thermostat, integrator);
molDyn.setFileType(fileType);
molDyn.setRestartFrequency(restartFrequency);

// Create the DISCOuNT object, linking it to the MD.
DiscountPh cphmd = new DiscountPh(mola, esvSystem, molDyn,
    dt, printInterval, saveInterval, initVelocities, 
    fileType, restartFrequency, dyn);

// Launch dynamics through the DISCOuNT controller.
cphmd.dynamic(nSteps, pH, temperature, titrationFrequency, titrationDuration, rotamerMoveRatio);
