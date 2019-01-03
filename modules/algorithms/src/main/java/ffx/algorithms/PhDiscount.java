/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;

import ffx.algorithms.MolecularDynamics.MonteCarloNotification;
import ffx.potential.AssemblyState;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Torsion;
import ffx.potential.extended.ExtConstants;
import ffx.potential.extended.ExtUtils;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig;
import ffx.potential.extended.TitrationUtils.Snapshots;
import ffx.potential.extended.TitrationUtils.Titration;
import ffx.potential.extended.TitrationUtils.TitrationConfig;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.utils.SystemTemperatureException;
import static ffx.potential.extended.ExtConstants.ns2sec;
import static ffx.potential.extended.TitrationUtils.propagateInactiveResidues;

/**
 * <p>PhDiscount class.</p>
 *
 * @author S. LuCore
 */
public class PhDiscount implements MonteCarloListener {
    private static final Logger logger = Logger.getLogger(PhDiscount.class.getName());

    private final TitrationConfig config;
    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final MolecularDynamics molDyn;
    private final ffx.algorithms.thermostats.Thermostat thermostat;
    private final String originalFilename;
    private List<MultiResidue> titratingMultiResidues = new ArrayList<>();
    private final RotamerLibrary library = new RotamerLibrary(RotamerLibrary.ProteinLibrary.Richardson, false);

    private final Random rng = new Random();
    private int snapshotIndex = 0;

    /* Fractional Protonation Mode: Extended Variables */
    private final ExtendedSystem esvSystem;
    private final ExtendedSystemConfig esvConfig;
    private Atom[] extendedAtoms;
    private int nAtomsExt;
    private int numESVs;

    /* Original parameters to MolecularDynamics constructor; necessary for relaunch. */
    private final double dt;
    private final double printInterval;
    private final double saveInterval;
    private final boolean initVelocities;
    private final String fileType;
    private final double writeRestartInterval;
    private final File dynLoader;

    /* Debug */
    private static final String keyPrefixes[] = new String[]{"phmd", "esv", "md", "ffe", "sys", "db", "sdl"};

    /**
     * Construct a "Discrete-Continuous" Monte-Carlo titration engine.
     * For traditional discrete titration, use Protonate.
     * For traditional continuous titration, run mdesv for populations.
     *
     * @param mola                 the molecular assembly
     * @param esvSystem            a {@link ffx.potential.extended.ExtendedSystem} object.
     * @param molDyn               a {@link ffx.algorithms.MolecularDynamics} object.
     * @param timeStep             a double.
     * @param printInterval        a double.
     * @param saveInterval         a double.
     * @param initVelocities       a boolean.
     * @param fileType             a {@link java.lang.String} object.
     * @param writeRestartInterval a double.
     * @param dyn                  a {@link java.io.File} object.
     */
    // java.lang.Double, java.lang.Integer, java.lang.Integer, java.lang.Boolean, java.lang.String, java.lang.Integer, null)
    public PhDiscount(MolecularAssembly mola, ExtendedSystem esvSystem, MolecularDynamics molDyn,
                      double timeStep, double printInterval, double saveInterval,
                      boolean initVelocities, String fileType, double writeRestartInterval, File dyn) {
        this.config = new TitrationConfig();
        this.mola = mola;
        this.ffe = mola.getPotentialEnergy();
        this.esvSystem = esvSystem;
        this.esvConfig = esvSystem.config;
        this.numESVs = esvSystem.size();
        this.molDyn = molDyn;
        this.dt = timeStep;
        this.printInterval = printInterval;
        this.saveInterval = saveInterval;
        this.initVelocities = initVelocities;
        this.fileType = fileType;
        this.writeRestartInterval = writeRestartInterval;
        this.dynLoader = dyn;
        this.thermostat = molDyn.getThermostat();
        this.originalFilename = FilenameUtils.removeExtension(mola.getFile().getAbsolutePath()) + "_dyn.pdb";
        SystemFilter.setVersioning(SystemFilter.Versioning.PREFIX_ABSOLUTE);

        config.print();
        ExtUtils.printConfigSet("All Config:", System.getProperties(), keyPrefixes);
        logger.info(format(" Running DISCOuNT-pH dynamics @ system pH %.2f\n", esvSystem.getConstantPh()));
        ffe.reInit();
        molDyn.reInit();
        molDyn.setMonteCarloListener(this, MonteCarloNotification.EACH_STEP);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Update the position of background atoms at each step, lest bonded term energies become inconsistent.
     */
    @Override
    public boolean mcUpdate(double temperature) {
        esvSystem.setTemperature(temperature);
        propagateInactiveResidues(titratingMultiResidues, true);
        return true;
    }

    /**
     * Top-level driver for DISCOUNT pHMD.
     *
     * @param physicalSteps    a int.
     * @param pH               a double.
     * @param temperature      a double.
     * @param attemptFrequency a int.
     * @param stepsPerMove     a int.
     */
    public void dynamic(int physicalSteps, double pH, double temperature,
                        int attemptFrequency, int stepsPerMove) {
        int movesAttempted = 0, movesAccepted = 0;

        /* Init step. */
        molDyn.dynamic(1, dt, printInterval, saveInterval, temperature,
                initVelocities, fileType, writeRestartInterval, dynLoader);

        /* Launch first round of MD with traditional call to dynamic(). */
        logger.info(format(" Launching DISCOUNT pHMD for %d physical steps.", physicalSteps));
        logger.info(format(" (Round %d) Launching fixed-protonation dynamics for %d steps.",
                movesAttempted + 1, attemptFrequency));

//        molDyn.dynamic(attemptFrequency, dt, printInterval, saveInterval,
//                temperature, initVelocities, fileType, writeRestartInterval, dynLoader);
        molDyn.redynamic(attemptFrequency, temperature);
        int stepsTaken = attemptFrequency;

        while (stepsTaken < physicalSteps) {
            /* Attempt Monte-Carlo through fractional-protonation dynamics. */
            logger.info(format(" (Attempt %d) Launching fractional-protonation dynamics for %d steps.",
                    ++movesAttempted, attemptFrequency));
            if (tryContinuousTitrationMove(stepsPerMove, temperature)) {
                movesAccepted++;
            }

            /* Return to fixed-protonation MD with redynamic() to avoid initialization. */
            logger.info(format(" (Round %d) Launching fixed-protonation dynamics for %d steps.",
                    movesAttempted + 1, attemptFrequency));
            if (stepsTaken + attemptFrequency < physicalSteps) {
                molDyn.redynamic(attemptFrequency, temperature);
                stepsTaken += attemptFrequency;
            } else {
                molDyn.redynamic(physicalSteps - stepsTaken, temperature);
                stepsTaken = physicalSteps;
            }
        }

        /* Profit{!,?,.} */
        logger.info(format(" (Summary) DISCOUNT pHMD completed %d physical steps and %d moves, of which %d were accepted.",
                physicalSteps, movesAttempted, movesAccepted));
    }

    /**
     * Attempt to print sources of catastrophic system heating.
     */
    private void crashDump(RuntimeException error) {
        writeSnapshot(".meltdown-");
        ffe.energy(false, true);
        mola.getDescendants(BondedTerm.class).stream()
                .filter(BondedTerm::isExtendedSystemMember)
                .forEach(term -> {
                    try {
                        ((Bond) term).log();
                    } catch (Exception ex) {
                    }
                    try {
                        ((Angle) term).log();
                    } catch (Exception ex) {
                    }
                    try {
                        ((Torsion) term).log();
                    } catch (Exception ex) {
                    }
                });
        if (ffe.getVanDerWaalsEnergy() > 1000) {
            extendedAtoms = esvSystem.getExtendedAtoms();
            nAtomsExt = esvSystem.getExtendedAtoms().length;
            for (int i = 0; i < nAtomsExt; i++) {
                Atom ai = extendedAtoms[i];
                for (int j = 0; j < nAtomsExt; j++) {
                    Atom aj = extendedAtoms[j];
                    if (!esvSystem.isExtended(i) && !esvSystem.isExtended(j)) {
                        continue;
                    }
                    if (ai == aj || ai.getBond(aj) != null) {
                        continue;
                    }
                    double dist = FastMath.sqrt(
                            FastMath.pow((aj.getX() - ai.getX()), 2) +
                                    FastMath.pow((aj.getY() - ai.getY()), 2) +
                                    FastMath.pow((aj.getZ() - ai.getZ()), 2));
                    if (dist < 0.8 * (aj.getVDWR() + ai.getVDWR())) {
                        logger.warning(String.format("Close vdW contact for atoms: \n   %s\n   %s", aj, ai));
                    }
                }
            }
        }
        throw error;
    }

    /**
     * Run continuous titration MD in implicit solvent as a Monte Carlo move.
     */
    private boolean tryContinuousTitrationMove(int titrationDuration, double targetTemperature) {
        long startTime = System.nanoTime();
        if (thermostat.getCurrentTemperature() > config.meltdownTemperature) {
            crashDump(new SystemTemperatureException(thermostat.getCurrentTemperature()));
        }
        propagateInactiveResidues(titratingMultiResidues, true);

        // Save the current state of the molecularAssembly. Specifically,
        //      Atom coordinates and MultiResidue states : AssemblyState
        //      Position, Velocity, Acceleration, etc    : DynamicsState (via MD::storeState)
        AssemblyState assemblyState = new AssemblyState(mola);
        molDyn.storeState();
        writeSnapshot(".pre-store");

        /* Assign initial titration lambdas. */
        switch (config.seedDistribution) {
            default:
            case FLAT:
                /* All lambdas have equal probability. */
                for (int i = 0; i < numESVs; i++) {
                    esvSystem.setLambda(i, rng.nextDouble());
                }
                break;
            case BETA:
                /* Draw from a mathematical distribution. */
                throw new UnsupportedOperationException();
            case BOLTZMANN:
                /* Draw from a Boltzmann distribution (via Rosenbluth sets). */
                throw new UnsupportedOperationException();
            case DIRAC_CURRENT:
                /* Assign only zero or unity representing current state. */
                for (int i = 0; i < numESVs; i++) {
                    Residue active = titratingMultiResidues.get(i).getActive();
                    double current = (active.getAminoAcid3() == Titration.lookup(active).protForm)
                            ? 1.0 : 0.0;
                    esvSystem.setLambda(i, current);
                }
                break;
            case DIRAC_POINTFIVE:
                /* Assign half-occupancy to all titratable protons. */
                for (int i = 0; i < numESVs; i++) {
                    esvSystem.setLambda(i, 0.5);
                }
                break;
        }

        /*
         * (1) Take pre-move energy.
         * (2) Hook the ExtendedSystem up to MolecularDynamics.
         * (3) Launch dynamics for nSteps = Monte-Carlo Frequency.
         * (4) Note that no callbacks to mcUpdate() occur during this period.
         * (5) Floor/ceil to discretize hydrogen occupancies.
         * (6) Take post-move energy and test on the combined Metropolis criterion.
         */
        final double Uo = currentTotalEnergy();

        // TODO: Need to ensure that we fully protonate the protein before entering continuous-protonation space.
        ffe.attachExtendedSystem(esvSystem);
        molDyn.attachExtendedSystem(esvSystem, 10);
        final double Uo_prime = currentTotalEnergy();
        logger.info(format(" %-40s %-s", "Trying continuous titration move.",
                format("Uo,Uo': %16.8f, %16.8f", Uo, Uo_prime)));

        molDyn.redynamic(titrationDuration, targetTemperature);
        final double Un_prime = currentTotalEnergy();
        for (int i = 0; i < esvSystem.size(); i++) {
            esvSystem.setLambda(i, Math.rint(esvSystem.getLambda(i)));
        }
        molDyn.detachExtendedSystem();
        ffe.detachExtendedSystem();
        final double Un = currentTotalEnergy();
        logger.info(format(" %-40s %-30s", "Move finished; detaching esvSystem.",
                format("Un',Un: %16.8f, %16.8f", Un_prime, Un)));

        /* Calculate acceptance probability from detailed balance equation. */
        final double beta = 1 / (ExtConstants.Boltzmann * thermostat.getCurrentTemperature());
        final double dgDiscrete = Un - Uo;
        final double dgContinuous = Un_prime - Uo_prime;
        final double crit = FastMath.exp(-beta * (dgDiscrete - dgContinuous));
        final double rand = rng.nextDouble();
        logger.info(format("   Un,Un',Uo',Uo:   %10.5f %10.5f %10.5f %10.5f", Un, Un_prime, Uo_prime, Uo));
        logger.info(format("   Crit,Rng:        %10.5f %10.5f", crit, rand));
        long took = System.nanoTime() - startTime;
        if (rand <= crit) {
            logger.info(format(" %-40s %-s", "Monte-Carlo accepted.",
                    format("Wallclock: %8.3f sec", took * ns2sec)));
            writeSnapshot(".post-acpt");
            return true;
        } else {
            logger.info(format(" %-40s %-s", "Monte-Carlo denied; reverting state.",
                    format("Wallclock: %8.3f sec", took * ns2sec)));
            writeSnapshot(".post-deny");
            assemblyState.revertState();
            ffe.reInit();
            try {
                molDyn.revertState();
            } catch (Exception ex) {
                Logger.getLogger(PhDiscount.class.getName()).log(Level.SEVERE, null, ex);
            }
            writeSnapshot(".post-rvrt");
            return false;
        }
    }

    private double currentTotalEnergy() {
        double x[] = new double[ffe.getNumberOfVariables() * 3];
        ffe.getCoordinates(x);
        ffe.energy(x);
        return ffe.getTotalEnergy();
    }

    private void writeSnapshot(String extension) {
        String basename = FilenameUtils.removeExtension(originalFilename);
        String filename;
        if (config.snapshots == Snapshots.INTERLEAVED) {
            if (basename.contains("_dyn")) {
                filename = basename.replace("_dyn", format("_dyn_%d.pdb", ++snapshotIndex));
            } else {
                filename = FilenameUtils.removeExtension(basename)
                        + format("_dyn_%d.pdb", ++snapshotIndex);
            }
        } else {
            if (!extension.startsWith(".")) {
                extension = "." + extension;
            }
            filename = basename + format("_%d", ++snapshotIndex) + extension;
        }
        File file = new File(filename);
        PDBFilter writer = new PDBFilter(file, mola, null, null);
        writer.writeFile(file, false);
    }

}
