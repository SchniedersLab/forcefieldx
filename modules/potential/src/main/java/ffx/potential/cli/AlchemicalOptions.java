/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential.cli;

import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

import ffx.potential.utils.PotentialsFunctions;
import picocli.CommandLine.Option;

/**
 * Alchemical options shared by scripts that do alchemy on a single topology and use the Pico CLI.
 */
public class AlchemicalOptions {

    /**
     * The logger for this class.
     */
    public static final Logger logger = Logger.getLogger(AlchemicalOptions.class.getName());

    // A regular expression used to parse ranges of atoms.
    public static final Pattern rangeregex = Pattern.compile("([0-9]+)-?([0-9]+)?");

    /**
     * -l or --lambda sets the initial lambda value.
     */
    @Option(names = {"-l", "--lambda"}, paramLabel = "-1",
            description = "Initial lambda value.")
    double initialLambda = -1.0;

    /**
     * -s1 or --start1 defines the first softcored atom for the first topology.
     */
    @Option(names = {"--s1", "--start1"}, paramLabel = "0",
            description = "Starting ligand atom for 1st topology.")
    int s1 = 0;

    /**
     * -f1 or --final1 defines the last softcored atom for the first topology.
     */
    @Option(names = {"--f1", "--final1"}, paramLabel = "-1",
            description = "Final ligand atom for the 1st topology.")
    int f1 = -1;

    /**
     * --la1 or -ligAtoms1 allows for multiple ranges and/or singletons of ligand atoms in the first topology, separated by periods.
     */
    @Option(names = {"--la1", "--ligAtoms1"},
            description = "Period-separated ranges of 1st topology ligand atoms (e.g. 40-50.72-83).")
    String ligAt1 = null;

    /**
     * -es1 or --noElecStart1 defines the first atom of the first topology to have no electrostatics.
     */
    @Option(names = {"--es1", "--noElecStart1"}, paramLabel = "1",
            description = "Starting no-electrostatics atom for 1st topology.")
    int es1 = 1;

    /**
     * -ef1 or --noElecFinal1 defines the last atom of the first topology to have no electrostatics.
     */
    @Option(names = {"--ef1", "--noElecFinal1"}, paramLabel = "-1",
            description = "Final no-electrostatics atom for 1st topology.")
    int ef1 = -1;

    /**
     * -as or --activeStart starts an active set of atoms for single-topology lambda gradients.
     */
    @Option(names = {"--as", "--activeStart"}, paramLabel = "1",
            description = "Starting active atom (single-topology only).")
    int actStart = 1;

    /**
     * -af or --activeFinal ends an active set of atoms for single-topology lambda gradients.
     */
    @Option(names = {"--af", "--activeFinal"}, paramLabel = "-1",
            description = "Final active atom (single-topology only).")
    int actFinal = -1;

    /**
     * Set active atoms for a MolecularAssembly.
     *
     * @param molecularAssembly
     */
    public void setActiveAtoms(MolecularAssembly molecularAssembly) {
        Atom[] atoms = molecularAssembly.getAtomArray();
        if (actFinal > 0) {
            // Apply active atom selection
            int nAtoms = atoms.length;
            if (actFinal > actStart && actStart > 0 && actFinal <= nAtoms) {
                // Make all atoms inactive.
                for (int i = 0; i <= nAtoms; i++) {
                    Atom ai = atoms[i - 1];
                    ai.setActive(false);
                }
                // Make requested atoms active.
                for (int i = actStart; i <= actFinal; i++) {
                    Atom ai = atoms[i - 1];
                    ai.setActive(true);
                }
            }
        }
    }

    /**
     * Set the alchemical atoms for this molecularAssembly.
     *
     * @param topology
     */
    public void setFirstSystemAlchemistry(MolecularAssembly topology) {
        setAlchemicalAtoms(topology, s1, f1, ligAt1);
    }

    /**
     * Sets the alchemical atoms for a MolecularAssembly.
     *
     * @param assembly Assembly to which the atoms belong.
     * @param start First atom to set lambda for.
     * @param fin Last atom to set lambda for.
     * @param ligAt Additional ranges of atoms to set lambda for.
     */
    public static void setAlchemicalAtoms(MolecularAssembly assembly, int start, int fin, String ligAt) {
        Atom[] atoms = assembly.getAtomArray();
        if (start > 0) {
            for (int i = start; i <= fin; i++) {
                Atom ai = atoms[i-1];
                ai.setApplyLambda(true);
                // ai.print();
            }
        }
        if (ligAt != null) {
            String[] ranges = ligAt.split("\\.");
            for (String range : ranges) {
                Matcher m = rangeregex.matcher(range);
                if (m.find()) {
                    int rangeStart = Integer.parseInt(m.group(1));
                    int rangeEnd = (m.group(2) != null) ? Integer.parseInt(m.group(2)) : rangeStart;
                    if (rangeStart > rangeEnd) {
                        logger.severe(String.format(" Range %s was invalid; start was greater than end", range));
                    }
                    // Don't need to worry about negative numbers; rangeregex just won't match.
                    for (int i = rangeStart; i <= rangeEnd; i++) {
                        Atom ai = atoms[i - 1];
                        ai.setApplyLambda(true);
                        ai.print();
                    }
                } else {
                    logger.warning(String.format(" Could not recognize %s as a valid range; skipping", range));
                }
            }
        }
    }

    /**
     * Set uncharged atoms for this molecularAssembly.
     *
     * @param topology
     */
    public void setFirstSystemUnchargedAtoms(MolecularAssembly topology) {
        setUnchargedAtoms(topology, es1, ef1);
    }

    /**
     * Set uncharged atoms for a MolecularAssembly.
     * @param assembly Assembly to decharge on.
     * @param eStart First atom to decharge.
     * @param eEnd Last atom to decharge.
     */
    public static void setUnchargedAtoms(MolecularAssembly assembly, int eStart, int eEnd) {
        Atom atoms[] = assembly.getAtomArray();
        // Apply the no electrostatics atom selection
        int noElecStart = eStart;
        noElecStart = (noElecStart < 1) ? 1 : noElecStart;

        int noElecStop = eEnd;
        noElecStop = (noElecStop > atoms.length) ? atoms.length : noElecStop;

        for (int i = noElecStart; i <= noElecStop; i++) {
            Atom ai = atoms[i - 1];
            ai.setElectrostatics(false);
            // ai.print();
        }
    }

    /**
     * If any softcore Atoms have been detected.
     *
     * @return Presence of softcore Atoms.
     */
    public boolean hasSoftcore() {
        return ((ligAt1 != null && ligAt1.length() > 0) || s1 > 0);
    }

    /**
     * Gets the initial value of lambda.
     *
     * @return Initial lambda.
     */
    public double getInitialLambda() {
        return getInitialLambda(true);
    }

    /**
     * Gets the initial value of lambda.
     *
     * @param quiet If true, do not warn about lambda not being in the range 0-1.
     * @return Initial lambda.
     */
    public double getInitialLambda(boolean quiet) {
        Level toLog = quiet ? Level.OFF : Level.WARNING;
        if (initialLambda < 0.0 || initialLambda > 1.0) {
            logger.log(toLog, String.format(" Initial alchemical lambda reset " +
                    "to 0.0 from %8.4g; must be between 0 and 1!", initialLambda));
            initialLambda = 0.0;
        }
        return initialLambda;
    }

    /**
     * Opens a File to a MolecularAssembly for alchemistry.
     *
     * @param potentialFunctions A utility object for opening Files into MolecularAssemblys.
     * @param topOptions TopologyOptions in case a dual-topology or greater is to be used.
     * @param threadsPer Number of threads to be used for this MolecularAssembly.
     * @param toOpen The name of the File to be opened.
     * @param topNum The index of this topology.
     * @return The processed MolecularAssembly.
     */
    public MolecularAssembly openFile(PotentialsFunctions potentialFunctions, Optional<TopologyOptions> topOptions, int threadsPer, String toOpen, int topNum) {
        potentialFunctions.openAll(toOpen, threadsPer);
        MolecularAssembly mola = potentialFunctions.getActiveAssembly();
        return processFile(topOptions, mola, topNum);
    }

    /**
     * Performs processing on a MolecularAssembly for alchemistry.
     *
     * @param topOptions TopologyOptions in case a dual-topology or greater is to be used.
     * @param mola The MolecularAssembly to be processed.
     * @param topNum The index of this topology, 0-indexed.
     * @return The processed MolecularAssembly.
     */
    public MolecularAssembly processFile(Optional<TopologyOptions> topOptions, MolecularAssembly mola, int topNum) {

        int remainder = (topNum % 2) + 1;
        switch (remainder) {
            case 1:
                setFirstSystemAlchemistry(mola);
                setFirstSystemUnchargedAtoms(mola);
                break;
            case 2:
                if (!topOptions.isPresent()) {
                    throw new IllegalArgumentException(" For >= 2 systems, topOptions must not be empty!");
                }
                TopologyOptions topology = topOptions.get();
                topology.setSecondSystemAlchemistry(mola);
                topology.setSecondSystemUnchargedAtoms(mola);
                break;
        }

        // Turn off checks for overlapping atoms, which is expected for lambda=0.
        ForceFieldEnergy energy = mola.getPotentialEnergy();
        energy.getCrystal().setSpecialPositionCutoff(0.0);

        return mola;
    }
}
