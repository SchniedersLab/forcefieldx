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
package ffx.algorithms.mc;

import java.io.File;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.algorithms.MonteCarloListener;
import ffx.algorithms.thermostats.Thermostat;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parsers.PDBFilter;

/**
 * Conformational Biased Monte Carlo (applied to ALL torsions of a peptide
 * side-chain).
 * <p>
 * This method is described by Frenkel/Smit in "Understanding Molecular
 * Simulation" Chapters 13.2,13.3 This uses the "conformational biasing" method
 * to select whole rotamer transformations that are frequently accepted.
 *
 * @author Stephen D. LuCore
 */
public class RosenbluthCBMC implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(RosenbluthCBMC.class.getName());

    /**
     * The MolecularAssembly to operate on.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * The ForceFieldEnergy in use.
     */
    private final ForceFieldEnergy forceFieldEnergy;
    /**
     * The temperature in use.
     */
    private final double temperature;
    /**
     * At each move, one of these residues will be chosen as the target.
     */
    private final List<Residue> targets;
    /**
     * Number of mcUpdate() calls (e.g. MD steps) between move proposals.
     */
    private final int mcFrequency;
    /**
     * Keeps track of calls to mcUpdate (e.g. MD steps).
     */
    private int steps = 0;
    /**
     * Size of the trial sets, k.
     */
    private final int trialSetSize;
    /**
     * Counters for proposed and accepted moves.
     */
    private int numMovesProposed = 0;
    /**
     * Writes PDBs of each trial set and original/proposed configurations.
     */
    private boolean writeSnapshots;
    /**
     * PDBFilter to write out result.
     */
    private PDBFilter writer = null;

    /**
     * RRMC constructor.
     *
     * @param targets           Residues to undergo RRMC.
     * @param mcFrequency       Number of MD steps between RRMC proposals.
     * @param trialSetSize      Larger values cost more but increase acceptance.
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param ffe               a {@link ffx.potential.ForceFieldEnergy} object.
     * @param thermostat        a {@link ffx.algorithms.thermostats.Thermostat} object.
     * @param writeSnapshots    a boolean.
     */
    public RosenbluthCBMC(MolecularAssembly molecularAssembly, ForceFieldEnergy ffe, Thermostat thermostat,
                          List<Residue> targets, int mcFrequency, int trialSetSize, boolean writeSnapshots) {
        this.targets = targets;
        this.mcFrequency = mcFrequency;
        this.trialSetSize = trialSetSize;
        this.molecularAssembly = molecularAssembly;
        this.forceFieldEnergy = ffe;
        this.writeSnapshots = writeSnapshots;

        if (thermostat != null) {
            temperature = thermostat.getTargetTemperature();
        } else {
            temperature = 298.15;
        }

        for (int i = targets.size() - 1; i >= 0; i--) {
            AminoAcid3 name = AminoAcid3.valueOf(targets.get(i).getName());
            if (name == AminoAcid3.GLY || name == AminoAcid3.PRO || name == AminoAcid3.ALA) {
                targets.remove(i);
            }
        }
        if (targets.size() < 1) {
            logger.severe(" Empty target list for CMBC.");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean mcUpdate(double temperature) {
        steps++;
        if (steps % mcFrequency == 0) {
            return cbmcStep();
        }
        return false;
    }

    /**
     * <p>cbmcStep.</p>
     *
     * @return a boolean.
     */
    public boolean cbmcStep() {
        numMovesProposed++;
        boolean accepted;

        // Select a target residue.
        int index = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(index);
        RosenbluthChiAllMove cbmcMove = new RosenbluthChiAllMove(
                molecularAssembly, target, trialSetSize, forceFieldEnergy, temperature,
                writeSnapshots, numMovesProposed, true);
        if (cbmcMove.getMode() == RosenbluthChiAllMove.MODE.CHEAP) {
            return cbmcMove.wasAccepted();
        }
        double Wn = cbmcMove.getWn();
        double Wo = cbmcMove.getWo();
        double criterion = min(1, Wn / Wo);
        double rng = ThreadLocalRandom.current().nextDouble();
        logger.info(format("    rng:    %5.2f", rng));
        if (rng < criterion) {
            cbmcMove.move();
            logger.info(format(" Accepted!  Energy: %.4f\n", cbmcMove.finalEnergy));
            accepted = true;
            write();
        } else {
            logger.info(" Denied.\n");
            accepted = false;
        }

        return accepted;
    }

    /**
     * Write out a PDB file.
     */
    private void write() {
        if (writer == null) {
            writer = new PDBFilter(molecularAssembly.getFile(), molecularAssembly, null, null);
        }
        String filename = molecularAssembly.getFile().getAbsolutePath();
        if (!filename.contains("_mc")) {
            filename = FilenameUtils.removeExtension(filename) + "_mc.pdb";
        }
        File file = new File(filename);
        writer.writeFile(file, false);
    }

    /**
     * <p>controlStep.</p>
     *
     * @return a boolean.
     */
    public boolean controlStep() {
        int index = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(index);
        RosenbluthChiAllMove cbmcMove = new RosenbluthChiAllMove(
                molecularAssembly, target, -1, forceFieldEnergy, temperature,
                false, numMovesProposed, true);
        return cbmcMove.wasAccepted();
    }

}
