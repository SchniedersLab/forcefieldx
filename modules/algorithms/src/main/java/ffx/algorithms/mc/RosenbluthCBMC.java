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
package ffx.algorithms.mc;

import java.io.File;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import ffx.algorithms.MonteCarloListener;
import ffx.algorithms.thermostats.Thermostat;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parsers.PDBFilter;
import static ffx.algorithms.mc.BoltzmannMC.BOLTZMANN;

/**
 * Conformational Biased Monte Carlo (applied to ALL torsions of a peptide
 * side-chain). As described by Frenkel/Smit in "Understanding Molecular
 * Simulation" Chapters 13.2,13.3 This uses the "conformational biasing" method
 * to select whole rotamer transformations that are frequently accepted.
 *
 * @author S. LuCore
 */
public class RosenbluthCBMC implements MonteCarloListener {

    private static final Logger logger = Logger.getLogger(RosenbluthCBMC.class.getName());

    private final MolecularAssembly mola;
    private final ForceFieldEnergy ffe;
    private final Thermostat thermostat;
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
    private int numMovesAccepted = 0;
    /**
     * Writes PDBs of each trial set and original/proposed configurations.
     */
    private boolean writeSnapshots = false;

    /**
     * RRMC constructor.
     *
     * @param targets        Residues to undergo RRMC.
     * @param mcFrequency    Number of MD steps between RRMC proposals.
     * @param trialSetSize   Larger values cost more but increase acceptance.
     * @param mola           a {@link ffx.potential.MolecularAssembly} object.
     * @param ffe            a {@link ffx.potential.ForceFieldEnergy} object.
     * @param thermostat     a {@link ffx.algorithms.thermostats.Thermostat} object.
     * @param writeSnapshots a boolean.
     */
    public RosenbluthCBMC(MolecularAssembly mola, ForceFieldEnergy ffe, Thermostat thermostat,
                          List<Residue> targets, int mcFrequency, int trialSetSize, boolean writeSnapshots) {
        this.targets = targets;
        this.mcFrequency = mcFrequency;
        this.trialSetSize = trialSetSize;
        this.mola = mola;
        this.ffe = ffe;
        this.thermostat = thermostat;
        this.writeSnapshots = writeSnapshots;
        for (int i = targets.size() - 1; i >= 0; i--) {
            AminoAcid3 name = AminoAcid3.valueOf(targets.get(i).getName());
            if (name == AminoAcid3.GLY || name == AminoAcid3.PRO || name == AminoAcid3.ALA) {
                targets.remove(i);
            }
        }
//        targets.removeIf(p -> AminoAcid3.valueOf(p.getName()) == AminoAcid3.ALA
//                || AminoAcid3.valueOf(p.getName()) == AminoAcid3.GLY
//                || AminoAcid3.valueOf(p.getName()) == AminoAcid3.PRO);
        if (targets.size() < 1) {
            logger.severe("Empty target list for CMBC.");
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
        double temperature;
        if (thermostat != null) {
            temperature = thermostat.getCurrentTemperature();
        } else {
            temperature = 298.15;
        }
        double beta = 1.0 / (BOLTZMANN * temperature);

        // Select a target residue.
        int which = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(which);
        RosenbluthChiAllMove cbmcMove = new RosenbluthChiAllMove(
                mola, target, trialSetSize, ffe, temperature,
                writeSnapshots, numMovesProposed, true);
        if (cbmcMove.getMode() == RosenbluthChiAllMove.MODE.CHEAP) {
            if (cbmcMove.wasAccepted()) {
                numMovesAccepted++;
            }
            return cbmcMove.wasAccepted();
        }
        double Wn = cbmcMove.getWn();
        double Wo = cbmcMove.getWo();
        double criterion = Math.min(1, Wn / Wo);
        double rng = ThreadLocalRandom.current().nextDouble();
        logger.info(String.format("    rng:    %5.2f", rng));
        if (rng < criterion) {
            cbmcMove.move();
            numMovesAccepted++;
            logger.info(String.format(" Accepted!  Energy: %.4f\n", cbmcMove.finalEnergy));
            accepted = true;
            write();
        } else {
            logger.info(String.format(" Denied.\n"));
            accepted = false;
        }

        return accepted;
    }

    private PDBFilter writer;

    private void write() {
        if (writer == null) {
            writer = new PDBFilter(mola.getFile(), mola, null, null);
        }
        String filename = FilenameUtils.removeExtension(mola.getFile().toString());
        filename = mola.getFile().getAbsolutePath();
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
        double temperature;
        if (thermostat != null) {
            temperature = thermostat.getCurrentTemperature();
        } else {
            temperature = 298.15;
        }
        double beta = 1.0 / (BOLTZMANN * temperature);
        int which = ThreadLocalRandom.current().nextInt(targets.size());
        Residue target = targets.get(which);
        RosenbluthChiAllMove cbmcMove = new RosenbluthChiAllMove(
                mola, target, -1, ffe, temperature,
                false, numMovesProposed, true);
        return cbmcMove.wasAccepted();
    }

}
