/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */
package ffx.algorithms;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.Rotamer;
import ffx.potential.RotamerLibrary;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;

/**
 * Optimize protein side-chain conformations using rotamers.
 *
 * @author Michael J. Schnieders
 */
public class RotamerOptimization implements Terminatable {

    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());

    public enum Algorithm {

        INDEPENDENT, GLOBAL, SLIDING_WINDOW
    }
    MolecularAssembly molecularAssembly;
    AlgorithmListener algorithmListener;
    ForceFieldEnergy potential;
    Polymer[] polymers;
    int startResID = -1;
    int finalResID = -1;
    Algorithm algorithm = null;
    boolean terminate = false;
    boolean done = true;

    public RotamerOptimization(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener) {
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
    }

    public RotamerOptimization(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener,
            int startResID, int finalResID, Algorithm algorithm) {
        this(molecularAssembly, algorithmListener);
        this.startResID = startResID;
        this.finalResID = finalResID;
        this.algorithm = algorithm;
    }

    public double optimize() {
        double e = Double.MAX_VALUE;
        done = false;
        terminate = false;
        if (finalResID >= startResID) {
            potential = (ForceFieldEnergy) molecularAssembly.getPotentialEnergy();
            polymers = molecularAssembly.getChains();
            switch (algorithm) {
                case INDEPENDENT:
                    e = independent();
                    break;
                case GLOBAL:
                    e = global();
                    break;
                case SLIDING_WINDOW:
                    logger.warning(" The sliding window rotamer optimization has not been implemented yet.");
                    break;
            }
        }
        terminate = false;
        done = true;
        return e;
    }

    public double optimize(int startResID, int finalResID, Algorithm algorithm) {
        this.startResID = startResID;
        this.finalResID = finalResID;
        this.algorithm = algorithm;
        return optimize();
    }

    private double independent() {
        double e = Double.MAX_VALUE;
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            logger.info(String.format(" Optimizing %s side-chain.", residue));
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            if (rotamers == null) {
                continue;
            }
            e = potential.energy(false, true);
            int bestRotamer = -1;
            for (int j = 0; j < rotamers.length; j++) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
                if (algorithmListener != null) {
                    algorithmListener.algorithmUpdate(molecularAssembly);
                }
                double newE = potential.energy(false, true);
                if (newE < e) {
                    bestRotamer = j;
                }
            }
            if (bestRotamer > -1) {
                Rotamer rotamer = rotamers[bestRotamer];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
            }
            if (algorithmListener != null) {
                algorithmListener.algorithmUpdate(molecularAssembly);
            }
            if (terminate) {
                logger.info(String.format("\n Terminating after residue %s.\n", residue));
                break;
            }
        }
        return e;
    }

    private double global() {
        ArrayList<Residue> residues = new ArrayList<Residue>();
        int permutations = 1;
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            residues.add(residue);
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            if (rotamers != null) {
                permutations *= rotamers.length;
            }
        }
        logger.info(String.format(" Number of permutations: %d.", permutations));
        ArrayList<Integer> optimum = new ArrayList<Integer>();
        double e = RotamerLibrary.rotamerOptimization(molecularAssembly, residues, Double.MAX_VALUE, optimum);
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymers[0].getResidue(i);
            AminoAcid3 name = AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(name);
            int j = optimum.remove(0);
            if (rotamers != null) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(name, residue, rotamer);
            }
        }
        return e;
    }

    private double slidingWindow() {
        double e = Double.MAX_VALUE;
        // TODO
        return e;
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, " Exception terminating rotamer optimization.\n", e);
                }
            }
        }
    }
}
