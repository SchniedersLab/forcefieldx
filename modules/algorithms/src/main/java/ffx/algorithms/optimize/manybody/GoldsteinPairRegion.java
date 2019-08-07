//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.optimize.manybody;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import static java.lang.Double.isFinite;
import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;

import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.Residue;

public class GoldsteinPairRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(GoldsteinPairRegion.class.getName());

    private RotamerOptimization rotamerOptimization;
    private Residue[] residues;
    private int i, riA, rjC;
    private int j, riB, rjD;
    private int[] possK;
    private int nK;
    private GoldsteinRotamerPairLoop[] goldsteinRotamerPairLoop;
    private SharedDouble sharedSumOverK = new SharedDouble();
    private ArrayList<Residue> blockedResidues;

    public GoldsteinPairRegion(int nThreads) {
        goldsteinRotamerPairLoop = new GoldsteinRotamerPairLoop[nThreads];
    }

    /**
     * Initializes a ParallelRegion to attempt the elimination
     * of riA,rjC by riB,rjD.
     *
     * @param residues         The residue array.
     * @param i                First residue of the pair.
     * @param riA              First member of the pair to attempt eliminating.
     * @param riB              First member of the pair to try eliminating by.
     * @param j                Second residue of the pair.
     * @param rjC              Second member of the pair to attempt eliminating.
     * @param rjD              Second member of the pair to try eliminating by.
     * @param bidiResNeighbors All interaction partners of a Residue, including prior residues
     */
    public void init(Residue[] residues, int i, int riA, int riB, int j, int rjC, int rjD,
                     int[][] bidiResNeighbors, RotamerOptimization rotamerOptimization) {
        this.residues = residues;
        this.i = i;
        this.riA = riA;
        this.riB = riB;
        this.j = j;
        this.rjC = rjC;
        this.rjD = rjD;
        this.rotamerOptimization = rotamerOptimization;
        int[] nI = bidiResNeighbors[i];
        int[] nJ = bidiResNeighbors[j];
        IntStream kStream = IntStream.concat(Arrays.stream(nI), Arrays.stream(nJ));
        possK = kStream.distinct().
                filter(k -> (k != i && k != j)).
                sorted().
                toArray();
        nK = possK.length;
    }

    public double getSumOverK() {
        return sharedSumOverK.get();
    }

    public ArrayList<Residue> getMissedResidues() {
        return blockedResidues;
    }

    public void start() {
        sharedSumOverK.set(0.0);
        blockedResidues = new ArrayList<>();
    }

    public void finish() {
        for (GoldsteinRotamerPairLoop rotamerPairLoop : goldsteinRotamerPairLoop) {
            blockedResidues.addAll(rotamerPairLoop.blockedResidues);
        }
    }

    @Override
    public void run() {
        int threadID = getThreadIndex();
        if (goldsteinRotamerPairLoop[threadID] == null) {
            goldsteinRotamerPairLoop[threadID] = new GoldsteinRotamerPairLoop();
        }
        try {
            execute(0, nK - 1, goldsteinRotamerPairLoop[threadID]);
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception in GoldsteinPairRegion.", e);
        }
    }

    private class GoldsteinRotamerPairLoop extends IntegerForLoop {

        double sumOverK;
        ArrayList<Residue> blockedResidues;

        @Override
        public void start() {
            sumOverK = 0.0;
            blockedResidues = new ArrayList<>();
        }

        @Override
        public void finish() {
            sharedSumOverK.addAndGet(sumOverK);
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            if (blockedResidues.isEmpty()) {
                double locSumOverK = rotamerOptimization.goldsteinPairSumOverK(
                        residues, lb, ub, i, riA, riB, j, rjC, rjD, blockedResidues, possK);
                // Should be redundant checks.
                if (isFinite(locSumOverK) && blockedResidues.isEmpty()) {
                    sumOverK += locSumOverK;
                } else {
                    sumOverK = 0;
                }
            } else {
                rotamerOptimization.logIfMaster(format(" Skipping %d to %d because we cannot eliminate", lb, ub), Level.FINE);
            }
        }
    }
}
