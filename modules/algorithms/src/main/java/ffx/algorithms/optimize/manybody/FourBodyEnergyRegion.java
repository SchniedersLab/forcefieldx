// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.algorithms.optimize.manybody;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;

import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerRegion;
import ffx.algorithms.optimize.RotamerOptimization;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import java.util.stream.DoubleStream;

/** Compute 4-Body energies. This code is experimental. */
public class FourBodyEnergyRegion extends WorkerRegion {

  private static final Logger logger = Logger.getLogger(FourBodyEnergyRegion.class.getName());
  private final Residue[] residues;
  private final RotamerOptimization rO;
  private final DistanceMatrix dM;
  private final EnergyExpansion eE;
  private final EliminatedRotamers eR;
  /**
   * A list of all residues being optimized. Note that Box and Window optimizations operate on
   * subsets of this list.
   */
  private final List<Residue> allResiduesList;
  /** Map of 3-body energy values to compute. */
  private final Map<Integer, Integer[]> fourBodyEnergyMap;
  /**
   * If a pair of residues have two atoms closer together than the superposition threshold, the
   * energy is set to NaN.
   */
  private final double superpositionThreshold;

  private Set<Integer> keySet;

  public FourBodyEnergyRegion(
      RotamerOptimization rotamerOptimization,
      DistanceMatrix dM,
      EnergyExpansion eE,
      EliminatedRotamers eR,
      Residue[] residues,
      List<Residue> allResiduesList,
      double superpositionThreshold) {
    this.rO = rotamerOptimization;
    this.dM = dM;
    this.eE = eE;
    this.eR = eR;
    this.residues = residues;
    this.allResiduesList = allResiduesList;
    this.superpositionThreshold = superpositionThreshold;

    this.fourBodyEnergyMap = eE.getFourBodyEnergyMap();
    logger.info(format(" Running quads: %d jobs.", fourBodyEnergyMap.size()));
  }

  @Override
  public void run() throws Exception {
    if (!keySet.isEmpty()) {
      execute(0, keySet.size() - 1, new QuadsEnergyLoop());
    }
  }

  @Override
  public void start() {
    keySet = fourBodyEnergyMap.keySet();
  }

  private class QuadsEnergyLoop extends WorkerIntegerForLoop {

    @Override
    public void run(int lb, int ub) {
      for (int key = lb; key <= ub; key++) {
        long time = -System.nanoTime();
        if (!fourBodyEnergyMap.containsKey(key)) {
          continue;
        }

        Integer[] job = fourBodyEnergyMap.get(key);
        int i = job[0];
        int ri = job[1];
        int j = job[2];
        int rj = job[3];
        int k = job[4];
        int rk = job[5];
        int l = job[6];
        int rl = job[7];

        if (eR.check(i, ri)
            || eR.check(j, rj)
            || eR.check(k, rk)
            || eR.check(l, rl)
            || eR.check(i, ri, j, rj)
            || eR.check(i, ri, k, rk)
            || eR.check(i, ri, l, rl)
            || eR.check(j, rj, k, rk)
            || eR.check(j, rj, l, rl)
            || eR.check(k, rk, l, rl)) {
          // Not implemented: 3-body or 4-body checks.
          continue;
        }

        Residue resi = residues[i];
        Residue resj = residues[j];
        Residue resk = residues[k];
        Residue resl = residues[l];

        int indexI = allResiduesList.indexOf(residues[i]);
        int indexJ = allResiduesList.indexOf(residues[j]);
        int indexK = allResiduesList.indexOf(residues[k]);
        int indexL = allResiduesList.indexOf(residues[l]);

        double rawDist = dM.getRawNBodyDistance(indexI, ri, indexJ, rj, indexK, rk, indexL, rl);
        double dIJ = dM.checkDistMatrix(indexI, ri, indexJ, rj);
        double dIK = dM.checkDistMatrix(indexI, ri, indexK, rk);
        double dIL = dM.checkDistMatrix(indexI, ri, indexL, rl);
        double dJK = dM.checkDistMatrix(indexJ, rj, indexK, rk);
        double dJL = dM.checkDistMatrix(indexJ, rj, indexL, rl);
        double dKL = dM.checkDistMatrix(indexK, rk, indexL, rl);

        double minDist = DoubleStream.of(dIJ, dIK, dIL, dJK, dJL, dKL).min().getAsDouble();

        String distString = "     large";
        if (rawDist < Double.MAX_VALUE) {
          distString = format("%10.3f", rawDist);
        }

        double resDist = dM.get4BodyResidueDistance(indexI, ri, indexJ, rj, indexK, rk, indexL, rl);
        String resDistString = "     large";
        if (resDist < Double.MAX_VALUE) {
          resDistString = format("%5.3f", resDist);
        }

        double fourBodyEnergy = 0.0;
        if (minDist < superpositionThreshold) {
          fourBodyEnergy = Double.NaN;
          logger.info(
              format(
                  " Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d:   set to NaN at %13.6f Ang (%s Ang by residue)  < %5.3f Ang.",
                  residues[i],
                  ri,
                  residues[j].toFormattedString(false, true),
                  rj,
                  residues[k].toFormattedString(false, true),
                  rk,
                  residues[l].toFormattedString(false, true),
                  rl,
                  minDist,
                  resDistString,
                  superpositionThreshold));
        } else if (dM.checkQuadDistThreshold(indexI, ri, indexJ, rj, indexK, rk, indexL, rl)) {
          // Set the 4-body energy to 0.0 for separation distances larger than the 4-body cutoff.
          fourBodyEnergy = 0.0;
          time += System.nanoTime();
          logger.info(
              format(
                  " Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                  resi.toFormattedString(false, true),
                  ri,
                  resj.toFormattedString(false, true),
                  rj,
                  resk.toFormattedString(false, true),
                  rk,
                  resl.toFormattedString(false, true),
                  rl,
                  rO.formatEnergy(fourBodyEnergy),
                  distString,
                  resDistString));
        } else {
          try {
            fourBodyEnergy = eE.compute4BodyEnergy(residues, i, ri, j, rj, k, rk, l, rl);
            time += System.nanoTime();
            logger.info(
                format(
                    " Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).",
                    resi.toFormattedString(false, true),
                    ri,
                    resj.toFormattedString(false, true),
                    rj,
                    resk.toFormattedString(false, true),
                    rk,
                    resl.toFormattedString(false, true),
                    rl,
                    rO.formatEnergy(fourBodyEnergy),
                    distString,
                    resDistString));
            if (abs(fourBodyEnergy) > 1.0) {
              StringBuilder sb = new StringBuilder();
              sb.append(
                  format(
                      " Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: %s at %s Ang (%s Ang by residue).\n",
                      resi.toFormattedString(false, true),
                      ri,
                      resj.toFormattedString(false, true),
                      rj,
                      resk.toFormattedString(false, true),
                      rk,
                      resl.toFormattedString(false, true),
                      rl,
                      rO.formatEnergy(fourBodyEnergy),
                      distString,
                      resDistString));
              sb.append(format("   Explain: (ref %d) \n", key));
              sb.append(
                  format("     Self %3d %3d:                  %.3f\n", i, ri, eE.getSelf(i, ri)));
              sb.append(
                  format("     Self %3d %3d:                  %.3f\n", j, rj, eE.getSelf(j, rj)));
              sb.append(
                  format("     Self %3d %3d:                  %.3f\n", k, rk, eE.getSelf(k, rk)));
              sb.append(
                  format("     Self %3d %3d:                  %.3f\n", l, rl, eE.getSelf(l, rl)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      i, ri, j, rj, eE.get2Body(i, ri, j, rj)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      i, ri, k, rk, eE.get2Body(i, ri, k, rk)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      i, ri, l, rl, eE.get2Body(i, ri, l, rl)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      j, rj, k, rk, eE.get2Body(j, rj, k, rk)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      j, rj, l, rl, eE.get2Body(j, rj, l, rl)));
              sb.append(
                  format(
                      "     Pair %3d %3d %3d %3d:          %.3f\n",
                      k, rk, l, rl, eE.get2Body(k, rk, l, rl)));
              sb.append(
                  format(
                      "     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n",
                      i, ri, j, rj, k, rk, eE.get3Body(residues, i, ri, j, rj, k, rk)));
              sb.append(
                  format(
                      "     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n",
                      i, ri, j, rj, l, rl, eE.get3Body(residues, i, ri, j, rj, l, rl)));
              sb.append(
                  format(
                      "     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n",
                      i, ri, k, rk, l, rl, eE.get3Body(residues, i, ri, k, rk, l, rl)));
              sb.append(
                  format(
                      "     Tri  %3d %3d %3d %3d %3d %3d:  %.3f\n",
                      j, rj, k, rk, l, rl, eE.get3Body(residues, j, rj, k, rk, l, rl)));
              sb.append(
                  format("     backbone:                      %.3f\n", rO.getBackboneEnergy()));
              sb.append(format("     quadEnergy:                 %.3f\n", fourBodyEnergy));
              sb.append("     --s--\n");
              sb.append("     Active residues:\n");
              for (Residue residue : residues) {
                if (residue.getSideChainAtoms().get(0).getUse()) {
                  sb.append(format("       %s\n", residue.toString()));
                }
              }
              sb.append("     --f--\n");
              logger.info(sb.toString());
            }
          } catch (ArithmeticException ex) {
            fourBodyEnergy = Double.NaN;
            time += System.nanoTime();
            logger.info(
                format(
                    " Quad %8s %-2d, %8s %-2d, %8s %-2d, %8s %-2d: NaN at %s Ang (%s Ang by residue).",
                    resi.toFormattedString(false, true),
                    ri,
                    resj.toFormattedString(false, true),
                    rj,
                    resk.toFormattedString(false, true),
                    rk,
                    resl.toFormattedString(false, true),
                    rl,
                    distString,
                    resDistString));
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
