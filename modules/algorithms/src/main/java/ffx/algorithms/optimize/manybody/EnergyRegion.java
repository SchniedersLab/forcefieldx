// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;
import ffx.potential.bonded.Residue;
import java.util.Arrays;

public class EnergyRegion extends ParallelRegion {

  private EnergyExpansion eE;
  /** Flag to control use of 3-body terms. */
  private boolean threeBodyTerm;

  private SharedDouble self;
  private SharedDouble twoBody;
  private SharedDouble threeBody;
  private EnergyLoop[] energyLoops;
  private Residue[] residues;
  private int[] rotamers;
  private int nResidues;

  public EnergyRegion(int nThreads) {
    self = new SharedDouble();
    twoBody = new SharedDouble();
    threeBody = new SharedDouble();
    energyLoops = new EnergyLoop[nThreads];
  }

  public double getSelf() {
    return self.get();
  }

  public double getThreeBody() {
    return threeBody.get();
  }

  public double getTwoBody() {
    return twoBody.get();
  }

  public void init(EnergyExpansion eE, Residue[] residues, int[] rotamers, boolean threeBodyTerm) {
    this.eE = eE;
    this.rotamers = rotamers;
    this.nResidues = residues.length;
    this.residues = Arrays.copyOf(residues, nResidues);
    this.threeBodyTerm = threeBodyTerm;
  }

  @Override
  public void run() throws Exception {
    int threadID = getThreadIndex();
    if (energyLoops[threadID] == null) {
      energyLoops[threadID] = new EnergyLoop();
    }
    execute(0, nResidues - 1, energyLoops[threadID]);
  }

  public void start() {
    self.set(0.0);
    twoBody.set(0.0);
    threeBody.set(0.0);
  }

  private class EnergyLoop extends IntegerForLoop {

    private double selfSum;
    private double pairSum;
    private double threeBodySum;

    @Override
    public void finish() {
      self.addAndGet(selfSum);
      twoBody.addAndGet(pairSum);
      threeBody.addAndGet(threeBodySum);
    }

    @Override
    public void run(int lb, int ub) {
      for (int a = lb; a <= ub; a++) {
        int ai = rotamers[a];
        selfSum += eE.getSelf(a, ai);
        for (int b = a + 1; b < nResidues; b++) {
          int bi = rotamers[b];
          pairSum += eE.get2Body(a, ai, b, bi);
          if (threeBodyTerm) {
            for (int c = b + 1; c < nResidues; c++) {
              int ci = rotamers[c];
              threeBodySum += eE.get3Body(residues, a, ai, b, bi, c, ci);
            }
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.dynamic();
    }

    @Override
    public void start() {
      selfSum = 0.0;
      pairSum = 0.0;
      threeBodySum = 0.0;
    }
  }
}
