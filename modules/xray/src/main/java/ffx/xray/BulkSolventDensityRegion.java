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
package ffx.xray;

import static java.util.Arrays.fill;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;

/**
 * This class implements a spatial decomposition based on partitioning a grid into octants. The
 * over-ridden "selectAtoms" method selects atoms that are not in the asymmetric unit, but are
 * within the supplied cutoff radius.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BulkSolventDensityRegion extends SpatialDensityRegion {

  private final BulkSolventList bulkSolventList;

  /**
   * Constructor for BulkSolventDensityRegion.
   *
   * @param gX a int.
   * @param gY a int.
   * @param gZ a int.
   * @param grid an array of double.
   * @param basisSize a int.
   * @param nSymm a int.
   * @param minWork a int.
   * @param threadCount a int.
   * @param crystal a {@link ffx.crystal.Crystal} object.
   * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
   * @param coordinates an array of double.
   * @param cutoff a double.
   * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
   */
  public BulkSolventDensityRegion(
      int gX,
      int gY,
      int gZ,
      double[] grid,
      int basisSize,
      int nSymm,
      int minWork,
      int threadCount,
      Crystal crystal,
      Atom[] atoms,
      double[][][] coordinates,
      double cutoff,
      ParallelTeam parallelTeam) {
    super(gX, gY, gZ, grid, basisSize, nSymm, minWork, threadCount, crystal, atoms, coordinates);

    // Asymmetric unit atoms never selected by this class.
    fill(select[0], false);
    bulkSolventList = new BulkSolventList(crystal, atoms, cutoff, parallelTeam);
  }

  /** {@inheritDoc} */
  @Override
  public void run() {
    int ti = getThreadIndex();
    int actualWork1 = actualWork - 1;
    SpatialDensityLoop loop = spatialDensityLoop[ti];
    try {
      execute(0, actualWork1, loop.setOctant(0));
      // Fractional chunks along the C-axis.
      if (nC > 1) {
        execute(0, actualWork1, loop.setOctant(1));
        // Fractional chunks along the B-axis.
        if (nB > 1) {
          execute(0, actualWork1, loop.setOctant(2));
          execute(0, actualWork1, loop.setOctant(3));
          // Fractional chunks along the A-axis.
          if (nA > 1) {
            execute(0, actualWork1, loop.setOctant(4));
            execute(0, actualWork1, loop.setOctant(5));
            execute(0, actualWork1, loop.setOctant(6));
            execute(0, actualWork1, loop.setOctant(7));
          }
        }
      }
    } catch (Exception e) {
      logger.severe(e.toString());
    }
  }

  /** {@inheritDoc} */
  @Override
  public void selectAtoms() {
    bulkSolventList.buildList(coordinates, select, false);
  }
}
