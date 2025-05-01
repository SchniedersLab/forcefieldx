// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.potential.nonbonded.pme;

import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;

import java.util.logging.Logger;

public class PMETimings {

  private static final Logger logger = Logger.getLogger(PMETimings.class.getName());

  public final long[] realSpaceEnergyTime;
  public final long[] realSpaceSCFTime;
  /** Timing variables. */
  private final int numThreads;

  public long realSpaceEnergyTotal, realSpaceSCFTotalTime;
  public long bornRadiiTotal, gkEnergyTotal;

  public PMETimings(int numThreads) {
    this.numThreads = numThreads;
    realSpaceEnergyTime = new long[numThreads];
    realSpaceSCFTime = new long[numThreads];
  }

  public void init() {
    for (int i = 0; i < numThreads; i++) {
      realSpaceEnergyTime[i] = 0;
      realSpaceSCFTime[i] = 0;
    }
    realSpaceEnergyTotal = 0;
    realSpaceSCFTotalTime = 0;
    bornRadiiTotal = 0;
    gkEnergyTotal = 0;
  }

  public void printRealSpaceTimings(int maxThreads,
                                    PermanentFieldRegion permanentFieldRegion,
                                    RealSpaceEnergyRegion realSpaceEnergyRegion) {
    long realSpacePermTotal = permanentFieldRegion.getRealSpacePermTime();

    double total = (realSpacePermTotal + realSpaceSCFTotalTime + realSpaceEnergyTotal) * NS2SEC;
    logger.info(format("\n Real Space: %7.4f (sec)", total));
    logger.info("           Electric Field");
    logger.info(" Thread    Direct  SCF     Energy     Counts");
    long minPerm = Long.MAX_VALUE;
    long maxPerm = 0;
    long minSCF = Long.MAX_VALUE;
    long maxSCF = 0;
    long minEnergy = Long.MAX_VALUE;
    long maxEnergy = 0;
    int minCount = Integer.MAX_VALUE;
    int maxCount = Integer.MIN_VALUE;

    for (int i = 0; i < maxThreads; i++) {
      int count = realSpaceEnergyRegion.getCount(i);
      long realSpacePermTime = permanentFieldRegion.getInitTime(i) + permanentFieldRegion.getPermTime(i);
      logger.info(format("    %3d   %7.4f %7.4f %7.4f %10d",
          i, realSpacePermTime * NS2SEC, realSpaceSCFTime[i] * NS2SEC,
          realSpaceEnergyTime[i] * NS2SEC, count));
      minPerm = min(realSpacePermTime, minPerm);
      maxPerm = max(realSpacePermTime, maxPerm);
      minSCF = min(realSpaceSCFTime[i], minSCF);
      maxSCF = max(realSpaceSCFTime[i], maxSCF);
      minEnergy = min(realSpaceEnergyTime[i], minEnergy);
      maxEnergy = max(realSpaceEnergyTime[i], maxEnergy);
      minCount = min(count, minCount);
      maxCount = max(count, maxCount);
    }
    logger.info(format(" Min      %7.4f %7.4f %7.4f %10d",
        minPerm * NS2SEC, minSCF * NS2SEC, minEnergy * NS2SEC, minCount));
    logger.info(format(" Max      %7.4f %7.4f %7.4f %10d",
        maxPerm * NS2SEC, maxSCF * NS2SEC, maxEnergy * NS2SEC, maxCount));
    logger.info(format(" Delta    %7.4f %7.4f %7.4f %10d",
        (maxPerm - minPerm) * NS2SEC,
        (maxSCF - minSCF) * NS2SEC,
        (maxEnergy - minEnergy) * NS2SEC,
        (maxCount - minCount)));
    logger.info(format(" Actual   %7.4f %7.4f %7.4f %10d",
        realSpacePermTotal * NS2SEC,
        realSpaceSCFTotalTime * NS2SEC,
        realSpaceEnergyTotal * NS2SEC,
        realSpaceEnergyRegion.getInteractions()));
  }
}
