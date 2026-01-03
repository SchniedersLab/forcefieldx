// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.numerics.clustering;

public final class SampleClusterData {

  public static final double[][] DISTANCES = new double[][] { { 0, 1, 9, 7, 11, 14 },
      { 1, 0, 4, 3, 8, 10 }, { 9, 4, 0, 9, 2, 8 }, { 7, 3, 9, 0, 6, 13 }, { 11, 8, 2, 6, 0, 10 },
      { 14, 10, 8, 13, 10, 0 } };
  public static final String[] NAMES = new String[] { "O1", "O2", "O3", "O4", "O5", "O6" };
  public static final double[] WEIGHTS = new double[] { 1, 2, 3, 4, 5, 6 };

  public static final String[] NAMES_WITH_DUPLICATE = new String[] { "O2", "O2", "O3", "O4", "O5",
      "O6" };

}
