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
package ffx.crystal;

import static ffx.crystal.SpaceGroupInfo.sohnckeGroup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The Spacegroup class defines the symmetry of a crystal. There are 230 distinct space groups in
 * three dimensions.
 *
 * @author Michael J. Schnieders
 * @see <ul>
 *     <li><a href="http://it.iucr.org/Ab/" target="_blank"> International Tables for
 *     Crystallography Volume A: Space-group symmetry </a>
 *     <li><a href="http://legacy.ccp4.ac.uk/html/symmetry.html" target="_blank"> CCP4 Symlib </a>
 *     </ul>
 * @since 1.0
 */
public class SpaceGroup {

  /** Space group number. */
  public final int number;
  /** Number of primitive symmetry equivalents. */
  public final int numPrimitiveSymEquiv;
  /** Space group name. */
  public final String shortName;
  /**
   * Point group name. There are 32 distinct points groups, or crystal classes in three dimensions.
   */
  public final String pointGroupName;
  /** Crystal system. */
  public final CrystalSystem crystalSystem;
  /** Laue group */
  public final LaueSystem laueSystem;
  /** Space group name under the PDB convention. */
  public final String pdbName;
  /** A List of SymOp instances. */
  public final List<SymOp> symOps;
  /** True for a Sohncke group (non-enantiogenic). */
  public final boolean respectsChirality;
  /** Real space ASU limit operators. */
  public final ASULimit[] asuLimitOperators;
  /** Lattice system. */
  public final LatticeSystem latticeSystem;
  /** Number of symmetry equivalents. */
  private final int numSymEquiv;
  /** Real space ASU limit values. */
  private final double[] asuLimits;

  /**
   * Immutable SpaceGroup instances are made available only through the factory method so this
   * constructor is private.
   *
   * @param number Space group number.
   * @param numSymEquiv Number of symmetry equivalents.
   * @param numPrimitiveSymEquiv Number of primitive symmetry equivalents.
   * @param shortName Short PDB name.
   * @param pointGroupName Point group name.
   * @param pdbName PDB space group name.
   * @param crystalSystem Crystal system.
   * @param laueSystem Laue System.
   * @param symOps Symmetry operators.
   * @param asuLimits Assymetric unit limit.
   * @param asuLimitOperators ASULimit instance.
   * @since 1.0
   */
  protected SpaceGroup(
      int number,
      int numSymEquiv,
      int numPrimitiveSymEquiv,
      String shortName,
      String pointGroupName,
      String pdbName,
      CrystalSystem crystalSystem,
      LatticeSystem latticeSystem,
      LaueSystem laueSystem,
      ASULimit[] asuLimitOperators,
      double[] asuLimits,
      SymOp... symOps) {
    this.number = number;
    this.numSymEquiv = numSymEquiv;
    this.numPrimitiveSymEquiv = numPrimitiveSymEquiv;
    this.shortName = shortName;
    this.pointGroupName = pointGroupName;
    this.crystalSystem = crystalSystem;
    this.latticeSystem = latticeSystem;
    this.laueSystem = laueSystem;
    this.asuLimitOperators = asuLimitOperators;
    this.asuLimits = asuLimits;
    this.pdbName = pdbName;
    this.respectsChirality = sohnckeGroup(number);
    this.symOps = new ArrayList<>(Arrays.asList(symOps));

    // ToDo: Crystal systems are subdivided into crystal classes. This info needs to be added to
    // each space group.
  }

  /**
   * Return the number of symmetry operators.
   *
   * @return the number of symmetry operators.
   * @since 1.0
   */
  public int getNumberOfSymOps() {
    return symOps.size();
  }

  /**
   * Return the ith symmetry operator.
   *
   * @param i the symmetry operator number.
   * @return the SymOp
   * @since 1.0
   */
  public SymOp getSymOp(int i) {
    return symOps.get(i);
  }

  /**
   * Check if the space group maintains chirality.
   *
   * @return Return true if chirality is respected.
   */
  public boolean respectsChirality() {
    return respectsChirality;
  }

}
