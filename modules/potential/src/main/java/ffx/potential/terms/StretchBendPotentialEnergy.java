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
package ffx.potential.terms;

import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.StretchBend;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;

/**
 * Stretch-Bend potential energy term using {@link ffx.potential.bonded.StretchBend} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class StretchBendPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(StretchBendPotentialEnergy.class.getName());


  /**
   * Internal list of StretchBend instances.
   */
  private final List<StretchBend> stretchBends = new ArrayList<>();

  /**
   * Create a StretchBendPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public StretchBendPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a StretchBendPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public StretchBendPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a StretchBendPotentialEnergy initialized with a list of stretch-bends and force group.
   *
   * @param name         Name for this term.
   * @param forceGroup   Force group identifier.
   * @param stretchBends List of StretchBend instances to add (null-safe).
   */
  public StretchBendPotentialEnergy(String name, int forceGroup, List<StretchBend> stretchBends) {
    super(name, forceGroup);
    if (stretchBends != null) {
      Collections.sort(stretchBends);
      this.stretchBends.addAll(stretchBends);
      logger.info(format("  Stretch-Bends:                     %10d", getNumberOfStretchBends()));
    }
  }

  /**
   * Get the number of terms in this potential energy term.
   *
   * @return The number of terms, which is the same as the number of stretch-bends.
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfStretchBends();
  }

  /**
   * Get an array of BondedTerms in this term.
   *
   * @return Array of BondedTerms, which are actually StretchBends in this case.
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getStretchBendArray();
  }

  /**
   * Create a StretchBendPotentialEnergy initialized with a collection of stretch-bends.
   *
   * @param name         Name for this term (may be null).
   * @param stretchBends Collection of StretchBend instances to add (null-safe).
   */
  public StretchBendPotentialEnergy(String name, Collection<StretchBend> stretchBends) {
    super(name);
    if (stretchBends != null) {
      this.stretchBends.addAll(stretchBends);
    }
  }

  /**
   * Add a StretchBend to this term.
   *
   * @param stretchBend StretchBend to add (ignored if null).
   * @return true if the stretch-bend was added.
   */
  public boolean addStretchBend(StretchBend stretchBend) {
    if (stretchBend == null) {
      return false;
    }
    return stretchBends.add(stretchBend);
  }

  /**
   * Add an array of StretchBends to this term.
   *
   * @param stretchBends Array of StretchBend instances to add.
   * @return true if the stretch-bends were added.
   */
  public boolean addStretchBends(StretchBend[] stretchBends) {
    if (stretchBends == null) {
      return false;
    }
    Collections.addAll(this.stretchBends, stretchBends);
    return true;
  }

  /**
   * Add a list of StretchBends to this term.
   *
   * @param stretchBends List of StretchBend instances to add.
   * @return true if the stretch-bends were added.
   */
  public boolean addStretchBends(List<StretchBend> stretchBends) {
    if (stretchBends == null) {
      return false;
    }
    this.stretchBends.addAll(stretchBends);
    return true;
  }

  /**
   * Remove a StretchBend from this term.
   *
   * @param stretchBend StretchBend to remove (ignored if null).
   * @return true if the stretch-bend was present and removed.
   */
  public boolean removeStretchBend(StretchBend stretchBend) {
    if (stretchBend == null) {
      return false;
    }
    return stretchBends.remove(stretchBend);
  }

  /**
   * Get the StretchBend at a given index.
   *
   * @param index Index in the internal list.
   * @return StretchBend at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public StretchBend getStretchBend(int index) {
    return stretchBends.get(index);
  }

  /**
   * Get an unmodifiable view of the StretchBends in this term.
   *
   * @return Unmodifiable List of StretchBends.
   */
  public List<StretchBend> getStretchBends() {
    return Collections.unmodifiableList(stretchBends);
  }

  /**
   * Get an array of StretchBends in this term.
   *
   * @return Array of StretchBends.
   */
  public StretchBend[] getStretchBendArray() {
    return stretchBends.toArray(new StretchBend[0]);
  }

  /**
   * Get the number of StretchBends in this term.
   *
   * @return The number of StretchBends.
   */
  public int getNumberOfStretchBends() {
    return stretchBends.size();
  }

  /**
   * Get a formatted string representing the energy expression for Stretch-Bend interactions.
   * @return String representing the Stretch-Bend energy expression.
   */
  public static String getStretchBendEnergyString() {
    return format("(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))", 180.0 / PI);
  }

  /**
   * Log the details of Stretch-Bend interactions.
   */
  @Override
  public void log() {
    if (getNumberOfStretchBends() <= 0) {
      return;
    }
    logger.info("\n Stretch-Bend Interactions:");
    for (StretchBend stretchBend : getStretchBends()) {
      logger.info(" Stretch-Bend \t" + stretchBend.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfStretchBends() <= 0) {
      return "";
    }
    return format("REMARK   3   %s %g (%d)\n", "STRETCH-BEND               : ", getEnergy(), getNumberOfStretchBends());
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f\n", "Stretch-Bend      ",
        getEnergy(), getNumberOfStretchBends(), getTime());
  }

}
