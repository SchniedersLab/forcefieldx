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
package ffx.potential.terms;

import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.parameters.BondType;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

/**
 * Bond potential energy term using {@link ffx.potential.bonded.Bond} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class BondPotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(BondPotentialEnergy.class.getName());


  /**
   * Internal list of Bond instances.
   */
  private final List<Bond> bonds = new ArrayList<>();

  /**
   * Create a BondPotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public BondPotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a BondPotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public BondPotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a BondPotentialEnergy initialized with a list of bonds and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Force group identifier.
   * @param bonds      List of Bond instances to add (null-safe).
   */
  public BondPotentialEnergy(String name, int forceGroup, List<Bond> bonds) {
    super(name, forceGroup);
    if (bonds != null) {
      Collections.sort(bonds);
      this.bonds.addAll(bonds);
      logger.info(format("  Bonds:                             %10d", getNumberOfBonds()));
    }
  }

  /**
   * Get the number of terms in this potential energy term.
   *
   * @return The number of terms, which is the same as the number of bonds.
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfBonds();
  }

  /**
   * Get an array of BondedTerms in this term.
   *
   * @return Array of BondedTerms, which are actually Bonds in this case.
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getBondArray();
  }

  /**
   * Create a BondPotentialEnergy initialized with a collection of bonds.
   *
   * @param name  Name for this term (may be null).
   * @param bonds Collection of Bond instances to add (null-safe).
   */
  public BondPotentialEnergy(String name, Collection<Bond> bonds) {
    super(name);
    if (bonds != null) {
      this.bonds.addAll(bonds);
    }
  }

  /**
   * Add a Bond to this term.
   *
   * @param bond Bond to add (ignored if null).
   * @return true if the bond was added.
   */
  public boolean addBond(Bond bond) {
    if (bond == null) {
      return false;
    }
    return bonds.add(bond);
  }

  /**
   * Add an array of Bonds to this term.
   *
   * @param bonds Array of Bond instances to add.
   * @return true if the bonds were added.
   */
  public boolean addBonds(Bond[] bonds) {
    if (bonds == null) {
      return false;
    }
    Collections.addAll(this.bonds, bonds);
    return true;
  }

  /**
   * Add a list of Bonds to this term.
   *
   * @param bonds List of Bond instances to add.
   * @return true if the bonds were added.
   */
  public boolean addBonds(List<Bond> bonds) {
    if (bonds == null) {
      return false;
    }
    this.bonds.addAll(bonds);
    return true;
  }

  /**
   * Remove a Bond from this term.
   *
   * @param bond Bond to remove (ignored if null).
   * @return true if the bond was present and removed.
   */
  public boolean removeBond(Bond bond) {
    if (bond == null) {
      return false;
    }
    return bonds.remove(bond);
  }

  /**
   * Get the Bond at a given index.
   *
   * @param index Index in the internal list.
   * @return Bond at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public Bond getBond(int index) {
    return bonds.get(index);
  }

  /**
   * Get an unmodifiable view of the Bonds in this term.
   *
   * @return Unmodifiable List of Bonds.
   */
  public List<Bond> getBonds() {
    return Collections.unmodifiableList(bonds);
  }

  /**
   * Get an array of Bonds in this term.
   *
   * @return Array of Bonds.
   */
  public Bond[] getBondArray() {
    return bonds.toArray(new Bond[0]);
  }

  /**
   * Get the number of Bonds in this term.
   *
   * @return The number of Bonds.
   */
  public int getNumberOfBonds() {
    return bonds.size();
  }

  public String getBondEnergyString() {
    BondType bondType = bonds.getFirst().getBondType();
    String energy;
    if (bondType.bondFunction == BondType.BondFunction.QUARTIC) {
      energy = format("""
              k*(d^2 + %.15g*d^3 + %.15g*d^4);
              d=r-r0;
              """,
          bondType.cubic / OpenMM_NmPerAngstrom,
          bondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
    } else {
      energy = """
          k*(d^2);
          d=r-r0;
          """;
    }
    return energy;
  }

  /**
   * Log the details of Bond interactions.
   */
  @Override
  public void log() {
    if (getNumberOfBonds() <= 0) {
      return;
    }
    logger.info("\n Bond Stretching Interactions:");
    for (Bond bond : getBonds()) {
      logger.info(" Bond \t" + bond.toString());
    }
  }

  @Override
  public String toPDBString() {
    if (getNumberOfBonds() <= 0) {
      return "";
    }
    StringBuilder sb = new StringBuilder();
    sb.append(format("REMARK   3   %s %g (%d)\n", "BOND STRETCHING            : ", getEnergy(), getNumberOfBonds()));
    sb.append(format("REMARK   3   %s %g\n", "BOND RMSD                  : ", getRMSD()));
    return sb.toString();
  }

  @Override
  public String toString() {
    return format("  %s %20.8f %12d %12.3f (%8.5f)\n", "Bond Stretching   ",
        getEnergy(), getNumberOfBonds(), getTime(), getRMSD());
  }

}
