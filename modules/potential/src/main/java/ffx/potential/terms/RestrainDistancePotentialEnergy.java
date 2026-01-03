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

import ffx.crystal.Crystal;
import ffx.numerics.switching.ConstantSwitch;
import ffx.numerics.switching.UnivariateFunctionFactory;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.RestrainDistance;
import ffx.potential.parameters.BondType;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Restrain-Distance potential energy term using {@link ffx.potential.bonded.RestrainDistance} instances.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainDistancePotentialEnergy extends EnergyTerm {

  private static final Logger logger = Logger.getLogger(RestrainDistancePotentialEnergy.class.getName());

  /**
   * Log the details of Restrain Distance interactions.
   */
  @Override
  public void log() {
    if (getNumberOfRestrainDistances() <= 0) {
      return;
    }
    logger.info("\n Restrain Distance Interactions:");
    for (RestrainDistance restrainDistance : getRestrainDistances()) {
      logger.info(" Restrain Distance \t" + restrainDistance.toString());
    }
  }

  /**
   * Internal list of RestrainDistance instances.
   */
  private final List<RestrainDistance> restrainDistances = new ArrayList<>();

  /**
   * Create a RestrainDistancePotentialEnergy with the provided name.
   *
   * @param name Name for this term.
   */
  public RestrainDistancePotentialEnergy(String name) {
    super(name);
  }

  /**
   * Create a RestrainDistancePotentialEnergy with the provided name and force group.
   *
   * @param name       Name for this term.
   * @param forceGroup Integer force group identifier.
   */
  public RestrainDistancePotentialEnergy(String name, int forceGroup) {
    super(name, forceGroup);
  }

  /**
   * Create a RestrainDistancePotentialEnergy initialized with a list of terms and force group.
   *
   * @param name              Name for this term.
   * @param forceGroup        Force group identifier.
   * @param restrainDistances List of RestrainDistance instances (null-safe).
   */
  public RestrainDistancePotentialEnergy(String name, int forceGroup, List<RestrainDistance> restrainDistances) {
    super(name, forceGroup);
    if (restrainDistances != null) {
      Collections.sort(restrainDistances);
      this.restrainDistances.addAll(restrainDistances);
      logger.info(String.format("  Restrain Distances:                 %10d", getNumberOfRestrainDistances()));
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfTerms() {
    return getNumberOfRestrainDistances();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BondedTerm[] getBondedTermsArray() {
    return getRestrainDistanceArray();
  }

  /**
   * Create a RestrainDistancePotentialEnergy initialized with a collection of terms.
   *
   * @param name              Name for this term (may be null).
   * @param restrainDistances Collection of RestrainDistance instances to add (null-safe).
   */
  public RestrainDistancePotentialEnergy(String name, Collection<RestrainDistance> restrainDistances) {
    super(name);
    if (restrainDistances != null) {
      this.restrainDistances.addAll(restrainDistances);
    }
  }

  /**
   * Add a RestrainDistance to this term.
   *
   * @param restrainDistance RestrainDistance to add (ignored if null).
   * @return true if it was added.
   */
  public boolean addRestrainDistance(RestrainDistance restrainDistance) {
    if (restrainDistance == null) {
      return false;
    }
    return restrainDistances.add(restrainDistance);
  }

  /**
   * Add an array of RestrainDistances to this term.
   *
   * @param restrainDistances Array of RestrainDistance instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainDistances(RestrainDistance[] restrainDistances) {
    if (restrainDistances == null) {
      return false;
    }
    Collections.addAll(this.restrainDistances, restrainDistances);
    return true;
  }

  /**
   * Add a list of RestrainDistances to this term.
   *
   * @param restrainDistances List of RestrainDistance instances to add.
   * @return true if they were added.
   */
  public boolean addRestrainDistances(List<RestrainDistance> restrainDistances) {
    if (restrainDistances == null) {
      return false;
    }
    this.restrainDistances.addAll(restrainDistances);
    return true;
  }

  /**
   * Remove a RestrainDistance from this term.
   *
   * @param restrainDistance RestrainDistance to remove (ignored if null).
   * @return true if it was present and removed.
   */
  public boolean removeRestrainDistance(RestrainDistance restrainDistance) {
    if (restrainDistance == null) {
      return false;
    }
    return restrainDistances.remove(restrainDistance);
  }

  /**
   * Get the RestrainDistance at a given index.
   *
   * @param index Index in the internal list.
   * @return RestrainDistance at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public RestrainDistance getRestrainDistance(int index) {
    return restrainDistances.get(index);
  }

  /**
   * Get an unmodifiable view of the RestrainDistances in this term.
   *
   * @return Unmodifiable List of RestrainDistances.
   */
  public List<RestrainDistance> getRestrainDistances() {
    return Collections.unmodifiableList(restrainDistances);
  }

  /**
   * Get an array of RestrainDistances in this term.
   *
   * @return Array of RestrainDistances.
   */
  public RestrainDistance[] getRestrainDistanceArray() {
    return restrainDistances.toArray(new RestrainDistance[0]);
  }

  /**
   * Get the number of RestrainDistances in this term.
   *
   * @return The number of RestrainDistances.
   */
  public int getNumberOfRestrainDistances() {
    return restrainDistances.size();
  }

  /**
   * Returns a list of restraint distances filtered by the specified bond function. If the bondFunction
   * is null, it returns all restrained bonds.
   *
   * @param bondFunction the type of bond function.
   * @return a {@link java.util.List} object.
   */
  public List<RestrainDistance> getRestrainDistances(@Nullable BondType.BondFunction bondFunction) {
    // If the bondFunction is null, return all restrained bonds.
    if (bondFunction == null) {
      return restrainDistances;
    }
    // Otherwise, return only the restraint bonds with the specified bond function.
    List<RestrainDistance> list = new ArrayList<>();
    for (RestrainDistance restrainDistance : restrainDistances) {
      if (restrainDistance.getBondType().bondFunction == bondFunction) {
        list.add(restrainDistance);
      }
    }
    if (!list.isEmpty()) {
      return list;
    }
    return null;
  }

  @Override
  public String toPDBString() {
    if (getNumberOfRestrainDistances() <= 0) {
      return "";
    }
    return String.format("REMARK   3   %s %g (%d)\n", "RESTRAIN DISTANCE          : ", getEnergy(), getNumberOfRestrainDistances());
  }

  @Override
  public String toString() {
    return String.format("  %s %20.8f %12d %12.3f\n", "Restrain Distance ",
        getEnergy(), getNumberOfRestrainDistances(), getTime());
  }

  /**
   * Method to parse the restrain-distance records.
   *
   * @param properties Configuration properties.
   */
  public static RestrainDistance[] configureRestrainDistances(
      CompositeConfiguration properties,
      Atom[] atoms, Crystal crystal, boolean lambdaTerm) {
    List<RestrainDistance> restrainDistanceList = new ArrayList<>();
    String[] bondRestraints = properties.getStringArray("restrain-distance");
    for (String bondRest : bondRestraints) {
      try {
        String[] toks = bondRest.split("\\s+");
        if (toks.length < 2) {
          throw new IllegalArgumentException(
              format(" restrain-distance value %s could not be parsed!", bondRest));
        }
        // Internally, everything starts with 0, but restrain distance starts at 1, so that 1 has to
        // be subtracted
        int at1 = Integer.parseInt(toks[0]) - 1;
        int at2 = Integer.parseInt(toks[1]) - 1;

        double forceConst = 100.0;
        double flatBottomRadius = 0;
        Atom a1 = atoms[at1];
        Atom a2 = atoms[at2];

        if (toks.length > 2) {
          forceConst = Double.parseDouble(toks[2]);
        }
        double dist;
        switch (toks.length) {
          case 2:
          case 3:
            double[] xyz1 = new double[3];
            xyz1 = a1.getXYZ(xyz1);
            double[] xyz2 = new double[3];
            xyz2 = a2.getXYZ(xyz2);
            // Current distance between restrained atoms
            dist = crystal.minDistOverSymOps(xyz1, xyz2);
            break;
          case 4:
            dist = Double.parseDouble(toks[3]);
            break;
          case 5:
          default:
            double minDist = Double.parseDouble(toks[3]);
            double maxDist = Double.parseDouble(toks[4]);
            dist = 0.5 * (minDist + maxDist);
            flatBottomRadius = 0.5 * Math.abs(maxDist - minDist);
            break;
        }

        UnivariateSwitchingFunction switchF;
        double lamStart = RestrainDistance.DEFAULT_RB_LAM_START;
        double lamEnd = RestrainDistance.DEFAULT_RB_LAM_END;
        if (toks.length > 5) {
          int offset = 5;
          if (toks[5].matches("^[01](?:\\.[0-9]*)?")) {
            offset = 6;
            lamStart = Double.parseDouble(toks[5]);
            if (toks[6].matches("^[01](?:\\.[0-9]*)?")) {
              offset = 7;
              lamEnd = Double.parseDouble(toks[6]);
            }
          }
          switchF = UnivariateFunctionFactory.parseUSF(toks, offset);
        } else {
          switchF = new ConstantSwitch();
        }

        RestrainDistance restrainDistance = createRestrainDistance(a1, a2, dist,
            forceConst, flatBottomRadius, lamStart, lamEnd, switchF, lambdaTerm, crystal);
        restrainDistanceList.add(restrainDistance);
      } catch (Exception ex) {
        logger.info(format(" Exception in parsing restrain-distance: %s", ex));
      }
    }
    if (!restrainDistanceList.isEmpty()) {
      return restrainDistanceList.toArray(new RestrainDistance[0]);
    } else {
      return null;
    }
  }

  /**
   * setRestrainDistance
   *
   * @param a1                a {@link ffx.potential.bonded.Atom} object.
   * @param a2                a {@link ffx.potential.bonded.Atom} object.
   * @param distance          a double.
   * @param forceConstant     the force constant in kcal/mole.
   * @param flatBottom        Radius of a flat-bottom potential in Angstroms.
   * @param lamStart          At what lambda does the restraint begin to take effect?
   * @param lamEnd            At what lambda does the restraint hit full strength?
   * @param switchingFunction Switching function to use as a lambda dependence.
   * @param lambdaTerm        Restrain distance is a function of lambda.
   * @param crystal           The unit cel and space group information.
   */
  private static RestrainDistance createRestrainDistance(Atom a1, Atom a2, double distance, double forceConstant,
                                                         double flatBottom, double lamStart, double lamEnd,
                                                         UnivariateSwitchingFunction switchingFunction,
                                                         boolean lambdaTerm, Crystal crystal) {
    boolean rbLambda = !(switchingFunction instanceof ConstantSwitch) && lambdaTerm;
    RestrainDistance restrainDistance = new RestrainDistance(a1, a2, crystal, rbLambda, lamStart, lamEnd, switchingFunction);
    int[] classes = {a1.getAtomType().atomClass, a2.getAtomType().atomClass};
    if (flatBottom != 0) {
      BondType bondType = new BondType(classes, forceConstant, distance,
          BondType.BondFunction.FLAT_BOTTOM_HARMONIC, flatBottom);
      restrainDistance.setBondType(bondType);
    } else {
      BondType bondType = new BondType(classes, forceConstant, distance,
          BondType.BondFunction.HARMONIC);
      restrainDistance.setBondType(bondType);
    }
    return restrainDistance;
  }
}
