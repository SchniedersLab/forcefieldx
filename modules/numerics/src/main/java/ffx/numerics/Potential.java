// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.numerics;

import java.util.Collections;
import java.util.List;

/**
 * The Potential interface defines methods required by an optimizer or molecular dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Potential extends OptimizationInterface {

  /**
   * getAcceleration.
   *
   * @param acceleration an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  double[] getAcceleration(double[] acceleration);

  /**
   * Returns the list of Constraints associated with this Potential. The default implementation
   * returns an empty list.
   *
   * @return All Constraints.
   */
  default List<Constraint> getConstraints() {
    return Collections.emptyList();
  }

  /**
   * Get the Potential Energy terms that is active.
   *
   * @return the STATE
   */
  STATE getEnergyTermState();

  /**
   * Set the Potential Energy terms that should be active.
   *
   * @param state include FAST varying energy terms, SLOW varying energy terms or BOTH.
   */
  void setEnergyTermState(STATE state);

  /**
   * Get the mass of each degree of freedom. This is required for molecular dynamics.
   *
   * @return The mass of each degree of freedom.
   */
  double[] getMass();

  /**
   * getPreviousAcceleration.
   *
   * @param previousAcceleration an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  double[] getPreviousAcceleration(double[] previousAcceleration);

  /**
   * Get the type of all variables.
   *
   * @return The VARIABLE_TYPE of each variable.
   */
  VARIABLE_TYPE[] getVariableTypes();

  /**
   * getVelocity.
   *
   * @param velocity an array of {@link double} objects.
   * @return an array of {@link double} objects.
   */
  double[] getVelocity(double[] velocity);

  /**
   * setAcceleration.
   *
   * @param acceleration an array of {@link double} objects.
   */
  void setAcceleration(double[] acceleration);

  /**
   * setPreviousAcceleration.
   *
   * @param previousAcceleration an array of {@link double} objects.
   */
  void setPreviousAcceleration(double[] previousAcceleration);

  /**
   * setVelocity.
   *
   * @param velocity an array of {@link double} objects.
   */
  void setVelocity(double[] velocity);

  /**
   * Writes additional restart information, if any (e.g. OST histogram and lambda restart files). The
   * recursive flag should generally only be true for the top-level Potential called.
   *
   * @param recursive Whether to have all underlying Potentials write additional restart info.
   */
  default void writeAdditionalRestartInfo(boolean recursive) {
    if (recursive) {
      getUnderlyingPotentials().forEach((Potential p) -> p.writeAdditionalRestartInfo(false));
    } // Else, no-op.
  }
  /**/

  /**
   * Recognized variables currently include Cartesian coordinates and OTHER.
   */
  enum VARIABLE_TYPE {
    X,
    Y,
    Z,
    OTHER
  }

  /**
   * Set the state of the Potential to include FAST varying energy terms, SLOW varying energy terms
   * or BOTH.
   */
  enum STATE {
    FAST,
    SLOW,
    BOTH
  }
}
