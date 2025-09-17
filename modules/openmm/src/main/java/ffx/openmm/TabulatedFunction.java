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
package ffx.openmm;

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_TabulatedFunction_getPeriodic;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_TabulatedFunction_getUpdateCount;

/**
 * A TabulatedFunction uses a set of tabulated values to define a mathematical function.
 * It can be used by various custom forces.
 * <p>
 * TabulatedFunction is an abstract class with concrete subclasses for more specific
 * types of functions. There are subclasses for:
 *
 * <ul>
 * <li>1, 2, and 3 dimensional functions. The dimensionality of a function means
 * the number of input arguments it takes.</li>
 * <li>Continuous and discrete functions. A continuous function is interpolated by
 * fitting a natural cubic spline to the tabulated values. A discrete function is
 * only defined for integer values of its arguments (that is, at the tabulated points),
 * and does not try to interpolate between them. Discrete function can be evaluated
 * more quickly than continuous ones.</li>
 * </ul>
 */
public abstract class TabulatedFunction {

  /**
   * The pointer is allocated and deallocated by classes that extend TabulatedFunction.
   */
  protected PointerByReference pointer;

  /**
   * Constructor for TabulatedFunction.
   *
   * @param pointer Pointer to the OpenMM TabulatedFunction.
   * @throws IllegalArgumentException if the pointer is null.
   */
  public TabulatedFunction(PointerByReference pointer) {
    if (pointer == null || pointer.getValue() == null) {
      throw new IllegalArgumentException("Pointer cannot be null.");
    }
    this.pointer = pointer;
  }

  /**
   * Destroy the tabulated function.
   */
  public abstract void destroy();

  /**
   * Get the pointer to the OpenMM TabulatedFunction.
   *
   * @return The pointer to the OpenMM TabulatedFunction.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get the periodicity status of the tabulated function.
   */
  public boolean getPeriodic() {
    int periodic = OpenMM_TabulatedFunction_getPeriodic(pointer);
    return periodic == OpenMM_True;
  }

  /**
   * Get the value of a counter that is updated every time setFunctionParameters()
   * is called. This provides a fast way to detect when a function has changed.
   */
  public int getUpdateCount() {
    return OpenMM_TabulatedFunction_getUpdateCount(pointer);
  }
}