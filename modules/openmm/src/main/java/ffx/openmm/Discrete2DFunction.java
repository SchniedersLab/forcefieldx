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

import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete2DFunction_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete2DFunction_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete2DFunction_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete2DFunction_setFunctionParameters;

/**
 * A Discrete2DFunction uses discrete values to define a function based on a discrete set of
 * tabulated values in two dimensions. This is useful for defining functions from experimental
 * data or other tabulated sources that depend on two independent variables where no
 * interpolation is desired.
 * <p>
 * The function is defined by a grid of (x, y, f(x,y)) values, and values are returned exactly
 * as tabulated without interpolation. The function can optionally be periodic in either or
 * both dimensions, meaning that values outside the tabulated range are computed by wrapping
 * around to the other end of the table.
 */
public class Discrete2DFunction extends TabulatedFunction {

  /**
   * Create a Discrete2DFunction.
   *
   * @param xsize  The number of table elements along the x direction.
   * @param ysize  The number of table elements along the y direction.
   * @param values The tabulated values of the function f(x,y).
   */
  public Discrete2DFunction(int xsize, int ysize, PointerByReference values) {
    super(OpenMM_Discrete2DFunction_create(xsize, ysize, values));
  }

  /**
   * Destroy the discrete 2D function.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_Discrete2DFunction_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param xsize  The number of table elements along the x direction (output).
   * @param ysize  The number of table elements along the y direction (output).
   * @param values The tabulated values of the function (output).
   */
  public void getFunctionParameters(IntByReference xsize, IntByReference ysize, PointerByReference values) {
    OpenMM_Discrete2DFunction_getFunctionParameters(pointer, xsize, ysize, values);
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param xsize  The number of table elements along the x direction (output).
   * @param ysize  The number of table elements along the y direction (output).
   * @param values The tabulated values of the function (output).
   */
  public void getFunctionParameters(IntBuffer xsize, IntBuffer ysize, PointerByReference values) {
    OpenMM_Discrete2DFunction_getFunctionParameters(pointer, xsize, ysize, values);
  }

  /**
   * Set the parameters for the tabulated function.
   *
   * @param xsize  The number of table elements along the x direction.
   * @param ysize  The number of table elements along the y direction.
   * @param values The tabulated values of the function.
   */
  public void setFunctionParameters(int xsize, int ysize, PointerByReference values) {
    OpenMM_Discrete2DFunction_setFunctionParameters(pointer, xsize, ysize, values);
  }

}