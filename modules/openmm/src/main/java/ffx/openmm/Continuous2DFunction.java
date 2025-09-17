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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous2DFunction_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous2DFunction_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous2DFunction_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous2DFunction_setFunctionParameters;

/**
 * This is a TabulatedFunction that computes a continuous two dimensional function.
 */
public class Continuous2DFunction extends TabulatedFunction {

  /**
   * Create a Continuous2DFunction f(x,y) based on a set of tabulated values.
   *
   * @param values   the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
   *                 and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
   *                 The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
   *                 values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
   * @param xsize    the number of table elements along the x direction
   * @param ysize    the number of table elements along the y direction
   * @param xmin     the value of x corresponding to the first element of values
   * @param xmax     the value of x corresponding to the last element of values
   * @param ymin     the value of y corresponding to the first element of values
   * @param ymax     the value of y corresponding to the last element of values
   * @param periodic whether the interpolated function is periodic
   */
  public Continuous2DFunction(PointerByReference values, int xsize, int ysize, double xmin, double xmax, double ymin, double ymax, boolean periodic) {
    super(OpenMM_Continuous2DFunction_create(xsize, ysize, values, xmin, xmax, ymin, ymax, periodic ? 1 : 0));
  }

  /**
   * Destroy the continuous 2D function.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_Continuous2DFunction_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param values the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
   *               and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
   *               The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
   *               values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param xmin   the value of x corresponding to the first element of values
   * @param xmax   the value of x corresponding to the last element of values
   * @param ymin   the value of y corresponding to the first element of values
   * @param ymax   the value of y corresponding to the last element of values
   */
  public void getFunctionParameters(PointerByReference values, IntByReference xsize, IntByReference ysize,
                                    DoubleByReference xmin, DoubleByReference xmax,
                                    DoubleByReference ymin, DoubleByReference ymax) {
    OpenMM_Continuous2DFunction_getFunctionParameters(pointer, xsize, ysize, values, xmin, xmax, ymin, ymax);
  }

  /**
   * Set the parameters for the tabulated function.
   *
   * @param values the tabulated values of the function f(x,y) at xsize uniformly spaced values of x between xmin
   *               and xmax, and ysize values of y between ymin and ymax.  A natural cubic spline is used to interpolate between the tabulated values.
   *               The function is assumed to be zero when x or y is outside its specified range.  The values should be ordered so that
   *               values[i+xsize*j] = f(x_i,y_j), where x_i is the i'th uniformly spaced value of x.  This must be of length xsize*ysize.
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param xmin   the value of x corresponding to the first element of values
   * @param xmax   the value of x corresponding to the last element of values
   * @param ymin   the value of y corresponding to the first element of values
   * @param ymax   the value of y corresponding to the last element of values
   */
  public void setFunctionParameters(PointerByReference values, int xsize, int ysize,
                                    double xmin, double xmax, double ymin, double ymax) {
    OpenMM_Continuous2DFunction_setFunctionParameters(pointer, xsize, ysize, values, xmin, xmax, ymin, ymax);
  }
}