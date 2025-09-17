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
import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous1DFunction_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous1DFunction_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous1DFunction_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Continuous1DFunction_setFunctionParameters;

/**
 * This is a TabulatedFunction that computes a continuous one dimensional function.
 */
public class Continuous1DFunction extends TabulatedFunction {

  /**
   * Create a Continuous1DFunction f(x) based on a set of tabulated values.
   *
   * @param values   the tabulated values of the function f(x) at uniformly spaced values of x between min
   *                 and max.  A natural cubic spline is used to interpolate between the tabulated values.
   *                 The function is assumed to be zero for x &lt; min or x &gt; max.
   * @param min      the value of x corresponding to the first element of values
   * @param max      the value of x corresponding to the last element of values
   * @param periodic whether the interpolated function is periodic
   */
  public Continuous1DFunction(PointerByReference values, double min, double max, boolean periodic) {
    super(OpenMM_Continuous1DFunction_create(values, min, max, periodic ? 1 : 0));
  }

  /**
   * Destroy the continuous 1D function.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_Continuous1DFunction_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param values the tabulated values of the function f(x) at uniformly spaced values of x between min
   *               and max.  A natural cubic spline is used to interpolate between the tabulated values.
   *               The function is assumed to be zero for x &lt; min or x &gt; max.
   * @param min    the value of x corresponding to the first element of values
   * @param max    the value of x corresponding to the last element of values
   */
  public void getFunctionParameters(PointerByReference values, DoubleByReference min, DoubleByReference max) {
    OpenMM_Continuous1DFunction_getFunctionParameters(pointer, values, min, max);
  }

  /**
   * Set the parameters for the tabulated function.
   *
   * @param values the tabulated values of the function f(x) at uniformly spaced values of x between min
   *               and max.  A natural cubic spline is used to interpolate between the tabulated values.
   *               The function is assumed to be zero for x &lt; min or x &gt; max.
   * @param min    the value of x corresponding to the first element of values
   * @param max    the value of x corresponding to the last element of values
   */
  public void setFunctionParameters(PointerByReference values, double min, double max) {
    OpenMM_Continuous1DFunction_setFunctionParameters(pointer, values, min, max);
  }
}