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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete3DFunction_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete3DFunction_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete3DFunction_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Discrete3DFunction_setFunctionParameters;

/**
 * This is a TabulatedFunction that computes a discrete three dimensional function f(x,y,z).
 * To evaluate it, x, y, and z are each rounded to the nearest integer and the table element with those
 * indices is returned. If any index is outside the range [0, size), the result is undefined.
 */
public class Discrete3DFunction extends TabulatedFunction {

  /**
   * Create a Discrete3DFunction f(x,y,z) based on a set of tabulated values.
   *
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param zsize  the number of table elements along the z direction
   * @param values the tabulated values of the function f(x,y,z), ordered so that
   *               values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
   */
  public Discrete3DFunction(int xsize, int ysize, int zsize, PointerByReference values) {
    super(OpenMM_Discrete3DFunction_create(xsize, ysize, zsize, values));
  }

  /**
   * Destroy the discrete 3D function.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_Discrete3DFunction_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param zsize  the number of table elements along the z direction
   * @param values the tabulated values of the function f(x,y,z), ordered so that
   *               values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
   */
  public void getFunctionParameters(IntByReference xsize, IntByReference ysize, IntByReference zsize, PointerByReference values) {
    OpenMM_Discrete3DFunction_getFunctionParameters(pointer, xsize, ysize, zsize, values);
  }

  /**
   * Get the parameters for the tabulated function.
   *
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param zsize  the number of table elements along the z direction
   * @param values the tabulated values of the function f(x,y,z), ordered so that
   *               values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
   */
  public void getFunctionParameters(IntBuffer xsize, IntBuffer ysize, IntBuffer zsize, PointerByReference values) {
    OpenMM_Discrete3DFunction_getFunctionParameters(pointer, xsize, ysize, zsize, values);
  }

  /**
   * Set the parameters for the tabulated function.
   *
   * @param xsize  the number of table elements along the x direction
   * @param ysize  the number of table elements along the y direction
   * @param zsize  the number of table elements along the z direction
   * @param values the tabulated values of the function f(x,y,z), ordered so that
   *               values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
   */
  public void setFunctionParameters(int xsize, int ysize, int zsize, PointerByReference values) {
    OpenMM_Discrete3DFunction_setFunctionParameters(pointer, xsize, ysize, zsize, values);
  }

}