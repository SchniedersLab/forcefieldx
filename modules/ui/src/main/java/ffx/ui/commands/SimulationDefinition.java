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
package ffx.ui.commands;

import java.io.Serializable;

/**
 * The SimulationDefinition class is a serializable wrapper that specifies an FFX simulation.
 *
 * @author Michael J. Schnieders
 */
public class SimulationDefinition implements Serializable {

  private static final long serialVersionUID = 1L;
  public boolean read = true;
  // System definition
  public int numatoms;
  public int numkeys;
  public String file;
  public String forcefield;
  public String[] keywords;
  public double[][] coordinates;
  public int[][] connectivity;
  public int[] types;
  public String[] name;
  public String[] story;
  public double[] charge;
  public double[] mass;
  public int[] atomic;

  /**
   * Constructor that allocates space for a simulation definition.
   *
   * @param a The number of atoms
   * @param k The number of keywords
   */
  public SimulationDefinition(int a, int k) {
    numatoms = a;
    numkeys = k;
    keywords = new String[k];
    coordinates = new double[3][a];
    connectivity = new int[4][a];
    types = new int[a];
    name = new String[a];
    story = new String[a];
    charge = new double[a];
    mass = new double[a];
    atomic = new int[a];
  }

  /** print */
  public void print() {
    System.out.println(this.toString());
  }

  /**
   * toString
   *
   * @return a {@link java.lang.String} object.
   */
  public String toString() {
    return new String("Atoms: " + numatoms + " Keywords: " + numkeys);
  }
}
