// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
 * The SimulationMessage class is used to pass simple messages between an FFXServer and its
 * FFXClient(s).
 *
 * @author Michael J. Schnieders
 */
public class SimulationMessage implements Serializable {

  private static final long serialVersionUID = 1L;
  /** Constant <code>SYSTEM=0</code> */
  public static int SYSTEM = 0;
  /** Constant <code>UPDATE=1</code> */
  public static int UPDATE = 1;
  /** Constant <code>CLOSING=2</code> */
  public static int CLOSING = 2;
  /** Constant <code>OK=3</code> */
  public static int OK = 3;

  private int message = 0;
  private int step = -1;
  private double time = -1.0;
  private int type = 0;

  /**
   * Constructor for SimulationMessage.
   *
   * @param m a int.
   */
  public SimulationMessage(int m) {
    message = m;
  }

  /**
   * Getter for the field <code>message</code>.
   *
   * @return a int.
   */
  public int getMessage() {
    return message;
  }

  /**
   * Setter for the field <code>message</code>.
   *
   * @param m a int.
   */
  public void setMessage(int m) {
    message = m;
  }

  /**
   * Getter for the field <code>step</code>.
   *
   * @return a int.
   */
  public int getStep() {
    return step;
  }

  /**
   * Setter for the field <code>step</code>.
   *
   * @param s a int.
   */
  public void setStep(int s) {
    step = s;
  }

  /**
   * Getter for the field <code>time</code>.
   *
   * @return a double.
   */
  public double getTime() {
    return time;
  }

  /**
   * Setter for the field <code>time</code>.
   *
   * @param t a double.
   */
  public void setTime(double t) {
    time = t;
  }

  /**
   * Getter for the field <code>type</code>.
   *
   * @return a int.
   */
  public int getType() {
    return type;
  }

  /**
   * Setter for the field <code>type</code>.
   *
   * @param t a int.
   */
  public void setType(int t) {
    type = t;
  }

  /** print */
  public void print() {
    System.out.println(toString());
  }

  /**
   * toString
   *
   * @return a {@link java.lang.String} object.
   */
  public String toString() {
    if (message == 0) {
      return new String("SYSTEM").intern();
    } else if (message == 1) {
      return new String("UPDATE").intern();
    } else if (message == 2) {
      return new String("CLOSING").intern();
    } else {
      return new String("OK").intern();
    }
  }
}
