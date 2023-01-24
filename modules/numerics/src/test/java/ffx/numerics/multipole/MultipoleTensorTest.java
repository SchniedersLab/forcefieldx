// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.numerics.multipole;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collection;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Parameterized Test of the MultipoleTensor class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class MultipoleTensorTest {

  private final int order;
  private final int tensorCount;
  private final String info;

  // Water Dimer O-O interaction
  protected static final double[] xyz = {2.97338529, 0.0, 0.03546452, 0.0};
  // Rotated multipole at site I.
  protected final static double[] Qi = {-0.51966,
      0.06979198988239577, 0.0, 0.0289581620819011,
      0.024871041044109393, -0.1170771231287098, 0.09220608208460039,
      0.0, -0.03374891685535346, 0.0};
  // Rotated multipole at site K.
  protected final static double[] Qk = {-0.51966,
      0.05872406108747119, 0.0, 0.047549780780788455,
      0.048623413357888695, -0.1170771231287098, 0.06845370977082109,
      0.0, -0.04662811558421081, 0.0};

  // Permanent multipole Coulomb energy and torque.
  protected final static double permanentEnergy = 0.08861401274679;
  protected final static double[] permTorqueI = {0.0, 0.00039200904382, 0.0};
  protected final static double[] permTorqueK = {0.0, -0.00400122461628, 0.0};
  // Permanent multipole Ewald energy and torque.
  protected final static double permanentEnergyEwald = 0.001275693067120;
  protected final static double[] permTorqueIEwald = {0.0, -0.000304652548761, 0.0};
  protected final static double[] permTorqueKEwald = {0.0, -0.000745949268609, 0.0};

  // Water Dimer O-O Coulomb mutual polarization Energy, grad and torque in vacuum
  protected final static double scaleMutual = 1.0;
  protected final static double[] Ui = {0.04886563833303603, 0.0, -0.0018979726219775425};
  protected final static double[] Uk = {-0.040839567654139396, 0.0, -5.982126263609587E-4};
  protected final static double polarizationEnergyCoulomb = -0.002576831234958;
  protected final static double[] polarGradICoulomb = {-0.003233364149133, 0.0, 0.000081374809911};
  protected final static double[] polarTorqueICoulomb = {0.0, 0.000015759704950, 0.0};
  protected final static double[] polarTorqueKCoulomb = {0.0, 0.000340863064066, 0.0};
  // Water Dimer O-O Thole correction for mutual Polarization energy, grad and torque in vacuum.
  protected final static double polarizationEnergyThole = 0.000000036672468;
  protected final static double[] polarGradIThole = {0.000001005045968, 0.0, 0.000000090221488};
  protected final static double[] polarTorqueIThole = {0.0, 0.000000072500602, 0.0};
  protected final static double[] polarTorqueKThole = {0.0, 0.000000151677615, 0.0};

  // Water Dimer O-O Coulomb mutual polarization Energy, grad and torque under Ewald.
  protected final static double[] UiEwald = {0.017212550663556914, 0.0, 0.0026849371538929636};
  protected final static double[] UkEwald = {0.0020051189004858657, 0.0, 0.005476556978973183};
  protected final static double polarizationEnergyEwald = -0.000078870483232;
  protected final static double[] polarGradIEwald = {-0.000345793382509, 0.0, 0.000025382496373};
  protected final static double[] polarTorqueIEwald = {0.0, -0.000003385814434, 0.0};
  protected final static double[] polarTorqueKEwald = {0.0, 0.000080790360336, 0.0};

  public MultipoleTensorTest(String info, int order) {
    this.info = info;
    this.order = order;
    this.tensorCount = MultipoleTensor.tensorCount(order);
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "Order 6", 6
            },
            {
                "Order 5", 5
            },
            {
                "Order 4", 4
            }
        });
  }

  /**
   * Test of tensorCount method, of class MultipoleTensor.
   *
   * @since 1.0
   */
  @Test
  public void tensorCountTest() {
    int result = MultipoleTensor.tensorCount(order);
    assertEquals(info, tensorCount, result);
  }

  /** Test of ti method, of class MultipoleTensor. */
  @Test
  public void tensorIndexTest() {
    int dx = 1;
    int dy = 0;
    int dz = 0;
    int expResult = 1;
    int result = MultipoleTensor.ti(dx, dy, dz, order);
    assertEquals(info, expResult, result);
  }

}
