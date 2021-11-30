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
package ffx.algorithms.groovy;

import static org.junit.Assert.assertEquals;

import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.misc.AlgorithmsTest;
import java.util.Arrays;
import java.util.Collection;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * Test NVT dynamics on a waterbox.
 *
 * @author Hernan V. Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsNVTTest extends AlgorithmsTest {

  private String info;
  private String filename;
  private double finalTemp;
  private double tempTolerance = 0.01;
  private double endTotalEnergy;
  private double energyTolerance = 0.01;

  public DynamicsNVTTest(String info, String filename, double finalTemp, double endTotalEnergy) {
    this.info = info;
    this.filename = filename;
    this.finalTemp = finalTemp;
    this.endTotalEnergy = endTotalEnergy;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "Water Box NVT", // info
                "ffx/algorithms/structures/waterbox_eq.xyz", // filename
                296.18798, // Final temperature.
                -24952.0595 // Final total energy
            }
        });
  }

  @Test
  public void testDynamicsNVT() {

    // Set-up the input arguments for the script.
    String[] args = {
        "-n", "10",
        "-t", "298.15",
        "-i", "VelocityVerlet",
        "-b", "Bussi",
        "-r", "0.001",
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    Dynamics dynamics = new Dynamics(binding).run();
    algorithmsScript = dynamics;
    MolecularDynamics molDyn = dynamics.getMolecularDynamics();

    // Assert that temperature is within tolerance at the end of the dynamics trajectory.
    assertEquals(info + " End temperature for NVT test", finalTemp, molDyn.getTemperature(),
        tempTolerance);

    // Assert that the end total energy is withing the tolerance at the end of the dynamics
    // trajectory.
    assertEquals(info + " End total energy for NVT test and set random seed",
        endTotalEnergy, molDyn.getTotalEnergy(), energyTolerance);
  }
}
