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

/** @author Hernan V Bernabe */
@RunWith(Parameterized.class)
public class DynamicsInitVelocitiesTest extends AlgorithmsTest {

  private String info;
  private String filename;
  private double endKineticEnergy;
  private double kineticEnergyTolerance = 5.0;
  private double endPotentialEnergy;
  private double potentialEnergyTolerance = 5.0;
  private double endTotalEnergy;
  private double totalEnergyTolerance = 5.0;

  public DynamicsInitVelocitiesTest(
      String info,
      String filename,
      double endKineticEnergy,
      double endPotentialEnergy,
      double endTotalEnergy) {

    this.info = info;
    this.filename = filename;
    this.endKineticEnergy = endKineticEnergy;
    this.endPotentialEnergy = endPotentialEnergy;
    this.endTotalEnergy = endTotalEnergy;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "Water Tiny Initialize Velocities", // info
                "ffx/algorithms/structures/watertiny.xyz", // filename
                64.8511, // endKineticEnergy
                -218.0760, // endPotentialEnergy
                -153.2249 // endTotalEnergy
            }
        });
  }

  @Test
  public void testDynamicsInitVelocities() {

    // Set-up the input arguments for the script.
    String[] args = {
        "-n", "10",
        "-t", "298.15",
        "-i", "VelocityVerlet",
        "-b", "Adiabatic",
        "-r", "0.001",
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    Dynamics dynamics = new Dynamics(binding).run();
    algorithmsScript = dynamics;

    MolecularDynamics molDyn = dynamics.getMolecularDynamics();

    // Assert that final energies are within the tolerance for the dynamics trajectory.
    assertEquals(
        info + "End kinetic energy for initializing velocities test",
        endKineticEnergy,
        molDyn.getKineticEnergy(),
        kineticEnergyTolerance);

    assertEquals(
        info + "End potential energy for initializing velocities test",
        endPotentialEnergy,
        molDyn.getPotentialEnergy(),
        potentialEnergyTolerance);

    assertEquals(
        info + "End total energy for initializing velocities test",
        endTotalEnergy,
        molDyn.getTotalEnergy(),
        totalEnergyTolerance);
  }
}
