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

import ffx.algorithms.dynamics.MolecularDynamicsOpenMM;
import ffx.algorithms.misc.AlgorithmsTest;
import java.util.Arrays;
import java.util.Collection;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/** @author Hernan V Bernabe */
@RunWith(Parameterized.class)
public class DynamicsOpenMMStochasticTest extends AlgorithmsTest {

  private final String info;
  private final String filename;
  private final double endKineticEnergy;
  private final double endPotentialEnergy;
  private final double endTotalEnergy;
  private final double tolerance = 5.0;

  public DynamicsOpenMMStochasticTest(
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
                "System OpenMM Stochastic (Starting Potential Energy = -35661.8041)",
                "ffx/algorithms/structures/waterbox_eq.xyz",
                11009.729434,
                -36039.826268,
                -25030.095644,
            }
        });
  }

  @Before
  public void before() {
    System.setProperty("platform", "omm");
  }

  @Test
  public void testDynamicsOpenMMStochastic() {
    if (!ffxOpenMM) {
      return;
    }

    // Set-up the input arguments for the script.
    String[] args = {
        "-n", "10",
        "-z", "1",
        "-t", "298.15",
        "-i", "Stochastic",
        "-r", "0.001",
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    Dynamics dynamics = new Dynamics(binding).run();
    algorithmsScript = dynamics;

    MolecularDynamicsOpenMM molDynOpenMM =
        (MolecularDynamicsOpenMM) dynamics.getMolecularDynamics();

    // Assert that the end energies are within the threshold for the dynamics trajectory.
    assertEquals(info + "End kinetic energy for OpenMM Langevin(Stochastic) integrator",
        endKineticEnergy, molDynOpenMM.getKineticEnergy(), tolerance);
    assertEquals(info + "End potential energy for OpenMM Langevin(Stochastic) integrator",
        endPotentialEnergy, molDynOpenMM.getPotentialEnergy(), tolerance);
    assertEquals(info + "End total energy for OpenMM Langevin(Stochastic) integrator",
        endTotalEnergy, molDynOpenMM.getTotalEnergy(), tolerance);
  }
}
