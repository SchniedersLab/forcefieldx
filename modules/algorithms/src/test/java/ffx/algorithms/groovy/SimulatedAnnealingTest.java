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

import ffx.algorithms.misc.AlgorithmsTest;
import ffx.algorithms.optimize.anneal.SimulatedAnnealing;
import java.util.Arrays;
import java.util.Collection;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/** @author Hernan Bernabe */
@RunWith(Parameterized.class)
public class SimulatedAnnealingTest extends AlgorithmsTest {

  private static final double tolerance = 0.01;
  private String info;
  private String filename;
  private double endKineticEnergy;
  private double endPotentialEnergy;
  private double endTotalEnergy;
  private double endTemperature;

  public SimulatedAnnealingTest(
      String info,
      String filename,
      double endKineticEnergy,
      double endPotentialEnergy,
      double endTotalEnergy,
      double endTemperature) {
    this.info = info;
    this.filename = filename;
    this.endKineticEnergy = endKineticEnergy;
    this.endPotentialEnergy = endPotentialEnergy;
    this.endTotalEnergy = endTotalEnergy;
    this.endTemperature = endTemperature;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "Acetamide with Stochastic integrator for Simulated Annealing Test", // info
                "ffx/algorithms/structures/acetamide_annealing.xyz", // filename
                2.7348819613719866, // endKineticEnergy
                -30.162204159002336, // endPotentialEnergy
                -27.427322197628616, // endTotalEnergy
                101.94414998339263 // endTemperature
            }
        });
  }

  @Test
  public void testSimulatedAnnealing() {

    // Set-up the input arguments for the script.
    String[] args = {
        "-n", "200",
        "-i", "STOCHASTIC",
        "-r", "0.001",
        "--tl", "100",
        "--tu", "400",
        "--tmS", "EXP",
        "-W", "10",
        "-w", "100",
        "-k", "100",
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);
    binding.setVariable("baseDir", registerTemporaryDirectory().toFile());

    Anneal anneal = new Anneal(binding).run();
    algorithmsScript = anneal;

    SimulatedAnnealing simulatedAnnealing = anneal.getAnnealing();

    // Assert that end energies are within the tolerance for the dynamics trajectory
    assertEquals(
        info + " Final kinetic energy",
        endKineticEnergy,
        simulatedAnnealing.getKineticEnergy(),
        tolerance);
    assertEquals(
        info + " Final potential energy",
        endPotentialEnergy,
        simulatedAnnealing.getPotentialEnergy(),
        tolerance);
    assertEquals(
        info + " Final total energy",
        endTotalEnergy,
        simulatedAnnealing.getTotalEnergy(),
        tolerance);
    assertEquals(
        info + " Final temperature",
        endTemperature,
        simulatedAnnealing.getTemperature(),
        tolerance);
  }
}
