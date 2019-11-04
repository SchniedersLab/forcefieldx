//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
//******************************************************************************
package ffx.algorithms.groovy;

import java.util.Arrays;
import java.util.Collection;

import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.junit.Assert.assertEquals;

import ffx.algorithms.dynamics.MolecularDynamicsOpenMM;
import ffx.algorithms.misc.PJDependentTest;
import ffx.algorithms.groovy.Dynamics;

import groovy.lang.Binding;

/**
 * @author Hernan V Bernabe
 */
@RunWith(Parameterized.class)
public class DynamicsOpenMMRESPANVETest extends PJDependentTest {

    private String info;
    private String filename;
    private double startingTotalEnergy;
    private double totalEnergyTolerance = 5.0;
    private boolean ffxOpenMM;

    private Binding binding;
    private Dynamics dynamics;

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "System OpenMM RESPA NVE", // info
                        "ffx/algorithms/structures/waterbox_eq.xyz", // filename
                        -24936.9565 // startingTotalEnergy
                }

        });
    }

    public DynamicsOpenMMRESPANVETest(String info, String filename, double startingTotalEnergy) {
        this.info = info;
        this.filename = filename;
        this.startingTotalEnergy = startingTotalEnergy;
    }

    @Before
    public void before() {
        binding = new Binding();
        dynamics = new Dynamics();
        dynamics.setBinding(binding);
    }

    @BeforeClass
    public static void beforeClass() {
        System.setProperty("platform", "omm");
        PJDependentTest.beforeClass();
    }

    @AfterClass
    public static void afterClass() {
        System.clearProperty("platform");
        PJDependentTest.afterClass();
    }

    @Test
    public void testDynamicsOpenMMRESPANVE() {
        if (!ffxOpenMM) {
            return;
        }

        // Set-up the input arguments for the script.
        String[] args = {"-n", "10", "-z", "1", "-t", "298.15", "-i", "RESPA", "-b", "Adiabatic", "-r", "0.001", "--mdE", "OpenMM", "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate script.
        dynamics.run();

        MolecularDynamicsOpenMM molDyn = (MolecularDynamicsOpenMM) dynamics.getMolecularDynamics();

        // Assert that the end total energy is within the threshold for the dynamics trajectory.
        assertEquals(info + "End total energy for OpenMM RESPA integrator under the NVE ensemble", startingTotalEnergy, molDyn.getTotalEnergy(), totalEnergyTolerance);
    }

}
