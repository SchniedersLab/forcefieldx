/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential;

import static org.junit.Assert.*;

import java.io.File;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import java.util.logging.Logger;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.utilities.Keyword;

/**
 * Test the PotentialEnergy class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class GradientTest {

    private static final Logger logger = Logger.getLogger(GradientTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                    {
                        false,
                        "Water Tiny (81 water molecules",
                        "ffx/potential/structures/watertiny.xyz"}});
    }
    private final File structure;
    private final MolecularAssembly molecularAssembly;
    private final ForceFieldEnergy energy;
    private final double tolerance = 1.0e-3;
    private final double gradientTolerance = 1.0e-4;
    private final boolean ci;

    public GradientTest(
            boolean ciOnly,
            String info, String filename) {

        String ffxCi = System.getProperty("ffx.ci", "false");
        if (ffxCi.equalsIgnoreCase("true")) {
            ci = true;
        } else {
            ci = false;
        }

        if (!ci && ciOnly) {
            structure = null;
            molecularAssembly = null;
            energy = null;
            return;
        }

        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());

        String name = structure.getName();
        int index = filename.lastIndexOf(".");
        name = filename.substring(0, index);
        molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);

        CompositeConfiguration properties = Keyword.loadProperties(structure);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);

        XYZFilter xyzFilter = new XYZFilter(structure, molecularAssembly, forceField, properties);
        boolean expectedReturn = true;
        boolean actualReturn = xyzFilter.readFile();
        assertEquals(info, expectedReturn, actualReturn);
        Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
        molecularAssembly.finalize(true);
        energy = new ForceFieldEnergy(molecularAssembly);
    }

    /**
     * Test of energy gradient, of class PotentialEnergy.
     */
    @Test
    public void testGradient() {
        boolean gradient = true;
        boolean print = true;
        energy.energy(gradient, print);
        gradient = false;
        print = true;
        double step = 0.00001;
        double analytic[] = new double[3];
        double numeric[] = new double[3];
        double xyz[] = new double[3];

        Atom[] atoms = molecularAssembly.getAtomArray();
        int n = atoms.length;

        for (int i = 0; i < n; i++) {
            Atom a0 = atoms[i];
            a0.getXYZGradient(analytic);
            a0.getXYZ(xyz);
            // Find numeric dX
            xyz[0] += step;
            a0.moveTo(xyz);
            double e = energy.energy(gradient, print);
            xyz[0] -= 2.0 * step;
            a0.moveTo(xyz);
            e -= energy.energy(gradient, print);
            numeric[0] = e / (2.0 * step);
            xyz[0] += step;
            // Find numeric dY
            xyz[1] += step;
            a0.moveTo(xyz);
            e = energy.energy(gradient, print);
            xyz[1] -= 2.0 * step;
            a0.moveTo(xyz);
            e -= energy.energy(gradient, print);
            numeric[1] = e / (2.0 * step);
            xyz[1] += step;
            // Find numeric dZ
            xyz[2] += step;
            a0.moveTo(xyz);
            e = energy.energy(gradient, print);
            xyz[2] -= 2.0 * step;
            a0.moveTo(xyz);
            e -= energy.energy(gradient, print);
            numeric[2] = e / (2.0 * step);
            xyz[2] += step;
            a0.moveTo(xyz);
            double dx = analytic[0] - numeric[0];
            double dy = analytic[1] - numeric[1];
            double dz = analytic[2] - numeric[2];
            double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (len > gradientTolerance) {
                logger.info("\n" + a0.toString() + String.format(" failed: %10.6f.", len) + String.format(
                        "\nAnalytic: (%12.4f, %12.4f, %12.4f)\n",
                        analytic[0], analytic[1], analytic[2]) + String.format("Numeric:  (%12.4f, %12.4f, %12.4f)\n",
                                                                               numeric[0], numeric[1], numeric[2]));
            } else {
                logger.info("\n" + a0.toString() + String.format(" passed: %10.6f.", len) + String.format(
                        "\nAnalytic: (%12.4f, %12.4f, %12.4f)\n",
                        analytic[0], analytic[1], analytic[2]) + String.format("Numeric:  (%12.4f, %12.4f, %12.4f)\n",
                                                                               numeric[0], numeric[1], numeric[2]));
            }
            assertEquals(a0.toString(), 0.0, len, gradientTolerance);
        }
    }

    public void testSoftCore() {
        boolean gradient = false;
        boolean print = true;
        double e = energy.energy(gradient, print);
        Atom atoms[] = molecularAssembly.getAtomArray();
        int n = atoms.length;
        // Make the 3 atoms soft
        for (int i = n; i > n - 3; i--) {
            Atom ai = atoms[i - 1];
            ai.setApplyLambda(true);
        }
        // Compute the energy with Lambda = 1.0;
        double lambda = 1.0;
        energy.setLambda(lambda);
        double e2 = energy.energy(gradient, print);
        assertEquals(e, e2, tolerance);
    }
}
