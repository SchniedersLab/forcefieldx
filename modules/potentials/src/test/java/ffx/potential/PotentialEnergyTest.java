/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
public class PotentialEnergyTest {

    private static final Logger logger = Logger.getLogger(PotentialEnergyTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"DHFR Benchmark",
                        "ffx/potential/structures/dhfr.xyz",
                        6423.84579926e0, 16569,
                        3749.73436272e0, 11584,
                        -21.85553039e0, 4031,
                        687.46861123e0, 7023,
                        198.72886589e0, 1566,
                        422.54712638e0, 6701,
                        48.26628393e0, 292,
                        -41.71473466e0, 147,
                        32676.57816012e0, 3480445,
                        -79670.54372856e0, 1463353,
                        -32033.24720761e0, 1463353} /*,
                {"SNARE P1",
                "ffx/potential/structures/1n7s.P1.xyz",
                1405.28569930e0, 20160,
                2976.77005458e0, 33020,
                25.68062976e0, 27696,
                10.00655326e0, 1288,
                541.44132652e0, 9948,
                1671.56977674e0, 45796,
                159.42575736e0, 1480,
                -2243.98305878e0, 1072,
                16013.08734189e0, 2966542,
                -49214.99281245e0, 1504176,
                -11244.77107599e0, 1504176},
                {"SNARE P212121",
                "ffx/potential/structures/1n7s.P212121.xyz",
                351.321424825e0, 5040,
                744.192513645e0, 8255,
                6.42015744e0, 6924,
                2.501638315e0, 322,
                135.36033163e0, 2487,
                417.892444185e0, 11449,
                39.85643934e0, 370,
                -560.995764695e0, 268,
                4003.2718354725e0, 2966542,
                -12303.7482031125e0, 1504176,
                -2811.1927689975e0, 1504176 } */});
    }
    private final String info;
    private final File structure;
    private final MolecularAssembly molecularAssembly;
    private final int nBonds;
    private final int nAngles;
    private final int nStretchBends;
    private final int nUreyBradleys;
    private final int nOutOfPlaneBends;
    private final int nTorsions;
    private final int nPiOrbitalTorsions;
    private final int nTorsionTorsions;
    private final int nVanDerWaals;
    private final int nPermanent;
    private final int nPolar;
    private double bondEnergy;
    private double angleEnergy;
    private double stretchBendEnergy;
    private double ureyBradleyEnergy;
    private double outOfPlaneBendEnergy;
    private double torsionEnergy;
    private double piOrbitalTorsionEnergy;
    private double torsionTorsionEnergy;
    private double vanDerWaalsEnergy;
    private double permanentEnergy;
    private Polarization polarization;
    private double polarizationEnergy;
    private final PotentialEnergy energy;
    private boolean mpoleTerm;
    private final double tolerance = 1.0e-3;
    private final double gradientTolerance = 1.0e-4;

    public PotentialEnergyTest(String info, String filename,
            double bondEnergy, int nBonds,
            double angleEnergy, int nAngles,
            double stretchBendEnergy, int nStretchBends,
            double ureyBradleyEnergy, int nUreyBradleys,
            double outOfPlaneBendEnergy, int nOutOfPlaneBends,
            double torsionEnergy, int nTorsions,
            double piOrbitalTorsionEnergy, int nPiOrbitalTorsions,
            double torsionTorsionEnergy, int nTorsionTorsions,
            double vanDerWaalsEnergy, int nVanDerWaals,
            double permanentEnergy, int nPermanent,
            double polarizationEnergy, int nPolar) {
        this.info = info;
        this.bondEnergy = bondEnergy;
        this.nBonds = nBonds;
        this.angleEnergy = angleEnergy;
        this.nAngles = nAngles;
        this.stretchBendEnergy = stretchBendEnergy;
        this.nStretchBends = nStretchBends;
        this.ureyBradleyEnergy = ureyBradleyEnergy;
        this.nUreyBradleys = nUreyBradleys;
        this.outOfPlaneBendEnergy = outOfPlaneBendEnergy;
        this.nOutOfPlaneBends = nOutOfPlaneBends;
        this.torsionEnergy = torsionEnergy;
        this.nTorsions = nTorsions;
        this.piOrbitalTorsionEnergy = piOrbitalTorsionEnergy;
        this.nPiOrbitalTorsions = nPiOrbitalTorsions;
        this.torsionTorsionEnergy = torsionTorsionEnergy;
        this.nTorsionTorsions = nTorsionTorsions;
        this.vanDerWaalsEnergy = vanDerWaalsEnergy;
        this.nVanDerWaals = nVanDerWaals;
        this.permanentEnergy = permanentEnergy;
        this.nPermanent = nPermanent;
        this.polarizationEnergy = polarizationEnergy;
        this.nPolar = nPolar;

        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());

        String name = structure.getName();
        int index = filename.lastIndexOf(".");
        name = filename.substring(0, index);
        molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);

        CompositeConfiguration properties = Keyword.loadProperties(structure);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties, null);
        ForceField forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);

        XYZFilter xyzFilter = new XYZFilter(molecularAssembly, forceField);
        boolean expectedReturn = true;
        boolean actualReturn = xyzFilter.readFile();
        assertEquals(info, expectedReturn, actualReturn);
        Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
        molecularAssembly.finalize(true);

        energy = new PotentialEnergy(molecularAssembly);
        mpoleTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, true);

        String polar = forceField.getString(ForceFieldString.POLARIZATION,
                "MUTUAL");
        if (polar.equalsIgnoreCase("MUTUAL")) {
            polarization = Polarization.MUTUAL;
        } else if (polar.equalsIgnoreCase("DIRECT")) {
            polarization = Polarization.DIRECT;
        } else {
            polarization = Polarization.NONE;
        }
    }

    /**
     * Test of energy method, of class PotentialEnergy.
     */
    @Test
    public void testEnergy() {
        boolean gradient = false;
        boolean print = true;
        double total = energy.energy(gradient, print);
        // Bond Energy
        assertEquals(info + " Bond Energy", bondEnergy, energy.bondEnergy, tolerance);
        assertEquals(info + " Bond Count", nBonds, energy.nBonds);
        // Angle Energy
        assertEquals(info + " Angle Energy", angleEnergy, energy.angleEnergy, tolerance);
        assertEquals(info + " Angle Count", nAngles, energy.nAngles);
        // Stretch-Bend Energy
        assertEquals(info + " Stretch-Bend Energy", stretchBendEnergy, energy.stretchBendEnergy, tolerance);
        assertEquals(info + " Stretch-Bend Count", nStretchBends, energy.nStretchBends);
        // Urey-Bradley Energy
        assertEquals(info + " Urey-Bradley Energy", ureyBradleyEnergy, energy.ureyBradleyEnergy, tolerance);
        assertEquals(info + " Urey-Bradley Count", nUreyBradleys, energy.nUreyBradleys);
        // Out-of-Plane Bend
        assertEquals(info + " Out-of-Plane Bend Energy", outOfPlaneBendEnergy, energy.outOfPlaneBendEnergy, tolerance);
        assertEquals(info + " Out-of-Plane Bend Count", nOutOfPlaneBends, energy.nOutOfPlaneBends);
        // Torsional Angle
        assertEquals(info + " Torsion Energy", torsionEnergy, energy.torsionEnergy, tolerance);
        assertEquals(info + " Torsion Count", nTorsions, energy.nTorsions);
        // Pi-Orbital Torsion
        assertEquals(info + " Pi-OrbitalTorsion Energy", piOrbitalTorsionEnergy, energy.piOrbitalTorsionEnergy, tolerance);
        assertEquals(info + " Pi-OrbitalTorsion Count", nPiOrbitalTorsions, energy.nPiOrbitalTorsions);
        // Torsion-Torsion
        assertEquals(info + " Torsion-Torsion Energy", torsionTorsionEnergy, energy.torsionTorsionEnergy, tolerance);
        assertEquals(info + " Torsion-Torsion Count", nTorsionTorsions, energy.nTorsionTorsions);
        // van Der Waals
        assertEquals(info + " van Der Waals Energy", vanDerWaalsEnergy, energy.vanDerWaalsEnergy, tolerance);
        //assertEquals(info + " van Der Waals Count", nVanDerWaals, energy.nVanDerWaals);
        // Permanent Multipoles
        if (mpoleTerm) {
            assertEquals(info + " Permanent Multipole Energy", permanentEnergy, energy.permanentMultipoleEnergy, tolerance);
        }
        //assertEquals(info + " van Der Waals Count", nVanDerWaals, energy.nVanDerWaals);
        // Polarization
        if (polarization == Polarization.MUTUAL) {
            assertEquals(info + " Polarization Energy", polarizationEnergy, energy.polarizationEnergy, tolerance);
            //assertEquals(info + " van Der Waals Count", nVanDerWaals, energy.nVanDerWaals);
        }
    }

    /**
     * Test of energy gradient, of class PotentialEnergy.
     */
    public void testGradient() {
        boolean gradient = true;
        boolean print = true;
        energy.energy(gradient, print);
        gradient = false;
        print = false;
        double step = 0.00001;
        double analytic[] = new double[3];
        double numeric[] = new double[3];
        double xyz[] = new double[3];

        ArrayList<Atom> atoms = molecularAssembly.getAtomList();
        int n = atoms.size();
        Random random = new Random();
        int i = random.nextInt(n);
        Atom a0 = atoms.get(i);

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
            logger.severe("\n" + a0.toString() + String.format(" failed: %10.6f.", len) + String.format(
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

    @Test
    public void testSoftCore() {

        Atom atoms[] = molecularAssembly.getAtomArray();
        int n = atoms.length;
        // Make the last water soft
        for (int i = n; i > n - 3; i--) {
            Atom ai = atoms[i - 1];
            ai.setSoftCore(true);
        }
        boolean gradient = false;
        boolean print = true;
        // Compute the energy with Lambda = 1.0;
        for (int i = 0; i <= 1; i++) {
            double lambda = 1.0 - i * 0.1;
            energy.setSoftCoreLambda(lambda);
            energy.energy(gradient, print);
        }
    }
}
