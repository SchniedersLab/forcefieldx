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
package ffx.potential.groovy;

import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;

/**
 * Test the Energy script with Generalized Kirkwood parameters.
 */
@RunWith(Parameterized.class)
public class GKEnergyTest extends ParentEnergyTest {

    public GKEnergyTest(String info, String filename, int nAtoms, double bondEnergy, int nBonds,
                        double angleEnergy, int nAngles, double stretchBendEnergy, int nStretchBends,
                        double ureyBradleyEnergy, int nUreyBradleys, double outOfPlaneBendEnergy,
                        int nOutOfPlaneBends, double torsionEnergy, int nTorsions, double improperTorsionEnergy,
                        int nImproperTorsions, double piOrbitalTorsionEnergy, int nPiOrbitalTorsions,
                        double torsionTorsionEnergy, int nTorsionTorsions, double stretchTorsionEnergy,
                        int nStretchTorsions, double angleTorsionEnergy, int nAngleTorsions, double vanDerWaalsEnergy,
                        int nVanDerWaals, double permanentEnergy, int nPermanent, double polarizationEnergy,
                        int nPolar, double gkEnergy, int nGK, boolean testOpenMM) {
        super(info, filename, nAtoms, bondEnergy, nBonds, angleEnergy, nAngles, stretchBendEnergy, nStretchBends,
                ureyBradleyEnergy, nUreyBradleys, outOfPlaneBendEnergy, nOutOfPlaneBends, torsionEnergy,
                nTorsions, improperTorsionEnergy, nImproperTorsions, piOrbitalTorsionEnergy,
                nPiOrbitalTorsions, torsionTorsionEnergy, nTorsionTorsions, stretchTorsionEnergy,
                nStretchTorsions, angleTorsionEnergy, nAngleTorsions, vanDerWaalsEnergy, nVanDerWaals,
                permanentEnergy, nPermanent, polarizationEnergy, nPolar, gkEnergy, nGK, testOpenMM);
    }

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(
                new Object[][]{
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "RNA (1mis) with GK from Corrigan et al (2023)",
                                "ffx/potential/structures/1mis.xyz",
                                522, // Atoms
                                60.88950457, // Bond
                                562,
                                136.35249695, // Angle
                                998,
                                2.43194433, // Stretch-Bend
                                848,
                                0.0, // Urey-Bradley
                                0,
                                0.10282386, // Out-of-Plane
                                294,
                                79.48824698, // Torsion
                                1492,
                                0.0, // Improper Torsion
                                0,
                                0.45934575, // Pi-Orbital Torsion
                                96,
                                0.0,  // Torsion-Torsion
                                0,
                                -3.02691282, // Stretch-Torsion
                                44,
                                -4.54420615, // Angle-Torsion
                                72,
                                218.31775844, // vdW
                                134421,
                                875.03795766, // Permanent
                                134421,
                                -238.72680230, // Polarization
                                134421,
                                -3236.93867799, // Total Solvation
                                // Generalized Kirkwood -3221.37466819
                                // Cavitation 272.96015049, Dispersion -288.52416028
                                136503,
                                false
                        },
                        {
                                "RNA (1mis) with GK from Corrigan et al (2021)",
                                "ffx/potential/structures/1mis.gk2021.xyz",
                                522, // Atoms
                                60.88950457, // Bond
                                562,
                                136.35249695, // Angle
                                998,
                                2.43194433, // Stretch-Bend
                                848,
                                0.0, // Urey-Bradley
                                0,
                                0.10282386, // Out-of-Plane
                                294,
                                79.48824698, // Torsion
                                1492,
                                0.0, // Improper Torsion
                                0,
                                0.45934575, // Pi-Orbital Torsion
                                96,
                                0.0,  // Torsion-Torsion
                                0,
                                -3.02691282, // Stretch-Torsion
                                44,
                                -4.54420615, // Angle-Torsion
                                72,
                                218.31775844, // vdW
                                134421,
                                875.03795766, // Permanent
                                134421,
                                -238.72700000, // Polarization
                                134421,
                                -3308.29464032, // Total Solvation
                                136503,
                                false
                        },
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2023)",
                                "ffx/potential/structures/1bpi.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.88505475, // Permanent
                                394854,
                                -242.15918841, // Polarization
                                394854,
                                -923.17087281, // Total Solvation
                                //Generalized Kirkwood -970.64995804
                                //Cavitation 444.79997963, Dispersion -397.32089439
                                398278,
                                false
                        },
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2023) - No Neck",
                                "ffx/potential/structures/1bpi.noneck.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.88505475, // Permanent
                                394854,
                                -242.15918841, // Polarization
                                394854,
                                -900.45662762, // Total Solvation
                                //Generalized Kirkwood -947.93571285
                                //Cavitation 444.79997963, Dispersion -397.32089439
                                398278,
                                false
                        },
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2023) - No Tanh",
                                "ffx/potential/structures/1bpi.notanh.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.88505475, // Permanent
                                394854,
                                -242.15918841, // Polarization
                                394854,
                                -1241.80712005, // Total Solvation
                                //Generalized Kirkwood  -1289.28620528
                                //Cavitation 444.79997963, Dispersion -397.32089439
                                398278,
                                false
                        },
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2023) - No Neck, No Tanh",
                                "ffx/potential/structures/1bpi.noneck.notanh.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.88505475, // Permanent
                                394854,
                                -242.15918841, // Polarization
                                394854,
                                -1325.30614458, // Total Solvation
                                //Generalized Kirkwood  -1372.78522981
                                //Cavitation 444.79997963, Dispersion -397.32089439
                                398278,
                                false
                        },
                        /*
                         * Test values for Solvation/GK/Atomic Multipoles/Polarization don't exactly match
                         * published values (Corrigan et. al. 2023) due to small multipole updates to the
                         * AMOEBA-BIO-2018 force field made in early 2023 to accommodate CpHMD work
                         *
                         * The error magnitude is <2 kcal for each term
                        */
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2023) - No Element HCT or Interstitial Space Corrections",
                                "ffx/potential/structures/1bpi.noelemhct.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.88505475, // Permanent
                                394854,
                                -242.15918841, // Polarization
                                394854,
                                -1330.12099084, // Total Solvation
                                //Generalized Kirkwood  -1377.60007607
                                //Cavitation 444.79997963, Dispersion -397.32089439
                                398278,
                                false
                        },
                        {
                                "Protein (1bpi) with GK from Corrigan et al (2021)",
                                "ffx/potential/structures/1bpi.gk2021.xyz",
                                892, // Atoms
                                281.75870048, // Bond
                                906,
                                235.11731039, // Angle
                                1626,
                                -12.11600646, // Stretch-Bend
                                1455,
                                0.0, // Urey-Bradley
                                0,
                                28.91012526, // Out-of-Plane
                                597,
                                69.15283653, // Torsion
                                2391,
                                0.0, // Improper Torsion
                                0,
                                27.49698981, // Pi-Orbital Torsion
                                109,
                                -36.43950083,  // Torsion-Torsion
                                6,
                                0.0, // Stretch-Torsion
                                0,
                                0.0, // Angle-Torsion
                                0,
                                470.1728659015, // vdW
                                394854,
                                -1189.885054752826, // Permanent
                                394854,
                                -242.15918841434635, // Polarization
                                394854,
                                -1166.1124974225615, // Total Solvation
                                398278,
                                false
                        }

                });
    }


}
