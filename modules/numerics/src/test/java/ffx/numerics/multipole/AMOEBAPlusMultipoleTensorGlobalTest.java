// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import static ffx.numerics.multipole.MultipoleTensorTest.*;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.junit.Assert.assertEquals;

@RunWith(Parameterized.class)
public class AMOEBAPlusMultipoleTensorGlobalTest {
    /**
     * Logger for the MultipoleTensor class.
     */
    private static final Logger logger = Logger.getLogger(AMOEBAPlusMultipoleTensorGlobalTest.class.getName());
    private final double tolerance = 1.0e-13;
    private final double fdTolerance = 1.0e-6;
    private final double[] r = new double[3];
    private final double[] r2 = new double[3];
    private final int order;
    private final String info;
    private final Operator operator;

    public AMOEBAPlusMultipoleTensorGlobalTest(
            String info,
            int order,
            Operator operator) {
        this.info = info;
        this.order = order;
        r[0] = QkXYZ[0] - QiXYZ[0];
        r[1] = QkXYZ[1] - QiXYZ[1];
        r[2] = QkXYZ[2] - QiXYZ[2];
        r2[0] = QjXYZ[0] - QkXYZ[0];
        r2[1] = QjXYZ[1] - QkXYZ[1];
        r2[2] = QjXYZ[2] - QkXYZ[2];
        this.operator = operator;
    }

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(
                new Object[][]{
                        {"Order 5 Aplus Direct", 5, Operator.AMOEBA_PLUS_OVERLAP_FIELD}
                });
    }

    @Test
    public void permanentMultipoleEnergyAndGradTest() {
        double delta = 1.0e-8;
        double delta2 = 2.0 * delta;
        double[] Fi = new double[3];
        double[] Fk = new double[3];
        double[] Ti = new double[3];
        double[] Tk = new double[3];

        // Oxygen-Oxygen Mpole-Mpole interaction
        PolarizableMultipole mI = new PolarizableMultipole(QiAmoebaP, QiInduced, QiInduced, Zi);
        PolarizableMultipole mK = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        MultipoleTensor multipoleTensor = new AmoebaPlusOverlapTensorGlobal(order, chargePenAlphaOx, chargePenAlphaOx);
        multipoleTensor.generateTensor(r);
        double eOverlap = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);

        // Analytic gradient.
        double aX = Fk[0];
        double aY = Fk[1];
        double aZ = Fk[2];

        r[0] += delta;
        multipoleTensor.generateTensor(r);
        double posX = multipoleTensor.multipoleEnergy(mI, mK);
        r[0] -= delta2;
        multipoleTensor.generateTensor(r);
        double negX = multipoleTensor.multipoleEnergy(mI, mK);
        r[0] += delta;
        r[1] += delta;
        multipoleTensor.generateTensor(r);
        double posY = multipoleTensor.multipoleEnergy(mI, mK);
        r[1] -= delta2;
        multipoleTensor.generateTensor(r);
        double negY = multipoleTensor.multipoleEnergy(mI, mK);
        r[1] += delta;
        r[2] += delta;
        multipoleTensor.generateTensor(r);
        double posZ = multipoleTensor.multipoleEnergy(mI, mK);
        r[2] -= delta2;
        multipoleTensor.generateTensor(r);
        double negZ = multipoleTensor.multipoleEnergy(mI, mK);
        r[2] += delta;

        double expect = aX;
        double actual = (posX - negX) / delta2;
        assertEquals(info + " OO Overlap Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OO Overlap Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OO Overlap Force Z", expect, actual, fdTolerance);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);

        // Oxygen-Oxygen Core-Multipole interaction
        AmoebaPlusDampTensorGlobal multipoleTensorDamp = new AmoebaPlusDampTensorGlobal(3, chargePenAlphaOx, chargePenAlphaOx);
        multipoleTensorDamp.generateTensor(r);
        double eProton = multipoleTensorDamp.coreInteractionAndGradient(mI, mK, Fi, Fk);
        double e = eOverlap + eProton;
        assertEquals(" Oxygen-Oxygen Perm Multipole Coulomb Energy", amoebaPlusMPoleEnergyOO, e, tolerance);

        // Analytic gradient.
        aX = Fk[0];
        aY = Fk[1];
        aZ = Fk[2];

        r[0] += delta;
        multipoleTensorDamp.generateTensor(r);
        posX = multipoleTensorDamp.coreInteraction(mI, mK);
        r[0] -= delta2;
        multipoleTensorDamp.generateTensor(r);
        negX = multipoleTensorDamp.coreInteraction(mI, mK);
        r[0] += delta;
        r[1] += delta;
        multipoleTensorDamp.generateTensor(r);
        posY = multipoleTensorDamp.coreInteraction(mI, mK);
        r[1] -= delta2;
        multipoleTensorDamp.generateTensor(r);
        negY = multipoleTensorDamp.coreInteraction(mI, mK);
        r[1] += delta;
        r[2] += delta;
        multipoleTensorDamp.generateTensor(r);
        posZ = multipoleTensorDamp.coreInteraction(mI, mK);
        r[2] -= delta2;
        multipoleTensorDamp.generateTensor(r);
        negZ = multipoleTensorDamp.coreInteraction(mI, mK);
        r[2] += delta;

        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " OO Core Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OO Core Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OO Core Force Z", expect, actual, fdTolerance);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);

        // Hydrogen-Oxygen Mpole-Mpole interaction
        mI = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        mK = new PolarizableMultipole(QjAmoebaP, QjInduced, QjInduced, Zj);
        multipoleTensor = new AmoebaPlusOverlapTensorGlobal(order, chargePenAlphaOx, chargePenAlphaHyd);
        multipoleTensor.generateTensor(r2);
        eOverlap = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);

        // Analytic gradient.
        aX = Fk[0];
        aY = Fk[1];
        aZ = Fk[2];

        r2[0] += delta;
        multipoleTensor.generateTensor(r2);
        posX = multipoleTensor.multipoleEnergy(mI, mK);
        r2[0] -= delta2;
        multipoleTensor.generateTensor(r2);
        negX = multipoleTensor.multipoleEnergy(mI, mK);
        r2[0] += delta;
        r2[1] += delta;
        multipoleTensor.generateTensor(r2);
        posY = multipoleTensor.multipoleEnergy(mI, mK);
        r2[1] -= delta2;
        multipoleTensor.generateTensor(r2);
        negY = multipoleTensor.multipoleEnergy(mI, mK);
        r2[1] += delta;
        r2[2] += delta;
        multipoleTensor.generateTensor(r2);
        posZ = multipoleTensor.multipoleEnergy(mI, mK);
        r2[2] -= delta2;
        multipoleTensor.generateTensor(r2);
        negZ = multipoleTensor.multipoleEnergy(mI, mK);
        r2[2] += delta;

        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " OH Overlap Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OH Overlap Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OH Overlap Force Z", expect, actual, fdTolerance);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);

        // Hydrogen-Oxygen Core-Multipole interaction
        // Order of args matter here first alpha -> mI, second alpha -> mK
        multipoleTensorDamp = new AmoebaPlusDampTensorGlobal(3, chargePenAlphaOx, chargePenAlphaHyd);
        multipoleTensorDamp.generateTensor(r2);
        // Order of args matter here first alpha -> mI, second alpha -> mK
        eProton = multipoleTensorDamp.coreInteractionAndGradient(mI, mK, Fi, Fk);
        e = eOverlap + eProton;
        assertEquals(" Hydrogen-Oxygen Perm Multipole Damped Coulomb Energy", amoebaPlusMPoleEnergyOH, e, tolerance);

        // Analytic gradient.
        aX = Fk[0];
        aY = Fk[1];
        aZ = Fk[2];

        r2[0] += delta;
        multipoleTensorDamp.generateTensor(r2);
        posX = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[0] -= delta2;
        multipoleTensorDamp.generateTensor(r2);
        negX = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[0] += delta;
        r2[1] += delta;
        multipoleTensorDamp.generateTensor(r2);
        posY = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[1] -= delta2;
        multipoleTensorDamp.generateTensor(r2);
        negY = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[1] += delta;
        r2[2] += delta;
        multipoleTensorDamp.generateTensor(r2);
        posZ = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[2] -= delta2;
        multipoleTensorDamp.generateTensor(r2);
        negZ = multipoleTensorDamp.coreInteraction(mI, mK);
        r2[2] += delta;

        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " OH Core Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OH Core Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OH Core Force Z", expect, actual, fdTolerance);
    }

    @Test
    public void inducedDipoleEnergyAndGradient() {
        double[] Fi = new double[3];
        double[] Ti = new double[3];
        double[] Tk = new double[3];
        double thole = .7;
        double aiak = 1.0 / (0.99595940317287035 * 0.99595940317287035);

        // Oxygen-Oxygen
        PolarizableMultipole mI = new PolarizableMultipole(QiAmoebaP, QiInduced, QiInduced, Zi);
        PolarizableMultipole mK = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        MultipoleTensor multipoleTensor = new TholeTensorGlobal(4, thole, aiak, true);
        multipoleTensor.generateTensor(r);
        // This test is for a Thole-direct induced dipole-multipole interaction
        double e = multipoleTensor.polarizationEnergyAndGradient(mI, mK, scaleMutual,
                scaleMutual, 0.0, Fi, Ti, Tk);
        assertEquals(" Oxygen-Oxygen Thole-Damped Polarization Energy", amoebaPlusIndDipoleEnergyOO, e, tolerance);
        double[] aplusIndGrad = new double[]
                {-7.8921873374791535E-004,
                        -6.3709119990021249E-004,
                        -2.4873367831826638E-021};
        assertEquals(info + " Polarization GradX", aplusIndGrad[0], Fi[0], tolerance);
        assertEquals(info + " Polarization GradY", aplusIndGrad[1], Fi[1], tolerance);
        assertEquals(info + " Polarization GradZ", aplusIndGrad[2], Fi[2], tolerance);
    }
}