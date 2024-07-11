// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
        r[0] = QkXYZ[0]-QiXYZ[0];
        r[1] = QkXYZ[1]-QiXYZ[1];
        r[2] = QkXYZ[2]-QiXYZ[2];
        r2[0] = QjXYZ[0]-QkXYZ[0];
        r2[1] = QjXYZ[1]-QkXYZ[1];
        r2[2] = QjXYZ[2]-QkXYZ[2];
        this.operator = operator;
    }

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(
                new Object[][]{
                        {"Order 5 Coulomb", 5, Operator.AMOEBA_PLUS_OVERLAP_FIELD}
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

        // Oxygen-Oxygen
        PolarizableMultipole mI = new PolarizableMultipole(QiAmoebaP, QiInduced, QiInduced, Zi);
        PolarizableMultipole mK = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        MultipoleTensor multipoleTensor = new AmoebaPlusOverlapTensorGlobal(order, chargePenAlphaOx, chargePenAlphaOx);
        multipoleTensor.generateTensor(r);
        double eOverlap = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);
        multipoleTensor = new AmoebaPlusDampTensorGlobal(2, chargePenAlphaOx, chargePenAlphaOx);
        multipoleTensor.generateTensor(r);
        double eProton = multipoleTensor.coreInteractionAndGradient(mI, mK, Fi, Fk);
        double e = eOverlap + eProton;
        assertEquals(" Oxygen-Oxygen Perm Multipole Coulomb Energy", amoebaPlusMPoleEnergyOO, e, tolerance);

        // Analytic gradient.
        double aX = Fk[0];
        double aY = Fk[1];
        double aZ = Fk[2];

        r[0] += delta;
        multipoleTensor.generateTensor(r);
        double posX = multipoleTensor.coreInteraction(mI, mK);
        r[0] -= delta2;
        multipoleTensor.generateTensor(r);
        double negX = multipoleTensor.coreInteraction(mI, mK);
        r[0] += delta;
        r[1] += delta;
        multipoleTensor.generateTensor(r);
        double posY = multipoleTensor.coreInteraction(mI, mK);
        r[1] -= delta2;
        multipoleTensor.generateTensor(r);
        double negY = multipoleTensor.coreInteraction(mI, mK);
        r[1] += delta;
        r[2] += delta;
        multipoleTensor.generateTensor(r);
        double posZ = multipoleTensor.coreInteraction(mI, mK);
        r[2] -= delta2;
        multipoleTensor.generateTensor(r);
        double negZ = multipoleTensor.coreInteraction(mI, mK);
        r[2] += delta;

        double expect = aX;
        double actual = (posX - negX) / delta2;
        assertEquals(info + " Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " Force Z", expect, actual, fdTolerance);

        // Hydrogen-Oxygen
        mI = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        mK = new PolarizableMultipole(QjAmoebaP, QjInduced, QjInduced, Zj);
        multipoleTensor = new AmoebaPlusOverlapTensorGlobal(order, chargePenAlphaOx, chargePenAlphaHyd);
        multipoleTensor.generateTensor(r2);
        eOverlap = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);
        // Order of args matter here first alpha -> mI, second alpha -> mK
        multipoleTensor = new AmoebaPlusDampTensorGlobal(2, chargePenAlphaOx, chargePenAlphaHyd);
        multipoleTensor.generateTensor(r2);
        // Order of args matter here first alpha -> mI, second alpha -> mK
        eProton = multipoleTensor.coreInteractionAndGradient(mI, mK, Fi, Fk);
        e = eOverlap + eProton;
        assertEquals(" Hydrogen-Oxygen Perm Multipole Damped Coulomb Energy", amoebaPlusMPoleEnergyOH, e, tolerance);

    }

    @Test
    public void inducedDipoleEnergyAndGradient(){
        // Oxygen-Oxygen
        PolarizableMultipole mI = new PolarizableMultipole(QiAmoebaP, QiInduced, QiInduced, Zi);
        PolarizableMultipole mK = new PolarizableMultipole(QkAmoebaP, QkInduced, QkInduced, Zk);
        double thole = .7;
        double aiak = 1.0 / (0.99595940317287035 * 0.99595940317287035);
        MultipoleTensor multipoleTensor = new TholeTensorGlobal(3, thole, aiak, true);
        multipoleTensor.generateTensor(r);
        double e = multipoleTensor.polarizationEnergy(mI, mK, scaleMutual);
        assertEquals(" Oxygen-Oxygen Thole-Damped Polarization Energy", amoebaPlusIndDipoleEnergyOO, e, tolerance);
    }
}
