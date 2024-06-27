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

import jdk.incubator.vector.DoubleVector;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;

import static ffx.numerics.multipole.MultipoleTensorTest.Qi;
import static ffx.numerics.multipole.MultipoleTensorTest.Qk;
import static ffx.numerics.multipole.MultipoleTensorTest.Ui;
import static ffx.numerics.multipole.MultipoleTensorTest.UiEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.Uk;
import static ffx.numerics.multipole.MultipoleTensorTest.UkEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueI;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueK;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueKEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permanentEnergy;
import static ffx.numerics.multipole.MultipoleTensorTest.permanentEnergyEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradICoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradIThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueICoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueIThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKCoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyCoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyThole;
import static ffx.numerics.multipole.MultipoleTensorTest.scaleMutual;
import static ffx.numerics.multipole.MultipoleTensorTest.xyz;
import static java.lang.String.format;
import static org.junit.Assert.assertEquals;

/**
 * Parameterized Test of the MultipoleTensorSIMD class based on previous MultipoleTensorSIMD tests.
 *
 * @author Matthew Speranza
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class MultipoleTensorGlobalSIMDTest {
    private final double tolerance = 1.0e-13;
    private final double fdTolerance = 1.0e-6;
    private final double[] r = new double[3];
    private DoubleVector xVector;
    private DoubleVector yVector;
    private DoubleVector zVector;
    private DoubleVector[] Fi;
    private DoubleVector[] Fk;
    private DoubleVector[] Ti;
    private DoubleVector[] Tk;
    private double[][] qI;
    private double[][] uI;
    private double[][] qK;
    private double[][] uK;
    private final int order;
    private final int tensorCount;
    private final String info;
    private final Operator operator;
    private final double beta = 0.545;
    private final double thole = 0.39;
    private final double AiAk = 1.061104559485911;

    public MultipoleTensorGlobalSIMDTest(
            String info,
            int order,
            Operator operator) {
        this.info = info;
        this.order = order;
        r[0] = MultipoleTensorTest.xyz[0];
        r[1] = MultipoleTensorTest.xyz[1];
        r[2] = MultipoleTensorTest.xyz[2];
        xVector = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, r[0]);
        yVector = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, r[1]);
        zVector = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, r[2]);
        this.tensorCount = MultipoleUtilities.tensorCount(order);
        this.operator = operator;

        // Fill in the arrays with the SIMD values.
        qI = new double[10][DoubleVector.SPECIES_PREFERRED.length()];
        uI = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        qK = new double[10][DoubleVector.SPECIES_PREFERRED.length()];
        uK = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        for (int i = 0; i < DoubleVector.SPECIES_PREFERRED.length(); i++) {
            for (int j = 0; j < 10; j++) {
                qI[j][i] = Qi[j];
                qK[j][i] = Qk[j];
            }
            for (int j = 0; j < 3; j++) {
                uI[j][i] = Ui[j];
                uK[j][i] = Uk[j];
            }
        }
        Fi = new DoubleVector[3];
        Fk = new DoubleVector[3];
        Ti = new DoubleVector[3];
        Tk = new DoubleVector[3];
        for(int i = 0; i < 3; i++){
            Fi[i] = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 0.0);
            Fk[i] = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 0.0);
            Ti[i] = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 0.0);
            Tk[i] = DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 0.0);
        }
    }

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Order 5 Coulomb", (Object) 5, Operator.COULOMB}});
    }

    @Test
    public void permanentMultipoleEnergyAndGradTest() {
        // Thole damping is not used for permanent AMOEBA electrostatics.
        if (operator == Operator.THOLE_FIELD) {
            return;
        }

        MultipoleTensorSIMD multipoleTensor = new CoulombTensorGlobalSIMD(order);
        PolarizableMultipoleSIMD mI = new PolarizableMultipoleSIMD(qI, uI, uI);
        PolarizableMultipoleSIMD mK = new PolarizableMultipoleSIMD(qK, uK, uK);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector e = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);
        double[] eArray = new double[DoubleVector.SPECIES_PREFERRED.length()];
        double[][] torI = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        double[][] torK = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        for(int i = 0; i < 3; i++){
            torI[i] = Ti[i].toArray();
            torK[i] = Tk[i].toArray();
        }
        e.intoArray(eArray, 0);
        for (int i = 0; i < DoubleVector.SPECIES_PREFERRED.length(); i++) {
            assertEquals(format("Permanent energy %s %d", info, i), permanentEnergy, eArray[i], tolerance);
            assertEquals(info + " Cartesian Multipole x-Torque I - Lane " + i, permTorqueI[0], torI[0][i], tolerance);
            assertEquals(info + " Cartesian Multipole x-Torque K - Lane " + i, permTorqueK[0], torK[0][i], tolerance);
            assertEquals(info + " Cartesian Multipole y-Torque I - Lane " + i, permTorqueI[1], torI[1][i], tolerance);
            assertEquals(info + " Cartesian Multipole y-Torque K - Lane " + i, permTorqueK[1], torK[1][i], tolerance);
            assertEquals(info + " Cartesian Multipole z-Torque I - Lane " + i, permTorqueI[2], torI[2][i], tolerance);
            assertEquals(info + " Cartesian Multipole z-Torque K - Lane " + i, permTorqueK[2], torK[2][i], tolerance);
        }

        // Analytic gradient.
        DoubleVector aX = Fk[0];
        DoubleVector aY = Fk[1];
        DoubleVector aZ = Fk[2];

        // Finite difference gradient.
        // X direction
        double delta = 1.0e-5;
        double delta2 = 2.0 * 1.0e-5;
        xVector = xVector.add(delta);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector posX = multipoleTensor.multipoleEnergy(mI, mK);
        xVector = xVector.sub(delta2);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector negX = multipoleTensor.multipoleEnergy(mI, mK);
        xVector = xVector.add(delta);

        // Y direction
        yVector = yVector.add(delta);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector posY = multipoleTensor.multipoleEnergy(mI, mK);
        yVector = yVector.sub(delta2);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector negY = multipoleTensor.multipoleEnergy(mI, mK);
        yVector = yVector.add(delta);

        // Z direction
        zVector = zVector.add(delta);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector posZ = multipoleTensor.multipoleEnergy(mI, mK);
        zVector = zVector.sub(delta2);
        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector negZ = multipoleTensor.multipoleEnergy(mI, mK);
        zVector = zVector.add(delta);

        double[] expectArray = aX.toArray();
        DoubleVector actual = posX.sub(negX).div(delta2);
        double[] actualArray = actual.toArray();
        double[] expectYArray = aY.toArray();
        DoubleVector actualY = posY.sub(negY).div(delta2);
        double[] actualYArray = actualY.toArray();
        double[] expectZArray = aZ.toArray();
        DoubleVector actualZ = posZ.sub(negZ).div(delta2);
        double[] actualZArray = actualZ.toArray();

        for(int i = 0; i < DoubleVector.SPECIES_PREFERRED.length(); i++){
            assertEquals(format("Permanent Grad %s - Lane %d", info, i), expectArray[i], actualArray[i], fdTolerance);
            assertEquals(format("Permanent Grad %s - Lane %d", info, i), expectYArray[i], actualYArray[i], fdTolerance);
            assertEquals(format("Permanent Grad %s - Lane %d", info, i), expectZArray[i], actualZArray[i], fdTolerance);
        }
    }

    @Test
    public void polarizationEnergyAndGradientTest(){
        MultipoleTensorSIMD multipoleTensor = new CoulombTensorGlobalSIMD(order);

        PolarizableMultipoleSIMD mI = new PolarizableMultipoleSIMD(qI, uI, uI);
        PolarizableMultipoleSIMD mK = new PolarizableMultipoleSIMD(qK, uK, uK);
        double[][] uIEwald = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        double[][] uKEwald = new double[3][DoubleVector.SPECIES_PREFERRED.length()];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < DoubleVector.SPECIES_PREFERRED.length(); j++) {
                uIEwald[i][j] = UiEwald[i];
                uKEwald[i][j] = UkEwald[i];
            }
        }
        if (operator == Operator.SCREENED_COULOMB) {
            mI.setInducedDipole(uIEwald, uIEwald);
            mK.setInducedDipole(uKEwald, uKEwald);
        }

        multipoleTensor.setR(xVector, yVector, zVector);
        multipoleTensor.generateTensor();
        DoubleVector e = multipoleTensor.polarizationEnergyAndGradient(mI, mK,
                DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 1.0),
                DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, 1.0),
                DoubleVector.broadcast(DoubleVector.SPECIES_PREFERRED, scaleMutual),
                Fi, Ti, Tk);
        double[] eArray = e.toArray();
        // Analytic gradient.
        double[] aX = Fi[0].toArray();
        double[] aY = Fi[1].toArray();
        double[] aZ = Fi[2].toArray();
        double[] axTorqueI = Ti[0].toArray();
        double[] axTorqueK = Tk[0].toArray();
        double[] ayTorqueI = Ti[1].toArray();
        double[] ayTorqueK = Tk[1].toArray();
        double[] azTorqueI = Ti[2].toArray();
        double[] azTorqueK = Tk[2].toArray();

        double energy = polarizationEnergyCoulomb;
        double[] gradI = polarGradICoulomb;
        double[] torqueI = polarTorqueICoulomb;
        double[] torqueK = polarTorqueKCoulomb;
        for(int i = 0; i < DoubleVector.SPECIES_PREFERRED.length(); i++) {
            assertEquals(info + " Polarization Energy", energy, eArray[i], tolerance);
            assertEquals(info + " Polarization GradX", gradI[0], aX[i], tolerance);
            assertEquals(info + " Polarization GradY", gradI[1], aY[i], tolerance);
            assertEquals(info + " Polarization GradZ", gradI[2], aZ[i], tolerance);
            assertEquals(info + " Polarization TorqueX I", torqueI[0], axTorqueI[i], tolerance);
            assertEquals(info + " Polarization TorqueY I", torqueI[1], ayTorqueI[i], tolerance);
            assertEquals(info + " Polarization TorqueZ I", torqueI[2], azTorqueI[i], tolerance);
            assertEquals(info + " Polarization TorqueX K", torqueK[0], axTorqueK[i], tolerance);
            assertEquals(info + " Polarization TorqueY K", torqueK[1], ayTorqueK[i], tolerance);
            assertEquals(info + " Polarization TorqueZ K", torqueK[2], azTorqueK[i], tolerance);
        }

        // Finite difference grad not tested due the ind. d. grad dependence on ind. d. values
    }
}
