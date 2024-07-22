package ffx.numerics.multipole;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import static ffx.numerics.multipole.MultipoleTensorTest.*;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import static org.junit.Assert.assertEquals;

@RunWith(Parameterized.class)
public class CombinedTensorTest {
    /**
     * Logger for the MultipoleTensor class.
     */
    private static final Logger logger = Logger.getLogger(AMOEBAPlusMultipoleTensorGlobalTest.class.getName());
    private final double tolerance = 1.0e-13;
    private final double fdTolerance = 1.0e-6;
    private final double[] r2 = new double[3];
    private final int order;
    private final String info;
    private final Operator operator;
    private final double ewaldA;
    private final double thole;
    private final double aiak;
    private final double alphaI;
    private final double alphaK;
    private double[] ewaldTerms;
    private double[] dampTerms;
    private double[] overlapTerms;
    private double[] tholeDirect;
    private double[] tholeTerms;
    private double R;

    public CombinedTensorTest(
            String info,
            int order,
            Operator operator) {
        this.info = info;
        this.order = order;
        this.operator = operator;
        this.ewaldA = 0.54459051633620059;
        this.thole = .7;
        this.aiak = 1.0 / 0.86460072979796410;
        this.alphaI = MultipoleTensorTest.chargePenAlphaOx;
        this.alphaK = MultipoleTensorTest.chargePenAlphaHyd;
        updateSources(MultipoleTensorTest.xyzEwaldAPlus);
    }

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(
                new Object[][]{
                        {"Order 5 Aplus Direct", 5, Operator.AMOEBA_PLUS_OVERLAP_FIELD}
                });
    }

    private void updateSources(double[] xyz){
        this.R = sqrt(xyz[0]*xyz[0]
                + xyz[1]*xyz[1]
                + xyz[2]*xyz[2]);
        // Use static tensor methods for terms fill with ones to not include coulomb terms already
        this.ewaldTerms = new double[order+1];
        Arrays.fill(ewaldTerms, 1.0);
        EwaldTensorGlobal.fillEwaldSource(order, ewaldA,
                EwaldTensorGlobal.initEwaldSource(order, ewaldA, new double[order+1]),
                R, ewaldTerms);
        this.dampTerms = new double[order+1];
        Arrays.fill(dampTerms, 1.0);
        AmoebaPlusDampTensorGlobal.dampSource(alphaI, R, dampTerms);
        this.overlapTerms = new double[order+1];
        Arrays.fill(overlapTerms, 1.0);
        AmoebaPlusOverlapTensorGlobal.overlapSource(alphaI, alphaK,
                alphaK*alphaK/(alphaK*alphaK - alphaI*alphaI),
                alphaI*alphaI/(alphaI*alphaI - alphaK*alphaK), R,
                overlapTerms);
        this.tholeDirect = new double[order+1];
        Arrays.fill(tholeDirect, 1.0);
        TholeTensorGlobal.tholeSource(thole, aiak, R, true, tholeDirect);
        this.tholeTerms = new double[order+1];
        Arrays.fill(tholeTerms, 1.0);
        TholeTensorGlobal.tholeSource(thole, aiak, R, false, tholeTerms);
    }

    private void buildOverlapSource(CombinedTensorGlobal multipoleTensor){
        multipoleTensor.resetSource();
        multipoleTensor.addTermsSeparate(ewaldTerms);
        multipoleTensor.addTerms(overlapTerms);
        multipoleTensor.addCoulombMultiplier(-1.0);
    }

    @Test
    public void permanentMultipoleEnergyAndGradTest() {
        double delta = 1.0e-8;
        double delta2 = 2.0 * delta;
        double[] Fi = new double[3];
        double[] Fk = new double[3];
        double[] Ti = new double[3];
        double[] Tk = new double[3];

        // Oxygen-Hydrogen A+ Real-Space Energy under Ewald
        double[] uI = new double[]{-9.5291976837210871E-002,  -5.1217206248926346E-002,   9.0954068003896951E-002};
        double[] uK = new double[]{-1.4736767767345298E-002,   1.3165041117687385E-002,  -7.4174987487104807E-005};
        PolarizableMultipole mI = new PolarizableMultipole(QlAPlusEwal,uI,uI,8.0);
        PolarizableMultipole mK = new PolarizableMultipole(QmAplusEwal,uK,uK,1.0);
        CombinedTensorGlobal multipoleTensor = new CombinedTensorGlobal(order);
        multipoleTensor.setR(xyzEwaldAPlus);

        //Amoeba+ = M-pol*(Ewald - Coulomb -CDamp)*M-pol
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor();
        double eOverlap = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);
        double overlapAns = 1.4311540110919852E-003;
        assertEquals("Real-Space PME Mpole-Mpole interaction", overlapAns, eOverlap, tolerance);

        // Analytic gradient.
        double aX = Fk[0];
        double aY = Fk[1];
        double aZ = Fk[2];

        xyzEwaldAPlus[0] += delta;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double posX = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[0] -= delta2;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double negX = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[0] += delta;
        xyzEwaldAPlus[1] += delta;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double posY = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[1] -= delta2;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double negY = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[1] += delta;
        xyzEwaldAPlus[2] += delta;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double posZ = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[2] -= delta2;
        updateSources(xyzEwaldAPlus);
        buildOverlapSource(multipoleTensor);
        multipoleTensor.generateTensor(xyzEwaldAPlus);
        double negZ = multipoleTensor.multipoleEnergy(mI, mK);
        xyzEwaldAPlus[2] += delta;

        double expect = aX;
        double actual = (posX - negX) / delta2;
        assertEquals(info + " OH Ewald Overlap Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OH Ewald Overlap Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OH Ewald Overlap Force Z", expect, actual, fdTolerance);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);

        // A+ CoreE = Core*(Ewald - Coulomb + Damp)*Mpole -> 2x for i->j j->i energies
        AmoebaPlusDampTensorGlobal newTensor = new AmoebaPlusDampTensorGlobal(3, alphaI, alphaK, ewaldA);
        newTensor.setR(xyzEwaldAPlus);
        newTensor.generateTensor();
        double eCore = newTensor.coreInteractionAndGradient(mI, mK, Fi, Fk);
        double answerCore = 1.8003552196793374E-006;
        assertEquals("Real-Space Core-Mpole Interaction", answerCore, eCore+overlapAns, tolerance);

        aX = Fk[0];
        aY = Fk[1];
        aZ = Fk[2];

        xyzEwaldAPlus[0] += delta;
        newTensor.generateTensor(xyzEwaldAPlus);
        posX = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[0] -= delta2;
        newTensor.generateTensor(xyzEwaldAPlus);
        negX = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[0] += delta;
        xyzEwaldAPlus[1] += delta;
        newTensor.generateTensor(xyzEwaldAPlus);
        posY = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[1] -= delta2;
        newTensor.generateTensor(xyzEwaldAPlus);
        negY = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[1] += delta;
        xyzEwaldAPlus[2] += delta;
        newTensor.generateTensor(xyzEwaldAPlus);
        posZ = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[2] -= delta2;
        newTensor.generateTensor(xyzEwaldAPlus);
        negZ = newTensor.coreInteraction(mI, mK);
        xyzEwaldAPlus[2] += delta;

        expect = aX;
        actual = (posX - negX) / delta2;
        assertEquals(info + " OH Ewald Core Force X", expect, actual, fdTolerance);
        expect = aY;
        actual = (posY - negY) / delta2;
        assertEquals(info + " OH Ewald Core Force Y", expect, actual, fdTolerance);
        expect = aZ;
        actual = (posZ - negZ) / delta2;
        assertEquals(info + " OH Ewald Core Force Z", expect, actual, fdTolerance);
        Arrays.fill(Fi, 0);
        Arrays.fill(Fk, 0);
    }

    @Test
    public void polarizationEnergyAndGradientTest() {
        // Oxygen-Hydrogen A+ Real-Space Direct Pol Energy under Ewald (aplussmall.xyz 1, 536)
        double[] Fi = new double[3];
        double[] Ti = new double[3];
        double[] Tk = new double[3];
        double[] uI = new double[]{-9.5291976837210871E-002,  -5.1217206248926346E-002,   9.0954068003896951E-002};
        double[] uK = new double[]{-1.4736767767345298E-002,   1.3165041117687385E-002,  -7.4174987487104807E-005};
        PolarizableMultipole mI = new PolarizableMultipole(QlAPlusEwal,uI,uI,8.0);
        PolarizableMultipole mK = new PolarizableMultipole(QmAplusEwal,uK,uK,1.0);
        CombinedTensorGlobal multipoleTensor = new CombinedTensorGlobal(order);
        multipoleTensor.setR(xyzEwaldAPlus);

        // A+ PolE = Dip*(Ewald-TholeDirect)*Mpole
        multipoleTensor = new CombinedTensorGlobal(order);
        multipoleTensor.addTermsSeparate(ewaldTerms);
        multipoleTensor.addTerms(tholeDirect);
        multipoleTensor.addCoulombMultiplier(-1.0);
        multipoleTensor.setR(xyzEwaldAPlus);
        multipoleTensor.generateTensor();
        double ePolOverlap = multipoleTensor.polarizationEnergyAndGradient(mI, mK, 1.0, 1.0, 0.0, Fi, Ti, Tk);
        double polAnswer = 4.8757773606071702E-006/2;
        assertEquals(" OH Real-Space PME PolE interaction", polAnswer, ePolOverlap, tolerance);

        double[] aplusIndGrad = new double[]{-1.1162913643840848E-005,
                6.4235992888231040E-006,   -6.2954494656288314E-006};
        assertEquals(info + " Polarization GradX", aplusIndGrad[0], Fi[0], tolerance);
        assertEquals(info + " Polarization GradY", aplusIndGrad[1], Fi[1], tolerance);
        assertEquals(info + " Polarization GradZ", aplusIndGrad[2], Fi[2], tolerance);
    }
}