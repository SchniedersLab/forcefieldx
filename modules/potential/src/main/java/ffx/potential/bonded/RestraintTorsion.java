package ffx.potential.bonded;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.parameters.TorsionType;

import java.util.function.DoubleUnaryOperator;

import static org.apache.commons.math3.util.FastMath.*;

public class RestraintTorsion extends BondedTerm implements LambdaInterface {

    private final Atom[] atoms;
    private final TorsionType torsionType;
    private final boolean lambdaTerm;
    private final DoubleUnaryOperator lamMapper;
    private final double units;

    private double lambda = 1.0;
    private double dEdL = 0.0;
    private double d2EdL2 = 0.0;

    public RestraintTorsion(Atom a1, Atom a2, Atom a3, Atom a4, TorsionType tType, boolean lamEnabled, boolean revLam, double units) {
        atoms = new Atom[]{a1, a2, a3, a4};
        this.torsionType = tType;
        lambdaTerm = lamEnabled;
        if (this.lambdaTerm) {
            if (revLam) {
                lamMapper = (double l) -> 1.0 - l;
            } else {
                lamMapper = (double l) -> l;
            }
        } else {
            lamMapper = (double l) -> 1.0;
        }
        this.units = units;
    }

    @Override
    public double energy(boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
        energy = 0.0;
        value = 0.0;
        dEdL = 0.0;
        var atomA = atoms[0];
        var atomB = atoms[1];
        var atomC = atoms[2];
        var atomD = atoms[3];
        var va = atomA.getXYZ();
        var vb = atomB.getXYZ();
        var vc = atomC.getXYZ();
        var vd = atomD.getXYZ();
        var vba = vb.sub(va);
        var vcb = vc.sub(vb);
        var vdc = vd.sub(vc);
        var vt = vba.X(vcb);
        var vu = vcb.X(vdc);
        var rt2 = vt.length2();
        var ru2 = vu.length2();
        var rtru2 = rt2 * ru2;
        if (rtru2 != 0.0) {
            var rr = sqrt(rtru2);
            var rcb = vcb.length();
            var cosine = vt.dot(vu) / rr;
            var sine = vcb.dot(vt.X(vu)) / (rcb * rr);
            value = toDegrees(acos(cosine));
            if (sine < 0.0) {
                value = -value;
            }
            var amp = torsionType.amplitude;
            var tsin = torsionType.sine;
            var tcos = torsionType.cosine;
            energy = amp[0] * (1.0 + cosine * tcos[0] + sine * tsin[0]);
            var dedphi = amp[0] * (cosine * tsin[0] - sine * tcos[0]);
            var cosprev = cosine;
            var sinprev = sine;
            var n = torsionType.terms;
            for (int i = 1; i < n; i++) {
                var cosn = cosine * cosprev - sine * sinprev;
                var sinn = sine * cosprev + cosine * sinprev;
                var phi = 1.0 + cosn * tcos[i] + sinn * tsin[i];
                var dphi = (1.0 + i) * (cosn * tsin[i] - sinn * tcos[i]);
                energy = energy + amp[i] * phi;
                dedphi = dedphi + amp[i] * dphi;
                cosprev = cosn;
                sinprev = sinn;
            }
            if (esvTerm) {
                esvDerivLocal = units * energy * dedesvChain * lambda;
            }
            energy = units * energy * esvLambda * lambda;
            dEdL = units * energy * esvLambda;
            if (gradient || lambdaTerm) {
                dedphi = units * dedphi * esvLambda;
                var vca = vc.sub(va);
                var vdb = vd.sub(vb);
                var dedt = vt.X(vcb).scaleI(dedphi / (rt2 * rcb));
                var dedu = vu.X(vcb).scaleI(-dedphi / (ru2 * rcb));
                var ga = dedt.X(vcb);
                var gb = vca.X(dedt).addI(dedu.X(vdc));
                var gc = dedt.X(vba).addI(vdb.X(dedu));
                var gd = dedu.X(vcb);
                int ia = atomA.getIndex() - 1;
                int ib = atomB.getIndex() - 1;
                int ic = atomC.getIndex() - 1;
                int id = atomD.getIndex() - 1;
                if (lambdaTerm) {
                    lambdaGrad.add(threadID, ia, ga);
                    lambdaGrad.add(threadID, ib, gb);
                    lambdaGrad.add(threadID, ic, gc);
                    lambdaGrad.add(threadID, id, gd);
                }
                if (gradient) {
                    grad.add(threadID, ia, ga.scaleI(lambda));
                    grad.add(threadID, ib, gb.scaleI(lambda));
                    grad.add(threadID, ic, gc.scaleI(lambda));
                    grad.add(threadID, id, gd.scaleI(lambda));
                }
            }
        }

        return energy;
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public void setLambda(double lambda) {
        this.lambda = lamMapper.applyAsDouble(lambda);
    }

    @Override
    public double getd2EdL2() {
        return d2EdL2;
    }

    @Override
    public double getdEdL() {
        return dEdL;
    }

    @Override
    public void getdEdXdL(double[] gradient) {
        // The chain rule term is at least supposedly zero.
    }

    @Override
    public String toString() {
        return String.format(" t-type %s, val %.3f, e %.3f", torsionType, value, energy);
    }
}
