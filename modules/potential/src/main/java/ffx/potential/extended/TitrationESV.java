package ffx.potential.extended;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import java.util.List;
import java.util.Optional;

/**
 * An extended system variable that allows continuous fractional protonation of an amino acid.
 * All atomic charges and bonded terms scale linearly between prot and deprot states.
 * @author slucore
 */
public final class TitrationESV extends ExtendedVariable {
    
    // System handles
    private final ForceFieldEnergy ffe;
    // Handles on scaled terms
    private Residue residueOne, residueZero;    // One*lamedh + Zero*(1-lamedh)
    private List<Atom> atomsOne, atomsZero;     // just those that are changing with lamedh
    private List<ROLS> rolsOne, rolsZero;
    // Derivatives
    private double dEdLamedh = 0.0;
    
    public TitrationESV(MultiResidue titrating, MolecularAssembly mola, double temperature, double dt) {
        super(mola, temperature, dt);
        Potential potential = mola.getPotentialEnergy();
        ffe = (potential instanceof ForceFieldEnergy) ? (ForceFieldEnergy) potential : null;
        if (ffe == null) {
            throw new UnsupportedOperationException("Lamedh ESVs not yet implemented for non-forcefield or hybrid potentials.");
        }
        setLamedh(1.0);
        List<Residue> members = titrating.getConsideredResidues();
        if (members.size() > 2) {
            throw new UnsupportedOperationException("TitrationESV not yet implemented for MultiResidues with >2 states.");
        }
        residueOne = titrating.getActive();
        residueZero = members.get(0).equals(residueOne) ? members.get(1) : members.get(0);
        for (MSNode node : residueOne.getChildList()) {
            if (node instanceof ROLS) {
                rolsOne.add((ROLS) node);
            }
        }
        for (MSNode node : residueZero.getChildList()) {
            if (node instanceof ROLS) {
                rolsZero.add((ROLS) node);
            }
        }
    }
    
    @Override
    public double getdEdLamedh() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public Optional<Double> getROLSScaling(ROLS rols) {
        if (rolsOne.contains(rols)) {
            return Optional.of(lamedh);
        } else if (rolsZero.contains(rols)) {
            return Optional.of(1.0 - lamedh);
        }
        return null;
    }
    
    public void reference_DualTopology_dEdL() {
//        double e1 = lambdaPow * dEdL_1 + dLambdaPow * energy1
//                + oneMinusLambdaPow * restraintdEdL_1 + dOneMinusLambdaPow * restraintEnergy1;
//        double e2 = oneMinusLambdaPow * dEdL_2 + dOneMinusLambdaPow * energy2
//                + lambdaPow * restraintdEdL_2 + dLambdaPow * restraintEnergy2;
//        double e1 = (L * dE) + (dL * E) + ((1-L) * d(Eres)) + (d(1-L) * Eres);    where E <- E1
//        double e2 = ((1-L) * dE) + (d(1-L) * E) + (L * d(Eres)) + (dL * Eres);    where E <- E2

    // When lambdaExponent = 1.0:
//        double e1 = (L * dE) + E + ((1-L) * d(Eres)) + -(Eres);    where E <- E1
//        double e2 = ((1-L) * dE) + -(E) + (L * d(Eres)) + Eres;    where E <- E2
//        return e1 + e2;
    }
    
    public void reference_DualTopology_setLambda() {
//        this.lamedh = lamedh;
//        oneMinusLambda = 1.0 - lamedh;
//
//        lambdaPow = pow(lamedh, lambdaExponent);
//        dLambdaPow = lambdaExponent * pow(lamedh, lambdaExponent - 1.0);
//
//        oneMinusLambdaPow = pow(oneMinusLambda, lambdaExponent);
//        dOneMinusLambdaPow = -lambdaExponent * pow(oneMinusLambda, lambdaExponent - 1.0);

    // When lambdaExponent = 1.0:
//        this.lamedh = lamedh;
//        oneMinusLambda = 1.0 - lamedh;
//        lambdaPow = lamedh;
//        oneMinusLambdaPow = oneMinusLambda;
//        dLambdaPow = 1.0;
//        dOneMinusLambdaPow = -1.0;
    }
    
}
