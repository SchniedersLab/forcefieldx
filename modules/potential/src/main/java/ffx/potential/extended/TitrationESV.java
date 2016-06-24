package ffx.potential.extended;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.Supplier;
import java.util.stream.Collectors;

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
        for (MSNode node : residueOne.getTerms().getChildList()) {
            if (node instanceof ROLS) {
                rolsOne.add((ROLS) node);
            }
        }
        for (MSNode node : residueZero.getTerms().getChildList()) {
            if (node instanceof ROLS) {
                rolsZero.add((ROLS) node);
            }
        }
        setAtoms();
    }
    
    @Override
    protected void setAtoms() {
        atoms = new ArrayList<>();
        atoms.addAll(residueOne.getAtomList());
        atoms.addAll(residueZero.getAtomList());
        atoms.parallelStream().forEach(a -> a.setApplyLamedh(true));
    }
    
    @Override
    public double getdEdLamedh() {
        return ffe.getdEdLdh()[index];
    }
    
    @Override
    public OptionalDouble getROLSScaling(ROLS rols) {
        if (rolsOne.contains(rols)) {
            return OptionalDouble.of(lamedh);
        } else if (rolsZero.contains(rols)) {
            return OptionalDouble.of(1.0 - lamedh);
        }
        return OptionalDouble.empty();
    }
    
}
