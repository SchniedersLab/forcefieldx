package ffx.potential.extended;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedUtils;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.extended.TitrationESV.TitrationUtils.Titr;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;
import static java.lang.String.format;
import java.util.Arrays;

import ffx.numerics.Potential;

/**
 * An extended system variable that allows continuous fractional protonation of an amino acid.
 * All atomic charges and bonded terms scale linearly between prot and deprot states.
 * @author slucore
 */
public final class TitrationESV extends ExtendedVariable {
    
    // System handles
    private static final Logger logger = Logger.getLogger(TitrationESV.class.getName());
    private final MultiResidue titrating;
    
    // Handles on scaled terms
    private Residue residueOne, residueZero;    // One*lamedh + Zero*(1-lamedh)
    private List<Atom> atomsOne, atomsZero;     // just those that are changing with lamedh
    private List<ROLS> rolsOne, rolsZero;
    private static final List<String> backboneNames = Arrays.asList("N","CA","C","O","HA","H");
    
    // Options
    private static final boolean HIEmode = false;
    
    public TitrationESV(MultiResidue titrating, 
            double temperature, double dt, double biasMag) {
        super(biasMag, 1.0);
        this.titrating = titrating;
        
        residueOne = titrating.getActive();
        residueZero = titrating.getInactive().get(0);
        rolsOne = new ArrayList<>();
        residueOne.getTerms().getChildList()
                .parallelStream().filter(node -> node instanceof ROLS)
                .forEachOrdered(node -> rolsOne.add((ROLS) node));
        rolsZero = new ArrayList<>();
        residueZero.getTerms().getChildList()
                .parallelStream().filter(node -> node instanceof ROLS)
                .forEachOrdered(node -> rolsZero.add((ROLS) node));
    }
    
    public TitrationESV(MultiResidue titrating, double temperature, double dt) {
        this(titrating, temperature, dt, 1.0);
    }
    
    @Override
    protected void finalize() {
        atomsOne = new ArrayList<>();
        atomsZero = new ArrayList<>();
        List<Atom> backboneAtoms = new ArrayList<>();
        backboneNames.parallelStream().forEach(name -> {
            backboneAtoms.add((Atom) residueOne.getAtomNode(name));
        });
        residueOne.getAtomList().parallelStream().filter(a -> !backboneAtoms.contains(a)).forEachOrdered(atomsOne::add);
        residueZero.getAtomList().parallelStream().filter(a -> !backboneAtoms.contains(a)).forEachOrdered(atomsZero::add);
        atoms.addAll(atomsOne);
        atoms.addAll(atomsZero);
        atoms.parallelStream().forEach(a -> a.setESV(this));
        describe();
    }
    
//    /**
//     * Query FFE for derivative terms, add bias.
//     */
//    @Override
//    public double getdEdLdh() {
//        return (ffe.getdEdLdh()[index] + getdBiasdLdh());
//    }
//    
//    /**
//     * Query FFE for derivative terms, add bias.
//     */
//    @Override
//    public double getd2EdLdh2() {
//        return (ffe.getd2EdLdh2()[index] + getd2BiasdLdh2());
//    }
    
    @Override
    public OptionalDouble getROLSScaling(ROLS rols) {
        for (ROLS termNode : rolsOne) {
            if (((MSNode) termNode).getChildList().contains(rols)) {
                return OptionalDouble.of(lambda);
            }
        }
        for (ROLS termNode : rolsZero) {
            if (((MSNode) termNode).getChildList().contains(rols)) {
                return OptionalDouble.of(1.0 - lambda);
            }
        }
        return OptionalDouble.empty();
    }
    
    /**
     * Eqs 5,6 from Wallace+Shen 2011 "Continuous constant pH M.D. in explicit..."
     * U_pH(ldh) = log(10)*kb*T*(pKa_model - pH)*ldh
     * U_mod(ldh) = potential of mean force for protonation (or -deprot) of model compound
     * U_star = sum(ldh) { U_pH(ldh) + U_mod_prot(ldh) + U_barr(ldh)
     * This method returns U_pH + U_mod_prot.
     */
    public double getPhEnergy(double pH, double temperature) {
        Titr titration = TitrationUtils.titrationLookup(this.titrating.getActive());
        double pKaModel = titration.pKa;
        double uph = Math.log(10)*ThermoConstants.kB*temperature*(pKaModel - pH)*lambda;
        double umod = 0.0;  // TODO PRIO find PMFs for monomers/trimers/pentapeptides
        umod = titration.refEnergy * lambda;
        return uph + umod;
    }
    
    @Override
    public String toString() {
        return format("Titration ESV %d: %s", index, titrating);
    }
    
    /**
     * List all the atoms and bonded terms associated with each end state.
     */
    public void describe() {
        StringBuilder sb = new StringBuilder();
        sb.append(this.toString() + format("\n"));
        sb.append(format("    Atoms (State 1): \n"));
        atomsOne.stream().forEachOrdered(a -> sb.append(format("   %s\n", a)));
        sb.append(format("    ROLS (State 1):  \n"));
        rolsOne.stream().forEachOrdered(r -> sb.append( format("      %s\n", r)));
        sb.append(format("    Atoms (State 0): \n"));
        atomsZero.stream().forEachOrdered(a -> sb.append(format("   %s\n", a)));
        sb.append(format("    ROLS (State 0):  \n"));
        rolsZero.stream().forEachOrdered(r -> sb.append( format("      %s\n", r)));
        logger.info(sb.toString());
    }
    
    /**
     * Helper methods to define titration-specific phenomena.
     */
    public static class TitrationUtils {

        public static MultiResidue titrationFactory(MolecularAssembly mola, Residue res) {
            ForceField ff = mola.getForceField();
            Potential potential = mola.getPotentialEnergy();
            if (!(potential instanceof ForceFieldEnergy)) {
                logger.severe("TitrationFactory only supported by ForceFieldEnergy potentials.");
            }
            ForceFieldEnergy ffe = (ForceFieldEnergy) potential;
            // Create new titration state.
            Titr t = TitrationUtils.titrationLookup(res);
            String targetName = (t.target != res.getAminoAcid3()) ? t.target.toString() : t.source.toString();
            int resNumber = res.getResidueNumber();
            Residue.ResidueType resType = res.getResidueType();
            Residue newRes = new Residue(targetName, resNumber, resType);
            // Wrap both states in a MultiResidue.
            MultiResidue multiRes = new MultiResidue(res, ff, ffe);
            Polymer polymer = findResiduePolymer(res, mola);
            polymer.addMultiResidue(multiRes);
            multiRes.addResidue(newRes);
            multiRes.setActiveResidue(res);
            propagateInactiveResidues(multiRes);
            ffe.reInit();
            logger.info(String.format(" Added Ldh-coupled titrating group: %s", res));
            return multiRes;
        }
        
        public static Titr titrationLookup(Residue res) {
            AminoAcid3 source = AminoAcid3.valueOf(res.getName());
            if (source == AminoAcid3.HIS || source == AminoAcid3.HID || source == AminoAcid3.HIE) {
                return (HIEmode ? Titr.ZtoH : Titr.UtoH);
            }
            for (Titr titr : Titr.values()) {
                if (titr.source == source || titr.target == source) {
                    return titr;
                }
            }
            logger.log(Level.SEVERE, "No titration lookup found for residue {0}", res);
            return null;
        }

        private static Polymer findResiduePolymer(Residue residue, MolecularAssembly mola) {
            Polymer polymers[] = mola.getChains();
            Optional<Polymer> polymer = IntStream.range(0,polymers.length).parallel()
                    .mapToObj(i -> polymers[i])
                    .filter(p -> p.getChainID().compareTo(residue.getChainID()) == 0)
                    .findAny();
            if (!polymer.isPresent()) {
                logger.log(Level.SEVERE, " Polymer not found for residue {0}", residue);
            }
            return polymer.get();
        }
        
        /**
         * Copies atomic coordinates from each active residue to its inactive
         * counterparts. Inactive hydrogen coordinates are updated by geometry
         * with the propagated heavies.
         */
        private static void propagateInactiveResidues(MultiResidue multiRes) {
            // Propagate all atom coordinates from active residues to their inactive counterparts.
            Residue active = multiRes.getActive();
            String activeResName = active.getName();
            List<Residue> inactives = multiRes.getInactive();
            for (Atom activeAtom : active.getAtomList()) {
                String activeName = activeAtom.getName();
                for (Residue inactive : inactives) {
                    Atom inactiveAtom = (Atom) inactive.getAtomNode(activeName);
                    if (inactiveAtom != null) {
                        // Propagate position and gradient.
                        double activeXYZ[] = activeAtom.getXYZ(null);
                        inactiveAtom.setXYZ(activeXYZ);
                        double grad[] = new double[3];
                        activeAtom.getXYZGradient(grad);
                        inactiveAtom.setXYZGradient(grad[0], grad[1], grad[2]);
                        // Propagate velocity, acceleration, and previous acceleration.
                        double activeVelocity[] = new double[3];
                        activeAtom.getVelocity(activeVelocity);
                        inactiveAtom.setVelocity(activeVelocity);
                        double activeAccel[] = new double[3];
                        activeAtom.getAcceleration(activeAccel);
                        inactiveAtom.setAcceleration(activeAccel);
                        double activePrevAcc[] = new double[3];
                        activeAtom.getPreviousAcceleration(activePrevAcc);
                        inactiveAtom.setPreviousAcceleration(activePrevAcc);
                    } else {
                        if (activeName.equals("C") || activeName.equals("O") || activeName.equals("N") || activeName.equals("CA")
                                || activeName.equals("H") || activeName.equals("HA")) {
                            // Backbone atoms aren't supposed to exist in inactive multiResidue components; so no problem.
                        } else if ((activeResName.equals("LYS") && activeName.equals("HZ3"))
                                || (activeResName.equals("TYR") && activeName.equals("HH"))
                                || (activeResName.equals("CYS") && activeName.equals("HG"))
                                || (activeResName.equals("HIS") && (activeName.equals("HD1") || activeName.equals("HE2")))
                                || (activeResName.equals("HID") && activeName.equals("HD1"))
                                || (activeResName.equals("HIE") && activeName.equals("HE2"))
                                || (activeResName.equals("ASH") && activeName.equals("HD2"))
                                || (activeResName.equals("GLH") && activeName.equals("HE2"))) {
                            // These titratable protons are handled below; so no problem.
                        } else {
                            // Now we have a problem.
                            logger.warning(String.format("Couldn't copy atom_xyz: %s: %s, %s",
                                    multiRes, activeName, activeAtom.toString()));
                        }
                    }
                }
            }

            // If inactive residue is a protonated form, move the stranded hydrogen to new coords (based on propagated heavies).
            // Also give the stranded hydrogen a maxwell velocity and remove its accelerations.
            for (Residue inactive : inactives) {
                List<Atom> resetMe = new ArrayList<>();
                switch (inactive.getName()) {
                    case "LYS": {
                        Atom HZ3 = (Atom) inactive.getAtomNode("HZ3");
                        Atom NZ = (Atom) inactive.getAtomNode("NZ");
                        Atom CE = (Atom) inactive.getAtomNode("CE");
                        Atom HZ1 = (Atom) inactive.getAtomNode("HZ1");
                        BondedUtils.intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
                        resetMe.add(HZ3);
                        break;
                    }
                    case "ASH": {
                        Atom HD2 = (Atom) inactive.getAtomNode("HD2");
                        Atom OD2 = (Atom) inactive.getAtomNode("OD2");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom OD1 = (Atom) inactive.getAtomNode("OD1");
                        BondedUtils.intxyz(HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0);
                        resetMe.add(HD2);
                        break;
                    }
                    case "GLH": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom OE2 = (Atom) inactive.getAtomNode("OE2");
                        Atom CD = (Atom) inactive.getAtomNode("CD");
                        Atom OE1 = (Atom) inactive.getAtomNode("OE1");
                        BondedUtils.intxyz(HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0);
                        resetMe.add(HE2);
                        break;
                    }
                    case "HIS": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom NE2 = (Atom) inactive.getAtomNode("NE2");
                        Atom CD2 = (Atom) inactive.getAtomNode("CD2");
                        Atom CE1 = (Atom) inactive.getAtomNode("CE1");
                        Atom HD1 = (Atom) inactive.getAtomNode("HD1");
                        Atom ND1 = (Atom) inactive.getAtomNode("ND1");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
                        BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                        resetMe.add(HE2);
                        resetMe.add(HD1);
                        break;
                    }
                    case "HID": {
                        Atom HD1 = (Atom) inactive.getAtomNode("HD1");
                        Atom ND1 = (Atom) inactive.getAtomNode("ND1");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                        resetMe.add(HD1);
                        break;
                    }
                    case "HIE": {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom NE2 = (Atom) inactive.getAtomNode("NE2");
                        Atom CD2 = (Atom) inactive.getAtomNode("CD2");
                        Atom CE1 = (Atom) inactive.getAtomNode("CE1");
                        BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
                        resetMe.add(HE2);
                        break;
                    }
                    case "CYS": {
                        Atom HG = (Atom) inactive.getAtomNode("HG");
                        Atom SG = (Atom) inactive.getAtomNode("SG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        Atom CA = (Atom) inactive.getAtomNode("CA");
                        BondedUtils.intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
                        resetMe.add(HG);
                        break;
                    }
                    case "TYR": {
                        Atom HH = (Atom) inactive.getAtomNode("HH");
                        Atom OH = (Atom) inactive.getAtomNode("OH");
                        Atom CZ = (Atom) inactive.getAtomNode("CZ");
                        Atom CE2 = (Atom) inactive.getAtomNode("CE2");
                        BondedUtils.intxyz(HH, OH, 0.97, CZ, 108.0, CE2, 0.0, 0);
                        resetMe.add(HH);
                        break;
                    }
                    default:
                }
                for (Atom a : resetMe) {
                    a.setXYZGradient(0, 0, 0);
                    a.setVelocity(ThermoConstants.singleRoomtempMaxwell(a.getMass()));
                    a.setAcceleration(new double[]{0, 0, 0});
                    a.setPreviousAcceleration(new double[]{0, 0, 0});
                }
            }
        }
    
        /**
         * All described as protonation reactions.
         */
        public enum Titr {
            ctoC( 8.18, +60.168, AminoAcid3.CYD, AminoAcid3.CYS),
            Dtod( 3.90, +53.188, AminoAcid3.ASP, AminoAcid3.ASH),
            Etoe( 4.25, +59.390, AminoAcid3.GLU, AminoAcid3.GLH),
            ktoK(10.53, -50.440, AminoAcid3.LYD, AminoAcid3.LYS),
            ytoY(10.07, +34.961, AminoAcid3.TYD, AminoAcid3.TYR),
            UtoH( 6.00, -42.923, AminoAcid3.HID, AminoAcid3.HIS),
            ZtoH( 6.00, +00.000, AminoAcid3.HIE, AminoAcid3.HIS),
            TerminusNH3toNH2 (8.23, +00.00, AminoAcid3.UNK, AminoAcid3.UNK),
            TerminusCOOHtoCOO(3.55, +00.00, AminoAcid3.UNK, AminoAcid3.UNK);

            public final double pKa, refEnergy;
            public final AminoAcid3 source, target;

            Titr(double pKa, double refEnergy, AminoAcid3 source, AminoAcid3 target) {
                this.pKa = pKa;
                this.refEnergy = refEnergy;
                this.source = source;
                this.target = target;
            }
        }
    }
    
}
