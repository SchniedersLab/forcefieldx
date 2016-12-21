package ffx.potential.extended;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtUtils.SB;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.TitrationESV.TitrationUtils.isTitratableHydrogen;

/**
 * A generalized extended system variable.
 * Treatment of ESVs: 
 *  a. Bonded terms interpolate linearly between end states.
 *      ESVs based on MultiResidue (e.g. TitrationESV) place a multiplier in the term objects themselves.
 *  b. PME and vdW scaling and derivatives are handled inside these classes' inner loops.
 *      Softcoring follows the same form used by OSRW lambda.
 * @author slucore
 */
public abstract class ExtendedVariable {
    
    // System handles
    private static final Logger logger = Logger.getLogger(ExtendedVariable.class.getName());
    private static int esvIndexer = 0;
    public final int index;
    private boolean ready = false;
    
    // Properties
    protected static final boolean esvPropagation = prop("esv-propagation", true);
    protected static final boolean scaleBondedTerms = prop("esv-scaleBonded", true);
    protected static final OptionalDouble biasOverride = prop("esv-bias", OptionalDouble.empty());
    private final double thetaMass = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
    private final double thetaFriction = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps
    private final ExtInclusionMode extInclusion = 
            prop(ExtInclusionMode.class, "esv-atomInclusion", ExtInclusionMode.TITRATABLE_H);
    public enum ExtInclusionMode {
        ALL_ATOMS, TITRATABLE_H;
    }
    
    // Lambda and derivative variables
    private double lambda;                        // ESVs travel on {0,1}
    private double theta;                           // Propagates lamedh particle via "lamedh=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final Random stochasticRandom = ThreadLocalRandom.current();
    private final double betat;
    private double discrBias, dDiscrBiasdL;
    
    // Atom lists and scaled terms
    private MultiResidue multiRes;
    private List<Atom> atoms = new ArrayList<>();
    private List<Atom> backbone = new ArrayList<>();
    private Residue resOne, resZro;                 // resOne*lamedh + resZero*(1-lamedh)
    private List<Atom> atomsOne, atomsZro;          // side-chain only; atomsZro stay unconnected to mola
    private List<BondedTerm> bondedOne, bondedZro;  // valence terms for each side; mola won't see zro by default
    private MSNode termNode;                  // modified to contain all applicable bonded terms
    private int moleculeNumber = 0;
    
    public ExtendedVariable(MultiResidue multiRes, double biasMag, double initialLamedh) {
        index = esvIndexer++;
        betat = biasOverride.isPresent() ? biasOverride.getAsDouble() : biasMag;
        lambda = initialLamedh;
        theta = Math.asin(Math.sqrt(lambda));

        this.multiRes = multiRes;
        resOne = multiRes.getActive();
        termNode = resOne.getTerms();
        resZro = multiRes.getInactive().get(0);
        moleculeNumber = resOne.getAtomList().get(0).getMoleculeNumber();
    }
    
    public ExtendedVariable(MultiResidue multiRes, double biasMag) {
        this(multiRes, biasMag, 1.0);
    }    
    public ExtendedVariable(MultiResidue multiRes) {
        this(multiRes, 0.0, 1.0);
    }
    
    public int getMoleculeNumber() {
        return moleculeNumber;
    }
        
    /**
     * Propagate lambda using Langevin dynamics.
     * Check that temperature goes to the value used below (when set as a constant) even when sim is decoupled.
     */
    public void propagate(double dEdEsv, double dt, double setTemperature) {
        if (!esvPropagation) {
            return;
        }
        double rt2 = 2.0 * ExtConstants.Boltzmann * setTemperature * thetaFriction / dt;
        double randomForce = sqrt(rt2) * stochasticRandom.nextGaussian() / ExtConstants.forceToKcal;
        double dEdL = -dEdEsv * sin(2.0 * theta);
        halfThetaVelocity = (halfThetaVelocity * (2.0 * thetaMass - thetaFriction * dt)
                + ExtConstants.forceToKcalSquared * 2.0 * dt * (dEdL + randomForce))
                / (2.0 * thetaMass + thetaFriction * dt);
        theta = theta + dt * halfThetaVelocity;

        if (theta > PI) {
            theta -= 2.0 * PI;
        } else if (theta <= -PI) {
            theta += 2.0 * PI;
        }

        double sinTheta = sin(theta);
        setLambda(sinTheta * sinTheta);
    }
    
    public void setLambda(double lambda) {
        this.lambda = lambda;
        theta = Math.asin(Math.sqrt(lambda));
        discrBias = betat - 4*betat*(lambda-0.5)*(lambda-0.5);
        dDiscrBiasdL = -8*betat*(lambda-0.5);
        updateBondedLambdas();
    }
    
    public final double getLambda() {
        return lambda;
    }
    
    public final int getIndex() {
        return index;
    }
    
    @Override
    public String toString() {
        return String.format("ESV%d", index);
    }
    
    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double totalBiasEnergy(double temperature);
    public abstract double totalBiasDeriv(double temperature);
    
    /**
     * From Shen&Huang 2016; drives ESVs to zero/unity.
     * bias = 4B*(L-0.5)^2
     */
    public double discretizationBiasEnergy() {
        return discrBias;
    }
    
    /**
     * dBiasdL = -8B*(L-0.5)
     */
    public double discretizationBiasDeriv() {
        return dDiscrBiasdL;
    }
    
    /**
     * Fill the atoms array; set esvLambda in all bonded terms.
     */
    public void readyup() {
        atomsOne = new ArrayList<>();
        atomsZro = new ArrayList<>();
        // Fill the atom Lists
        for (String bbName : ExtConstants.backboneNames) {
            Atom bb = (Atom) resOne.getAtomNode(bbName);
            if (bb != null) {
                backbone.add(bb);
                bb.applyPersistentIndex();
            }
        }
        // Fill the atom lists.
        for (Atom a1 : resOne.getAtomList()) {
            if (!backbone.contains(a1)) {
                a1.setESV(this);
                a1.applyPersistentIndex();
                atomsOne.add(a1);
                a1.setEsvState(1);
            }
        }
        for (Atom a0 : resZro.getAtomList()) {
            if (!backbone.contains(a0) && !atomsOne.contains(a0)) {
                a0.setESV(this);
                a0.applyPersistentIndex();
                if ((extInclusion == ExtInclusionMode.TITRATABLE_H && isTitratableHydrogen(a0))
                        || extInclusion == ExtInclusionMode.ALL_ATOMS) {
                    atomsZro.add(a0);
                    a0.setEsvState(0);
                }                
            }
        }
        
        atoms.addAll(atomsOne);
        atoms.addAll(atomsZro);
        
        // Fill bonded term list and set all esvLambda values.
        bondedOne = resOne.getDescendants(BondedTerm.class);
//        logfn("getDescendants() wayA,wayB: %d %d", resOne.getDescendants(BondedTerm.class).size(), bondedOne.size());
        List<BondedTerm> remove = new ArrayList<>();
        for (BondedTerm bt1 : bondedOne) {
            List<Atom> termAtoms = Arrays.asList(bt1.getAtomArray());
            boolean found = false;
            for (Atom atom : atomsOne) {    // Remove BondedTerm object from the active residue if they never impact an ESV atom.
                if (termAtoms.contains(atom)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                remove.add(bt1);
            }
        }
        bondedOne.removeAll(remove);
        
        bondedZro = new ArrayList<>();
        for (BondedTerm term : resZro.getDescendants(BondedTerm.class)) {
            if (!bondedOne.contains(term)) {
                List<Atom> termAtoms = Arrays.asList(term.getAtomArray());
                boolean found = false;
                for (Atom atom : atomsZro) {    // Only include bonded terms if they impact an ESV'd atom.
                    if (termAtoms.contains(atom)) {
                        found = true;
                        break;
                    }
                }
                if (true) {
                    bondedZro.add(term);
                }
            }
        }
        MSNode bondedZroNode = new MSNode(format("Extended (%d)", bondedZro.size()));
        multiRes.getActive().getTerms().add(bondedZroNode);
        
//        printTreeFromNode(resOne);
//        printTreeFromNode(resZro);
        updateBondedLambdas();        
        ready = true;
        describe();
    }
    
    /**
     * List all the atoms and bonded terms associated with each end state.
     */
    public void describe() {
        SB.logfn(this.toString());
//        SB.logfn("    Backbone");
//        for (Atom atom : backbone) {
//            SB.logfn(" %s", atom);
//        }
        SB.logfn("    State 1: Atoms");
        for (Atom atom : atomsOne) {
            SB.logfn(" %s", atom);
        }
        SB.logfn("    State 1: Bonded");
        for (MSNode term : resOne.getTerms().getChildList()) {
            SB.logfn("      %s", term);
        }
        SB.logfn("    State 0: Atoms");
        for (Atom atom : atomsZro) {
            SB.logfn(" %s", atom);
        }
        SB.logfn("    State 0: Bonded");
        for (MSNode term : resZro.getTerms().getChildList()) {
            SB.logfn("      %s", term);
        }
        SB.print();
    }
    
    public boolean isReady() {
        return ready;
    }
    
    public void updateBondedLambdas() {
        if (!scaleBondedTerms) {
            return;
        }
        double lambda = getLambda();
        for (BondedTerm bt1 : bondedOne) {
            bt1.setEsvLambda(lambda);
        }
        for (BondedTerm bt0 : bondedZro) {
            bt0.setEsvLambda(1.0 - lambda);
        }
    }
    public List<Atom> getAtomListZro() {
        return atomsZro;
    }
    
    
}
