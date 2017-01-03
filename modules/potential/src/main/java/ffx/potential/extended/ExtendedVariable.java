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
import ffx.potential.nonbonded.MultiplicativeSwitch;

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
    protected static final boolean esvPropagation = prop("esv-propagation", false);
    protected static final boolean scaleBondedTerms = prop("esv-scaleBonded", false);
    protected static final OptionalDouble biasOverride = prop("esv-bias", OptionalDouble.empty());
    private final double thetaMass = prop("esv-thetaMass", 1.0e-18);            // from OSRW, reasonably 100 a.m.u.
    private final double thetaFriction = prop("esv-thetaFriction", 1.0e-19);    // from OSRW, reasonably 60/ps
    private final ExtInclusionMode extInclusion = 
            prop(ExtInclusionMode.class, "esv-atomInclusion", ExtInclusionMode.TITRATABLE_H);
    public enum ExtInclusionMode {
        ALL_ATOMS, TITRATABLE_H;
    }
    
    // Lambda and derivative variables
    private double lambda;                          // ESVs travel on {0,1}
    private double theta;                           // Propagates lambda particle via "lambda=sin(theta)^2"
    private double halfThetaVelocity = 0.0;         // from OSRW, start theta with zero velocity
    private final Random stochasticRandom = ThreadLocalRandom.current();
    /**
     * Magnitude of the discretization bias in kcal/mol.
     */
    private final double discrBiasBeta;
    /**
     * Discretization bias and its (chain rule) derivative.
     */
    private double discrBias, dDiscrBiasdL;
    /**
     * Switched lambda value and its (chain rule) derivative.
     */
    private double lSwitch, dlSwitch;
    /**
     * Sigmoidal switching function. Maps lambda -> S(lambda) which has a flatter deriv near zero/unity.
     */
    private final MultiplicativeSwitch switchingFunction;
    
    // Atom lists and scaled terms
    private MultiResidue multiRes;
    private List<Atom> backbone = new ArrayList<>();
    private Residue resOne, resZro;                 // resOne*lamedh + resZero*(1-lamedh)
    private List<Atom> atomsOne, atomsZro;          // side-chain only; atomsZro stay unconnected to mola
    private List<Atom> atomsShared, atomsUnshared;
    private List<BondedTerm> bondedOne, bondedZro;  // valence terms for each side; mola won't see zro by default
    private MSNode termNode;                        // modified to contain all applicable bonded terms
    private int moleculeNumber = 0;
    
    public ExtendedVariable(MultiResidue multiRes, double biasMag, double initialLambda) {
        index = esvIndexer++;
        discrBiasBeta = biasOverride.isPresent() ? biasOverride.getAsDouble() : biasMag;
        this.switchingFunction = new MultiplicativeSwitch(0.0, 1.0);
        setLambda(initialLambda);

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
    
    /**
     * Called by readyUp() to populate atomsOne, atomsZro, atomsShared, atomsUnshared.
     * Use this to move titration-specific initialization to TitrationESV.
     */
//    public abstract void fillAtomLists();
    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double getTotalBias(double temperature, boolean print);
    /**
     * Should include at least the discretization bias; add any type-specific biases (eg pH).
     */
    public abstract double getTotalBiasDeriv(double temperature, boolean print);
        
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
        this.lSwitch = switchingFunction.taper(lambda);
        this.dlSwitch = switchingFunction.dtaper(lambda);
        theta = Math.asin(Math.sqrt(lambda));
        discrBias = discrBiasBeta - 4*discrBiasBeta*(lambda-0.5)*(lambda-0.5);
        dDiscrBiasdL = -8*discrBiasBeta*(lambda-0.5);
        updateBondedLambdas();
    }
    
    /**
     * The unswitched lambda value, ie input to S(L).
     * This is probably not what you want.
     */
    public final double getLambda() {
        return lambda;      // L
    }

    public final double getLambdaSwitch() {
        return lSwitch;     // S(L)
    }
    
    public final double getSwitchDeriv() {
        return dlSwitch;    // dS(L)dL
    }
    
    public final int getIndex() {
        return index;
    }
    
    @Override
    public String toString() {
        return String.format(" ESV%d:(%4.2f->%4.2f)", 
                index, getLambda(), getLambdaSwitch());
    }
    
    /**
     * From Shen&Huang 2016; drives ESVs to zero/unity.
     * bias = 4B*(L-0.5)^2
     */
    public double getDiscrBias() {
        return discrBias;
    }
    
    /**
     * dBiasdL = -8B*(L-0.5)
     */
    public double getDiscrBiasDeriv() {
        return dDiscrBiasdL;
    }
    
    /**
     * Fill the atom arrays; apply persistent indexing; set atom esv properties;
     * fill the bonded term arrays; set esv lambda on bonded terms.
     */
    public void readyup() {        
        // Fill the atom lists.
        atomsOne = new ArrayList<>();
        atomsZro = new ArrayList<>();
        atomsShared = new ArrayList<>();
        atomsUnshared = new ArrayList<>();
        for (String bbName : ExtConstants.backboneNames) {
            Atom bb = (Atom) resOne.getAtomNode(bbName);
            if (bb != null) {
                backbone.add(bb);
                bb.applyPersistentIndex();
            }
        }
        for (Atom a1 : resOne.getAtomList()) {
            if (!backbone.contains(a1)) {
                a1.setESV(this);
                a1.applyPersistentIndex();
                if ((extInclusion == ExtInclusionMode.TITRATABLE_H && isTitratableHydrogen(a1))
                        || extInclusion == ExtInclusionMode.ALL_ATOMS) {
                    atomsOne.add(a1);
                    a1.setEsvState(1);
                    atomsUnshared.add(a1);
                } else {
                    atomsShared.add(a1);
                }
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
                    atomsUnshared.add(a0);
                }
            }
        }
        
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
        describe(ExtUtils.DebugHandler.VERBOSE());
    }
    
    /**
     * List all the atoms and bonded terms associated with each end state.
     */
    public void describe(boolean verbose) {
        if (!ready) {
            return;
        }
        if (!verbose) {
            SB.logfn(this.toString());
        // Print the switching function.
            SB.logfn(" Switching on (%.2f,%.2f) with e0,de0,e1,de1: %.2f %.2f %.2f %.2f",
                    0.0, 1.0,
                    switchingFunction.taper(0.0), switchingFunction.dtaper(0.0),
                    switchingFunction.taper(1.0), switchingFunction.dtaper(1.0));
    //        SB.logfn("    Backbone");
    //        for (Atom atom : backbone) {
    //            SB.logfn(" %s", atom);
    //        }
            SB.logfn("    Shared Atoms");
            for (Atom atom : atomsShared) {
                SB.logfn(" %s", atom);
            }
            SB.logfn("    Unshared Atoms");
            for (Atom atom : atomsUnshared) {
                SB.logfn(" %s", atom);
            }
            SB.logfn("    Bonded");
            for (MSNode term : resOne.getTerms().getChildList()) {
                SB.logfn("      %s", term);
                if (term.getName().trim().contains("Extended")) {
                    for (MSNode ext : term.getChildList()) {
                        SB.logfn("        %s", ext);
                    }
                }
            }
        } else {
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
        }
        SB.print();
    }
    
    public List<Atom> getUnsharedAtoms() {
        List<Atom> ret = new ArrayList<>();
        ret.addAll(atomsUnshared);
        return ret;
    }
    
    public boolean isReady() {
        return ready;
    }
    
    public void updateBondedLambdas() {
        if (!scaleBondedTerms) {
            return;
        }
        double Sl = getLambdaSwitch();
        double dSldL = getSwitchDeriv();
        for (BondedTerm bt1 : bondedOne) {
            bt1.setEsvLambda(Sl, dSldL);
        }
        for (BondedTerm bt0 : bondedZro) {
            bt0.setEsvLambda(1.0 - Sl, -dSldL);
        }
    }
    public List<Atom> getAtomListZro() {
        return atomsZro;
    }
    
    public int getMoleculeNumber() {
        return moleculeNumber;
    }
    
}
