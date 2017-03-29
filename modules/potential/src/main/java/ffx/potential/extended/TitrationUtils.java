package ffx.potential.extended;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedUtils;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.Mutation;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.Keyword;

import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.extended.SBLogger.SB;

/**
 * Helper methods to define titration-specific phenomena.
 */
@SuppressWarnings("serial")
public class TitrationUtils {
    private TitrationUtils() {} // utility class
    private static final Logger logger = Logger.getLogger(TitrationUtils.class.getName());
    
    /**
     * Advanced options to both DiscreteMCMD and DiscountPh.
     */
    public static class TitrationConfig {
        public final ContinuousSeedDistribution seedDistribution
                                                        = prop(ContinuousSeedDistribution.class, "cphmd-seedMode", ContinuousSeedDistribution.FLAT);
        public MCOverride mcOverride                    = prop(MCOverride.class, "cphmd-override", MCOverride.NONE);
        public final Snapshots snapshots                = prop(Snapshots.class, "cphmd-snapshots", Snapshots.NONE);
        public final HistidineMode histidineMode        = prop(HistidineMode.class, "cphmd-histidineMode", HistidineMode.HIE_ONLY);
        public final OptionalDouble referenceOverride   = prop("cphmd-referenceOverride", OptionalDouble.empty());
        public final double meltdownTemperature         = prop("cphmd-meltdownTemp", 6000.0);
        public final double warningTemperature          = prop("cphmd-warningTemp", 1000.0);
        public final boolean logTimings                 = prop("cphmd-logTimings", false);
        public final boolean titrateTermini             = prop("cphmd-termini", false);             // TODO Finish termini support.
        public final boolean zeroReferences             = prop("cphmd-zeroReferences", true);
        public final int debugLogLevel                  = prop("cphmd-debugLog", 0);
        public final boolean useConformationalBias      = prop("cphmd-cbmcRotamerMoves", false);    // untested
        public final boolean inactivateBackground       = prop("cphmd-inactivateBackground", false);
        public final boolean zeroReferenceEnergies =
                prop("cphmd-zeroReferences", false, "Zeroing all reference energies!");
        public final OptionalDouble refOverride = 
                prop("cphmd-refOverride", OptionalDouble.empty(), "Reference protonation energies overridden!");
        
        public void print() {
			ExtUtils.printConfigSet("Titration Config:", System.getProperties(), "cphmd");
        }
    }
    
    /*
     * Properties of the TitrationUtils class itself; must be set at JVM launch.
     */
    public static final boolean heavyStrandedDynamics      = prop("cphmd-heavyStrandedDynamics", true);
    public static final boolean threeStateHistidines       = prop("cphmd-threeState", false);   // not yet implemented
    public static final boolean threeStateCarboxylics      = prop("cphmd-threeState", false);   // not yet implemented
    
    /**
     * How DISCOUNT initializes lambda values at outset of continuous dynamics.
     */
    public enum ContinuousSeedDistribution {
        FLAT, BETA, BOLTZMANN, DIRAC_CURRENT, DIRAC_POINTFIVE;
    }
    
    /**
     * Global override of MC acceptance criteria.
     */
    public enum MCOverride {
        NONE, ACCEPT, REJECT, ONCE;
    }

    /**
     * Writes .s-[num] and .f-[num] files representing before/after MC move
     * structures. Note: The 'after' snapshot represents the change that was
     * PROPOSED, regardless of accept/reject.
     */
    public enum Snapshots {
        INTERLEAVED, SEPARATE, NONE;
    }

    public enum HistidineMode {
        HIE_ONLY, HID_ONLY, SINGLE, DOUBLE;
    }
    
    public static void renumberAtoms(MolecularAssembly mola) {
        Atom[] atoms = mola.getAtomArray();
        int setter = 0;
        for (Atom atom : atoms) {
            atom.setXyzIndex(setter++);
        }
    }
    
    /**
     * @see TitrationUtils::initEsvPreloadProperties
     */
    public static void initDiscountPreloadProperties(Double cutoffs) {
        initEsvPreloadProperties(cutoffs);
        /* Available options: Key  */
        System.getProperty("cphmd-seedMode");
        System.getProperty("cphmd-override");
        System.getProperty("cphmd-snapshotsType");
        System.getProperty("cphmd-histidineMode");
        System.getProperty("cphmd-debugLog");
        System.getProperty("cphmd-referenceOverride");
        System.getProperty("cphmd-tempMonitor");
        System.getProperty("cphmd-logTimings");
        System.getProperty("cphmd-termini");
        System.getProperty("cphmd-zeroReferences");
    }
    public static void initDiscountPreloadProperties() {
        initDiscountPreloadProperties(null);
    }
    
    /**
     * Note that this must (generally) be called before loading the input file
     * or instantiating titration classes.
     */
    public static void initEsvPreloadProperties(Double cutoffs) {
        // Active Potential
        System.setProperty("forcefield", "AMOEBA_PROTEIN_2013");
        System.setProperty("esvterm", "true");
        System.setProperty("lambdaterm", "true");
        System.setProperty("bondterm", "true");
        System.setProperty("angleterm", "true");
        System.setProperty("strbndterm", "true");
        System.setProperty("ureyterm", "true");
        System.setProperty("opbendterm", "true");
        System.setProperty("torsionterm", "true");
        System.setProperty("pitorsterm", "true");
        System.setProperty("tortorterm", "true");
        System.setProperty("improperterm", "true");

        // Optional Potential
        System.setProperty("vdwterm", "true");         // van der Waals
        System.setProperty("esv-vdw", "true");
        System.setProperty("mpoleterm", "true");       // permanent real space
        System.setProperty("pme-qi", "true");
        System.setProperty("esv-pme", "true");
        System.setProperty("recipterm", "true");       // permanent reciprocal space

        // Inactive Potential
        System.setProperty("polarizeterm", "false");   // polarization
        System.setProperty("polarization", "NONE");
        System.setProperty("gkterm", "false");
        System.setProperty("restrainterm", "false");
        System.setProperty("comrestrainterm", "false");
        System.setProperty("lambda_torsions", "false");

        // Potential Settings
        System.setProperty("permanent-lambda-alpha","2.0");
        System.setProperty("permanent-lambda-exponent","3.0");
        System.setProperty("polarization-lambda-start","0.0");      // polarize on the whole range [0,1]
        System.setProperty("polarization-lambda-exponent","0.0");   // polarization not softcored, only prefactored
        System.setProperty("ligand-vapor-elec", "false");           // cancels when reference is solution phase
        System.setProperty("no-ligand-condensed-scf", "false");     // don't need condensed phase polarization
        System.setProperty("intramolecular-softcore", "true");
        System.setProperty("intermolecular-softcore", "true");
        if (cutoffs != null) {
            System.setProperty("vdw-cutoff", String.valueOf(cutoffs));
            System.setProperty("ewald-cutoff", String.valueOf(cutoffs));
        }

        // ESV Settings
        System.setProperty("esv-biasTerm", "true");             // include discretization and pH biases
        System.setProperty("esv-scaleBonded", "true");          // include effects on bonded terms
        System.setProperty("esv-backgroundBonded", "true");     // hook up BG bonded terms to FG node
        System.setProperty("esv-scaleUnshared", "true");        // use multipole scaling in all cases (eliminates softcoring)
    }
    public static void initEsvPreloadProperties() {
        initEsvPreloadProperties(null);
    }
    
    public static void activateResidue(Residue addDoF) {
        List<Atom> atomList = addDoF.getAtomList();
        for (Atom atom : atomList) {
            atom.setActive(true);
        }
    }
    
    public static void inactivateResidue(Residue killDoF) {
        List<Atom> atomList = killDoF.getAtomList();
        for (Atom atom : atomList) {
            atom.setActive(false);
        }
    }
    
    /**
     * Perform the requested titration on the given MultiResidue.
     * Remember to call reInit() on affected ForceFieldEnergy and MolecularDynamics objects!
     *
     * @param multiRes
     * @param titration
     */
    public static TitrationType performTitration(MultiResidue multiRes, Titration titration, boolean inactivateBackground) {
        AminoAcid3 current = multiRes.getActive().getAminoAcid3();
        final TitrationType type;
        final AminoAcid3 target;
        if (current == titration.protForm) {
            type = TitrationType.DEPROT;
            target = titration.deprotForm;
        } else if (current == titration.deprotForm) {
            type = TitrationType.PROT;
            target = titration.protForm;
        } else {
            throw new IllegalStateException();
        }
        boolean success = multiRes.requestSetActiveResidue(target);
        if (!success) {
            logger.severe(String.format("Couldn't perform requested titration for MultiRes: %s", multiRes.toString()));
        }
        
        List<Atom> oldAtoms = multiRes.getActive().getAtomList();
        List<Atom> newAtoms = multiRes.getActive().getAtomList();

        // identify which atoms were actually inserted/removed
        List<Atom> removedAtoms = new ArrayList<>();
        List<Atom> insertedAtoms = new ArrayList<>();
        for (Atom oldAtom : oldAtoms) {
            boolean found = false;
            for (Atom newAtom : newAtoms) {
                if (newAtom == oldAtom ||
                        (newAtom.getResidueNumber() == oldAtom.getResidueNumber()
                        && newAtom.getName().equals(oldAtom.getName()))) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                removedAtoms.add(oldAtom);
            }
        }
        for (Atom newAtom : newAtoms) {
            boolean found = false;
            for (Atom oldAtom : oldAtoms) {
                if (newAtom == oldAtom ||
                        (newAtom.getResidueNumber() == oldAtom.getResidueNumber()
                        && newAtom.getName().equals(oldAtom.getName()))) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                insertedAtoms.add(newAtom);
            }
        }
        if (insertedAtoms.size() + removedAtoms.size() > 1) {
            logger.warning("Protonate: removed + inserted atom count > 1.");
        }
        
        if (inactivateBackground) {
            activateResidue(multiRes.getActive());
            for (Residue res : multiRes.getInactive()) {
                inactivateResidue(res);
            }
        }
        
        return type;
    }
    
    public static MolecularAssembly openFullyProtonated(File structure) {
        String name = format("%s-prot", FilenameUtils.removeExtension(structure.getName()));
        MolecularAssembly mola = new MolecularAssembly(name);
        mola.setFile(structure);

        List<Mutation> mutations = new ArrayList<>();
        List<Residue> residues = mola.getResidueList();
        for (Residue res : residues) {
            char chain = res.getChainID();
            int resID = res.getResidueNumber();
            Titration titration = Titration.lookup(res);
            if (res.getAminoAcid3() != titration.protForm) {
                String protName = titration.protForm.name();
                mutations.add(new PDBFilter.Mutation(chain, resID, protName));
            }
        }
        
        PotentialsUtils utils = new PotentialsUtils();
        return utils.openWithMutations(structure, mutations);
    }
    public static MolecularAssembly openFullyProtonated(String filename) {
        return openFullyProtonated(new File(filename));
    }
    
    /**
     * Create a MultiResidue from the given Residue by adding its alternated protonation
     * state(s) as alternate possibilities.
     */
    public static MultiResidue titratingMultiresidueFactory(MolecularAssembly mola, Residue res) {
        ForceField ff = mola.getForceField();
        Potential potential = mola.getPotentialEnergy();
        if (!(potential instanceof ForceFieldEnergy)) {
            SB.warning("TitrationFactory only supported by ForceFieldEnergy potentials.");
            throw new IllegalStateException();
        }
        ForceFieldEnergy ffe = (ForceFieldEnergy) potential;
		
        /* Create new titration state. */
        Titration titration = Titration.lookup(res);
        String targetName = (titration.protForm != res.getAminoAcid3())
                ? titration.protForm.toString() : titration.deprotForm.toString();
        int resNumber = res.getResidueNumber();
        Residue.ResidueType resType = res.getResidueType();
        Residue newRes = new Residue(targetName, resNumber, resType);
		
        /* Wrap both states in a MultiResidue. */
        MultiResidue multiRes = new MultiResidue(res, ff, ffe);
        Polymer polymer = findResiduePolymer(res, mola);
        polymer.addMultiResidue(multiRes);
        multiRes.addResidue(newRes);
		
		/* Begin in protonated state by default. */
		multiRes.setActiveResidue(titration.protForm);
        propagateInactiveResidues(multiRes, false);
        ffe.reInit();
        return multiRes;
    }
	
    /**
     * Locate to which Polymer in a MolecularAssembly the given Residue belongs.
     */
    public static Polymer findResiduePolymer(Residue residue, MolecularAssembly mola) {
        if (residue.getChainID() == null) {
            logger.severe("No chain ID for residue " + residue);
        }
        Polymer polymers[] = mola.getChains();
        Polymer location = null;
        for (Polymer p : polymers) {
            if (p.getChainID().equals(residue.getChainID())) {
                location = p;
            }
        }
        if (location == null) {
            logger.severe("Couldn't find polymer for residue " + residue);
        }
        return location;
    }
    
    /**
     * Identify titratable residues and choose them all.
     */
    public static List<Residue> chooseTitratables(MolecularAssembly searchMe) {
        List<Residue> chosen = new ArrayList<>();
        Polymer polymers[] = searchMe.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                Residue res = residues.get(j);
                Titration[] avail = Titration.multiLookup(res);
                if (avail != null) {
                    chosen.add(residues.get(j));
                }
            }
        }
        return chosen;
    }
	
    /**
     * Choose titratables with intrinsic pKa inside (pH-window,pH+window).
     *
     * @param pH
     * @param window
     */
    public static List<Residue> chooseTitratables(double pH, double window, MolecularAssembly searchMe) {
        List<Residue> chosen = new ArrayList<>();
        Polymer polymers[] = searchMe.getChains();
        for (int i = 0; i < polymers.length; i++) {
            ArrayList<Residue> residues = polymers[i].getResidues();
            for (int j = 0; j < residues.size(); j++) {
                Residue res = residues.get(j);
                Titration[] avail = Titration.multiLookup(res);
                for (Titration titration : avail) {
                    double pKa = titration.pKa;
                    if (pKa >= pH - window && pKa <= pH + window) {
                        chosen.add(residues.get(j));
                    }
                }
            }
        }
        return chosen;
    }
	
    public static List<Residue> chooseTitratables(String residueIDs, MolecularAssembly searchMe) {
        String[] tokens =
                  (residueIDs.split(".").length > 1) ? residueIDs.split(".")
                : (residueIDs.split(",").length > 1) ? residueIDs.split(",")
                : new String[]{residueIDs};
        return chooseTitratables(Arrays.asList(tokens), searchMe);

    }
    
    /**
     * Select titrating residues by amino acid.
     */
    public static List<Residue> chooseTitratables(AminoAcid3 aa, MolecularAssembly searchMe) {
        List<Residue> chosen = new ArrayList<>();
        Polymer polymers[] = searchMe.getChains();
        for (Polymer polymer : polymers) {
            ArrayList<Residue> residues = polymer.getResidues();
            for (Residue res : residues) {
                if (res.getAminoAcid3() == aa) {
                    Titration[] avail = Titration.multiLookup(res);
                    if (avail != null) {
                        chosen.add(res);
                    }
                }
            }
        }
        return chosen;
    }
	
    public static List<Residue> chooseTitratables(List<String> crIDs, MolecularAssembly searchMe) {
        List<Residue> chosen = new ArrayList<>();
        for (String crID : crIDs) {
            char chain = crID.charAt(0);
            int num = Integer.parseInt(crID.substring(1));
            boolean found = false;
            List<Residue> allRes = searchMe.getResidueList();
            for (Residue res : allRes) {
                if (res.getChainID() == chain && res.getResidueNumber() == num) {
                    chosen.add(res);
                    found = true;
                    break;
                }
            }
            if (!found) {
                logger.severe(format("Couldn't find residue for crID %c,%d.", chain, num));
            }
        }
        return chosen;
    }
	
    public static List<Residue> chooseTitratables(char chain, int resID, MolecularAssembly searchMe) {
        List<Residue> chosen = new ArrayList<>();
        Polymer polymers[] = searchMe.getChains();
        for (Polymer polymer : polymers) {
            if (polymer.getChainID() == chain) {
                ArrayList<Residue> residues = polymer.getResidues();
                for (Residue residue : residues) {
                    if (residue.getResidueNumber() == resID) {
                        chosen.add(residue);
                        logger.info(String.format(" Chosen: %s", residue));
                    }
                }
            }
        }
        return chosen;
    }
    
    /**
     * @see TitrationUtils::propagateInactiveResidues(MultiResidue, boolean)
     */
    public static void propagateInactiveResidues(List<MultiResidue> multiResidues, boolean propagateDynamics) {
        for (MultiResidue multiRes : multiResidues) {
            propagateInactiveResidues(multiRes, propagateDynamics);
        }
    }
    /**
     * @see TitrationUtils::propagateInactiveResidues(MultiResidue, boolean)
     */
    public static void propagateInactiveResidues(List<MultiResidue> multiResidues) {
        for (MultiResidue multiRes : multiResidues) {
            propagateInactiveResidues(multiRes, true);
        }
    }
    /**
     * @see TitrationUtils::propagateInactiveResidues(MultiResidue, boolean)
     */
    public static void propagateInactiveResidues(MultiResidue multiResidue) {
        propagateInactiveResidues(multiResidue, true);
    }    
    
    /**
     * Copies atomic coordinates from each active residue to its inactive
     * counterparts. Inactive hydrogen coordinates are updated by geometry
     * with the propagated heavies.
     */
    public static void propagateInactiveResidues(MultiResidue multiRes, boolean propagateDynamics) {
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
                    double[] activeXYZ = activeAtom.getXYZ(null);
                    inactiveAtom.setXYZ(activeXYZ);
                    double[] grad = new double[3];
                    activeAtom.getXYZGradient(grad);
                    inactiveAtom.setXYZGradient(grad[0], grad[1], grad[2]);
                    if (propagateDynamics) {
                        // Propagate velocity, acceleration, and previous acceleration.
                        double[] activeVelocity = new double[3];
                        activeAtom.getVelocity(activeVelocity);
                        inactiveAtom.setVelocity(activeVelocity);
                        double[] activeAccel = new double[3];
                        activeAtom.getAcceleration(activeAccel);
                        inactiveAtom.setAcceleration(activeAccel);
                        double[] activePrevAcc = new double[3];
                        activeAtom.getPreviousAcceleration(activePrevAcc);
                        inactiveAtom.setPreviousAcceleration(activePrevAcc);
                    }
                } else {
                    if (activeName.equals("C") || activeName.equals("O") || activeName.equals("N") || activeName.equals("CA") || activeName.equals("H") || activeName.equals("HA")) {
                        // Backbone atoms aren't supposed to exist in inactive multiResidue components; so no problem.
                    } else if (isTitratableHydrogen(activeAtom)) {
                        /** i.e.
                        ((activeResName.equals("LYS") && activeName.equals("HZ3"))
                        || (activeResName.equals("TYR") && activeName.equals("HH"))
                        || (activeResName.equals("CYS") && activeName.equals("HG"))
                        || (activeResName.equals("HIS") && (activeName.equals("HD1") || activeName.equals("HE2")))
                        || (activeResName.equals("HID") && activeName.equals("HD1"))
                        || (activeResName.equals("HIE") && activeName.equals("HE2"))
                        || (activeResName.equals("ASH") && activeName.equals("HD2"))
                        || (activeResName.equals("GLH") && activeName.equals("HE2")))   */
                        // These titratable protons are handled below; so no problem.
                    } else {
                        // Now we have a problem.
                        SB.warning("Couldn't propagate inactive MultiResidue atom: %s: %s, %s", multiRes, activeName, activeAtom);
                    }
                }
            }
        }
        rebuildStrandedProtons(multiRes);
    }

    public static boolean isTitratableHydrogen(Atom atom) {
        String name = atom.getName();
        switch (atom.getResidueName()) {
            case "LYS":
                if (name.equals("HZ3")) {
                    return true;
                }
                break;
            case "TYR":
                if (name.equals("HH")) {
                    return true;
                }
                break;
            case "CYS":
                if (name.equals("HG")) {
                    return true;
                }
                break;
            case "HIS":
                if (name.equals("HD1") || name.equals("HE2")) {
                    return true;
                }
                break;
            case "HID":
                if (name.equals("HD1")) {
                    return true;
                }
                break;
            case "HIE":
                if (name.equals("HE2")) {
                    return true;
                }
                break;
            case "ASH":
                if (name.equals("HD2")) {
                    return true;
                }
                break;
            case "GLH":
                if (name.equals("HE2")) {
                    return true;
                }
                break;
        }
        return false;
    }

    /**
     * Rebuild stranded titratable protons from ideal geometry.
     * "Stranded protons" are titrating H+ atoms on inactive MultiRes members; when
     * propagating new coordinates to inactive residues, no coords/velocity exist for them.
     */
    private static void rebuildStrandedProtons(MultiResidue multiRes) {
        // If inactive residue is a protonated form, move the stranded hydrogen to new coords (based on propagated heavies).
        // Also give the stranded hydrogen a maxwell velocity and remove its accelerations.
        List<Residue> inactives = multiRes.getInactive();
        for (Residue inactive : inactives) {
            List<Atom> resetMe = new ArrayList<>();
            switch (inactive.getName()) {
                case "LYS":
                    {
                        Atom HZ3 = (Atom) inactive.getAtomNode("HZ3");
                        Atom NZ = (Atom) inactive.getAtomNode("NZ");
                        Atom CE = (Atom) inactive.getAtomNode("CE");
                        Atom HZ1 = (Atom) inactive.getAtomNode("HZ1");
                        BondedUtils.intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
                        resetMe.add(HZ3);
                        break;
                    }
                case "ASH":
                    {
                        Atom HD2 = (Atom) inactive.getAtomNode("HD2");
                        Atom OD2 = (Atom) inactive.getAtomNode("OD2");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom OD1 = (Atom) inactive.getAtomNode("OD1");
                        BondedUtils.intxyz(HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0);
                        resetMe.add(HD2);
                        break;
                    }
                case "GLH":
                    {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom OE2 = (Atom) inactive.getAtomNode("OE2");
                        Atom CD = (Atom) inactive.getAtomNode("CD");
                        Atom OE1 = (Atom) inactive.getAtomNode("OE1");
                        BondedUtils.intxyz(HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0);
                        resetMe.add(HE2);
                        break;
                    }
                case "HIS":
                    {
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
                case "HID":
                    {
                        Atom HD1 = (Atom) inactive.getAtomNode("HD1");
                        Atom ND1 = (Atom) inactive.getAtomNode("ND1");
                        Atom CG = (Atom) inactive.getAtomNode("CG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
                        resetMe.add(HD1);
                        break;
                    }
                case "HIE":
                    {
                        Atom HE2 = (Atom) inactive.getAtomNode("HE2");
                        Atom NE2 = (Atom) inactive.getAtomNode("NE2");
                        Atom CD2 = (Atom) inactive.getAtomNode("CD2");
                        Atom CE1 = (Atom) inactive.getAtomNode("CE1");
                        BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
                        resetMe.add(HE2);
                        break;
                    }
                case "CYS":
                    {
                        Atom HG = (Atom) inactive.getAtomNode("HG");
                        Atom SG = (Atom) inactive.getAtomNode("SG");
                        Atom CB = (Atom) inactive.getAtomNode("CB");
                        Atom CA = (Atom) inactive.getAtomNode("CA");
                        BondedUtils.intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
                        resetMe.add(HG);
                        break;
                    }
                case "TYR":
                    {
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
                if (heavyStrandedDynamics) {
                    // Use of heavy atom dynamics properties is in testing.
                    a.setXYZGradient(0, 0, 0);
                    double[] heavyVelocity = new double[3];
                    double[] heavyAccel = new double[3];
                    double[] heavyPrevAccel = new double[3];
                    Atom heavy = a.getBonds().get(0).get1_2(a);
                    heavy.getVelocity(heavyVelocity);
                    heavy.getAcceleration(heavyAccel);
                    heavy.getPreviousAcceleration(heavyPrevAccel);
                    a.setVelocity(heavyVelocity);
                    a.setAcceleration(heavyAccel);
                    a.setPreviousAcceleration(heavyPrevAccel);
                } else {
                    // PREVIOUSLY: draw vel from maxwell and set accel to zero
                    a.setXYZGradient(0, 0, 0);
                    a.setVelocity(ExtUtils.singleRoomtempMaxwell(a.getMass()));
                    a.setAcceleration(new double[]{0, 0, 0});
                    a.setPreviousAcceleration(new double[]{0, 0, 0});
                }
            }
        }
    }
    

    /**
     * Amino acid protonation reactions.
     * Constructors below specify intrinsic pKa and reference free energy of protonation,
     * obtained via (OSRW) metadynamics on capped monomers.
     */
    public enum Titration {
        ctoC(   8.18,  60.168, AminoAcid3.CYD, AminoAcid3.CYS),
        Dtod(   3.90,  53.188, AminoAcid3.ASP, AminoAcid3.ASH),
        Etoe(   4.25,  59.390, AminoAcid3.GLU, AminoAcid3.GLH),
        ktoK(  10.53, -50.440, AminoAcid3.LYD, AminoAcid3.LYS),
        ytoY(  10.07,  34.961, AminoAcid3.TYD, AminoAcid3.TYR),
        UtoH(   6.00, -42.923, AminoAcid3.HID, AminoAcid3.HIS),
        ZtoH(   6.00,  00.000, AminoAcid3.HIE, AminoAcid3.HIS),
        TerminalNH3toNH2(   8.23,  00.00, AminoAcid3.UNK, AminoAcid3.UNK),
        TerminalCOOHtoCOO(  3.55,  00.00, AminoAcid3.UNK, AminoAcid3.UNK);
        
        public final double pKa;
        public final double refEnergy;
        public final AminoAcid3 protForm;
        public final AminoAcid3 deprotForm;

        /**
         * Invoked by Enum; use the factory method to obtain instances.
         */
        private Titration(double pKa, double refEnergy, AminoAcid3 deprotForm, AminoAcid3 protForm) {
            this.pKa = pKa;
            this.refEnergy = refEnergy;
            this.deprotForm = deprotForm;
            this.protForm = protForm;
        }
        
        public static Titration lookup(Residue res) {
            Titration[] titrations = multiLookup(res);
            if (titrations.length > 1) {
                SB.warning("Titration::lookup returned more results than expected. Did you mean to invoke multi-state?");
            }
            return (titrations != null) ? titrations[0] : null;
        }
        
        /**
         * Return a Titration object for the given Residue.
         * TODO: Add support for multi-state titrations (HIS,ASP,GLU).
         */
        public static Titration[] multiLookup(Residue res) {
            List<Titration> titrations = new ArrayList<>();
            AminoAcid3 current = AminoAcid3.valueOf(res.getName());
            if (threeStateHistidines) {
                if (current == AminoAcid3.HIS || current == AminoAcid3.HID || current == AminoAcid3.HIE) {
                    return new Titration[]{ZtoH, UtoH};
                }
            }
            if (threeStateCarboxylics) {
                if (current == AminoAcid3.ASP || current == AminoAcid3.ASH) {
//                    return new Titration[]{Dtod, Dtod2};
                }
                if (current == AminoAcid3.GLU || current == AminoAcid3.GLH) {
//                    return new Titration[]{Etoe, Etoe2};
                }
            }
            for (Titration titration : Titration.values()) {
                if (current == titration.protForm || current == titration.deprotForm) {
                    return new Titration[]{titration};
                }
            }
            SB.warning("No titration lookup found for residue %s", res);
            return null;
        }
    }
    
    public enum TitrationType {
        PROT, DEPROT;
    }

}
