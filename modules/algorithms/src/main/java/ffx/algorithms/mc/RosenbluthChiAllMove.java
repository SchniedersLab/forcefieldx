/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.algorithms.mc;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Torsion;
import ffx.potential.parsers.PDBFilter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;

/**
 * Represents a Boltzmann-drawn spin of all residue torsions.
 * For use with RosenbluthCBMC (configurational-bias Monte Carlo).
 * Biases each torsion by drawing test set from Boltzmann distr on torsion energy.
 * Selects from amongst each test set on Boltzmann weight of remaining energy.
 * Calculates Rosenbluth factors along the way; acceptance criterion = Wn/Wo.
 * ----
 * Note:
 * Contains much of the infrastructure necessary to generalize toward full 
 * polymer-construction-type CBMC (i.e. changing angles and bonds as well).
 * This implementation leaves bonds/angles fixed and considers only 
 * torsion energy as dependent (Ubond in Frenkel/Smit 13.2).
 * @author S. LuCore
 */
public class RosenbluthChiAllMove implements MCMove {
    private static final Logger logger = Logger.getLogger(RosenbluthChiAllMove.class.getName());
    private final double BOLTZMANN = 0.0019872041; // In kcal/(mol*K)
    
    private final Residue target;
    private final ResidueState origState;
    private Rotamer proposedMove;
    private double Wn = 0.0;
    private double Wo = 0.0;
    private final int testSetSize;
    private final ForceFieldEnergy ffe;
    private final double beta;
    private final ThreadLocalRandom rand = ThreadLocalRandom.current();
    private final StringBuilder report = new StringBuilder();
    private boolean verbose = false;
    private SnapshotWriter snapshotWriter = null;
    private int moveNumber = 0;
    
    public RosenbluthChiAllMove(Residue target, int testSetSize, 
            ForceFieldEnergy ffe, double temperature) {
        this.target = target;
        this.testSetSize = testSetSize;
        this.ffe = ffe;
        this.beta = 1 / (BOLTZMANN * temperature);
        origState = target.storeState();
        engage();
    }
    
    public RosenbluthChiAllMove(MolecularAssembly mola, Residue target, 
            int testSetSize, ForceFieldEnergy ffe, double temperature, 
            boolean writeSnapshots, int moveNumber, boolean verbose) {
        this.target = target;
        this.testSetSize = testSetSize;
        this.ffe = ffe;
        this.beta = 1 / (BOLTZMANN * temperature);
        this.moveNumber = moveNumber;
        this.verbose = verbose;
        origState = target.storeState();
        if (writeSnapshots) {
            snapshotWriter = new SnapshotWriter(mola, false);
        }
        engage();
    }
    
    private void engage() {
        report.append(String.format(" Rosenbluth CBMC Move: %4d\n", moveNumber));
        report.append(String.format("    residue:   %s\n", target.toString()));
        writeSnapshot("orig", false);
        
        AminoAcid3 name = AminoAcid3.valueOf(target.getName());
        double chi[] = RotamerLibrary.measureRotamer(target, false);
        HashMap<Integer,BackBondedList> map = createBackBondedMap(name);
        
        // For each chi, create a test set from Boltzmann distr on torsion energy (Ubond).
        // Select from among this set based on Boltzmann weight of REMAINING energy (Uext).
        // The Uext partition function of each test set (wn) becomes a factor of the overall Rosenbluth (Wn).
        double wn[] = new double[chi.length];   // factors of Wn
        double finalChi[] = new double[chi.length];
        for (int i = 0; i < chi.length; i++) {
            Torsion tors = map.get(i).torsion;
            TrialSet trialSet = boltzmannTorsionSet(tors, i, testSetSize, "bkn");
            wn[i] = trialSet.sumExtBolt();
            if (i == 0) {
                Wn = wn[i];
            } else {
                Wn *= wn[i];
            }
            
            // Choose a proposal move from amongst this trial set (bn).
            double rng = rand.nextDouble(wn[i]);
            double running = 0.0;
            for (int j = 0; j < trialSet.uExt.length; j++) {
                double uExtBolt = FastMath.exp(-beta * trialSet.uExt[j]);
                running += uExtBolt;
                if (rng < running) {
                    finalChi[i] = trialSet.theta[j];    // Yes, I mean i then j.
                    double prob = uExtBolt / wn[i] * 100;
                    report.append(String.format("       Chose %d   %7.4f\t%7.4f\t  %4.1f%%\n", 
                            j, trialSet.uExt[j], uExtBolt, prob));
                    break;
                }
            }
        }
        proposedMove = createRotamer(name, finalChi);
        
        // Reprise the above procedure for the old configuration.
        // The existing conformation forms first member of each test set (wo).
        // Overall Rosenbluth Wo is product of uExt partition functions.
        double wo[] = new double[chi.length];   // factors of Wo
        for (int i = 0; i < chi.length; i++) {
            Torsion tors = map.get(i).torsion;
            double ouDep = tors.energy(false);      // original-conf uDep
            double ouExt = totalEnergy() - ouDep;   // original-conf uExt
            double ouExtBolt = FastMath.exp(-beta * ouExt);
            TrialSet trialSet = boltzmannTorsionSet(tors, i, testSetSize - 1, "bko");
            wo[i] = ouExtBolt + trialSet.sumExtBolt();
            if (i == 0) {
                Wo = wo[i];
            } else {
                Wo *= wo[i];
            }
        }
        
        target.revertState(origState);
        updateAll();
        if (verbose) {
            logger.info(report.toString());
        }
        
        // Logging yet to be added.
        if (false) {
//            report.append(String.format("    chi0:      %s\n", chi0.toString()));
//            writeSnapshot("uIndO");
//            writeSnapshot("uIndN");
//            report.append(String.format("    theta:     %3.2f\n", ((RosenbluthChi0Move) proposal).theta));
//            report.append(String.format("    criterion: %1.4f\n", criterion));
//            report.append(String.format("       Wn/Wo:     %.2f\n", Wn/Wo));
//            report.append(String.format("       uIndN,O:  %7.2f\t%7.2f\n", uIndN, uIndO));
//            report.append(String.format("       dInd(E):  %7.2f\t%7.2f\n", dInd, dIndE));
//            report.append(String.format("    rng:       %1.4f\n", rng));
        }
    }
    
    public double getWn() { return Wn; }
    public double getWo() { return Wo; }
    
    /**
     * Uses the accept-reject method (F/S Algorithm46) to draw new
     * chi values for the given torsion.
     */
    private TrialSet boltzmannTorsionSet(Torsion tors, int chiIndex, int setSize, String snapSuffix) {
        report.append(String.format("    TrialSet for Chi%d\t\t(Theta uDep uExt)\n", chiIndex));
        TrialSet trialSet = new TrialSet(setSize);
        double origChi[] = RotamerLibrary.measureRotamer(target, false);
        int i = 0;
        while (i < setSize) {
            double theta = rand.nextDouble(360.0) - 180;
            double newChi[] = new double[origChi.length];
            System.arraycopy(origChi, 0, newChi, 0, origChi.length);
            newChi[chiIndex] = theta;
            Rotamer newState = createRotamer(target, newChi);
            RotamerLibrary.applyRotamer(target, newState);
            double uBond = tors.energy(false);
            double criterion = FastMath.exp(-beta * uBond);
            double rng = rand.nextDouble();
            if (rng < criterion) {
                trialSet.theta[i] = theta;
                trialSet.rotamer[i] = newState;
                trialSet.uDep[i] = uBond;
                trialSet.uExt[i] = totalEnergy() - uBond;
                i++;
                writeSnapshot(snapSuffix, true);
                report.append(String.format("       %3s %d:      %5.2f\t%5.2f\t%5.2f\n",
                        snapSuffix, i, theta, trialSet.uDep[i], trialSet.uExt[i]));
            }
        }
        target.revertState(origState);
        updateAll();
        return trialSet;
    }
    
    private Rotamer createRotamer(AminoAcid3 name, double chi[]) {
        // Need to add sigma values to construct a new Rotamer with these chis.
        double values[] = new double[chi.length * 2];
        for (int k = 0; k < chi.length; k++) {
            int kk = 2*k;
            values[kk] = chi[k];
            values[kk+1] = 0.0;
        }
        return new Rotamer(name, values);
    }
    
    private Rotamer createRotamer(Residue res, double chi[]) {
        return createRotamer(AminoAcid3.valueOf(res.getName()), chi);
    }
    
    /**
     * Calculates all 'back-bonded' (ie toward the peptide backbone) 
     * energy dependent on a given chi.
     * @return 
     */
    private double backBondedEnergy(BackBondedList bbl) {
        double sum = 0.0;
        sum += bbl.bond.energy(false);
        sum += bbl.angle.energy(false);
        sum += bbl.torsion.energy(false);
        return sum;
    }
    
    /**
     * Yields a random vector on the surface of the unit sphere.
     * Algorithm 42 from Frenkel/Smit.
     */
    private double[] vectorOnASphere() {
        ThreadLocalRandom rand = ThreadLocalRandom.current();
        double ranA, ranB, ranC, ransq;
        do {
            ranA = 1 - 2*rand.nextDouble();
            ranB = 1 - 2*rand.nextDouble();
            ransq = ranA*ranA + ranB*ranB;
        } while (ransq >= 1);
        ranC = 2*FastMath.sqrt(1-ransq);
        double vec[] = new double[3];
        vec[0] = ranA*ranC;     // x
        vec[1] = ranB*ranC;     // y
        vec[2] = 1 - 2*ransq;   // z
        return vec;
    }
    
    /**
     * Performs the move associated with this MCMove.
     * Also updates chi values in associated Torsion objects.
     */
    @Override
    public void move() {
        RotamerLibrary.applyRotamer(target, proposedMove);
        updateAll();
    }
    
    /**
     * Reverts the last applied move() call.
     */
    @Override
    public void revertMove() {
        target.revertState(origState);
        updateAll();
    }
    
    private void updateAll() {
        updateBonds();
        updateAngles();
        updateTorsions();
    }
    
    private void updateBonds() {
        List<ROLS> bonds = target.getBondList();
        for (ROLS rols : bonds) {
            ((Bond) rols).update();
        }
//        bonds.stream().forEach(b -> ((Bond) b).update());
    }
    
    private void updateAngles() {
        List<ROLS> angles = target.getAngleList();
        for (ROLS rols : angles) {
            ((Angle) rols).update();
        }
//        angles.stream().forEach(a -> ((Angle) a).update());
    }
    
    private void updateTorsions() {
        List<ROLS> torsions = target.getTorsionList();
        for (ROLS rols : torsions) {
            ((Torsion) rols).update();
        }
//        torsions.stream().forEach(t -> ((Torsion) t).update());
    }
    
    @Override
    public String toString() {
        return String.format("Rosenbluth Rotamer Move:\n   Res:   %s\n   Rota: %s",
                target.toString(), proposedMove.toString());
    }
    
    private double totalEnergy() {
        double x[] = new double[ffe.getNumberOfVariables()*3];
        ffe.getCoordinates(x);
        return ffe.energy(x);
    }
    
    private class TrialSet {
        public final Rotamer rotamer[];
        public final double uDep[];
        public final double uExt[];
        public final double theta[];
        
        public TrialSet(int setSize) {
            rotamer = new Rotamer[setSize];
            uDep = new double[setSize];
            uExt = new double[setSize];
            theta = new double[setSize];
        }
        
        public double sumExtBolt() {
            double sum = 0.0;
            for (int i = 0; i < uExt.length; i++) {
                sum += FastMath.exp(-beta * uExt[i]);
            }
            return sum;
        }
    }
    
    private class BackBondedList {
        public final Bond bond;
        public final Angle angle;
        public final Torsion torsion;
        
        public BackBondedList(Bond bond, Angle angle, Torsion tors) {
            this.bond = bond;
            this.angle = angle;
            this.torsion = tors;
        }
    }
    
    /**
     * Maps the back-bonded terms affected by key atoms in an amino acid.
     * Here, 'key atom' refers to each new rotamer-torsion-completing atom.
     * e.g. VAL has 1 key atom (CG1), ARG has 4 key atoms (CG,CD,NE,CZ).
     * 'Back-bonded' means we only map terms that lead toward the backbone.
     */
    private HashMap<Integer,BackBondedList> createBackBondedMap(AminoAcid3 name) {
        HashMap<Integer,BackBondedList> map = new HashMap<>();
        List<Atom> chain = new ArrayList<>();
        Atom N = (Atom) target.getAtomNode("N");
        Atom CA = (Atom) target.getAtomNode("CA");
        Atom CB = (Atom) target.getAtomNode("CB");
        List<Atom> keyAtoms = new ArrayList<>();
        switch (name) {
            case VAL: {
                Atom CG1 = (Atom) target.getAtomNode("CG1");
                keyAtoms.add(CG1);
                keyAtoms.add(CB);
                break;
            }
            case LEU: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                keyAtoms.add(CG);
                keyAtoms.add(CD1);
                break;
            }
            case ILE: {
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                Atom CG1 = (Atom) target.getAtomNode("CG1");
                keyAtoms.add(CD1);
                keyAtoms.add(CG1);
                break;
            }
            case SER: {
                Atom OG = (Atom) target.getAtomNode("OG");
                Atom HG = (Atom) target.getAtomNode("HG");
                keyAtoms.add(OG);
                keyAtoms.add(HG);
                break;
            }
            case THR: {
                Atom OG1 = (Atom) target.getAtomNode("OG1");
                Atom HG1 = (Atom) target.getAtomNode("HG1");
                keyAtoms.add(OG1);
                keyAtoms.add(HG1);
                break;
            }
            case CYX:
            case CYD: {
                Atom SG = (Atom) target.getAtomNode("SG");
                keyAtoms.add(SG);
                break;
            }
            case PHE: {
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                Atom CG = (Atom) target.getAtomNode("CG");
                keyAtoms.add(CG);
                break;
            }
            case PRO: {
                // Not allowed yet.
                Atom CD = (Atom) target.getAtomNode("CD");
                Atom CG = (Atom) target.getAtomNode("CG");
                keyAtoms.add(CG);
                keyAtoms.add(CD);
                break;
            }
            case TYR: {
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                Atom CE2 = (Atom) target.getAtomNode("CE2");
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CZ = (Atom) target.getAtomNode("CZ");
                Atom OH = (Atom) target.getAtomNode("OH");
                Atom HH = (Atom) target.getAtomNode("HH");
                // SPECIAL CASE: have to create map manualy.
                Bond b1 = CG.getBond(CB);
                Angle a1 = CG.getAngle(CB, CA);
                Torsion t1 = CG.getTorsion(CB, CA, N);
                Bond b2 = CD1.getBond(CG);
                Angle a2 = CD1.getAngle(CG, CB);
                Torsion t2 = CD1.getTorsion(CG, CB, CA);
                Bond b3 = HH.getBond(OH);
                Angle a3 = HH.getAngle(OH, CZ);
                Torsion t3 = HH.getTorsion(OH, CZ, CE2);
                BackBondedList bbl1 = new BackBondedList(b1,a1,t1);
                BackBondedList bbl2 = new BackBondedList(b2,a2,t2);
                BackBondedList bbl3 = new BackBondedList(b3,a3,t3);
                map.put(0, bbl1);
                map.put(1, bbl2);
                map.put(2, bbl3);
                return map;                     // Note the return here.
            }
            case TYD: {
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                Atom CG = (Atom) target.getAtomNode("CG");
                keyAtoms.add(CG);
                keyAtoms.add(CD1);
                break;
            }
            case TRP: {
                Atom CD1 = (Atom) target.getAtomNode("CD1");
                Atom CG = (Atom) target.getAtomNode("CG");
                keyAtoms.add(CG);
                keyAtoms.add(CD1);
                break;
            }
            case HIS:
            case HID: 
            case HIE: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom ND1 = (Atom) target.getAtomNode("ND1");
                keyAtoms.add(CG);
                keyAtoms.add(ND1);
                break;
            }
            case ASP: {
                Atom CG = (Atom) target.getAtomNode("CG");
                keyAtoms.add(CG);
                break;
            }
            case ASH: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom OD1 = (Atom) target.getAtomNode("OD1");
                keyAtoms.add(CG);
                keyAtoms.add(OD1);
                break;
            }
            case ASN: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom OD1 = (Atom) target.getAtomNode("OD1");
                keyAtoms.add(CG);
                keyAtoms.add(OD1);
                break;
            }
            case GLU: 
            case GLH: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CD = (Atom) target.getAtomNode("CD");
                Atom OE1 = (Atom) target.getAtomNode("OE1");
                keyAtoms.add(CG);
                keyAtoms.add(CD);
                keyAtoms.add(OE1);
                break;
            }
            case GLN: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CD = (Atom) target.getAtomNode("CD");
                Atom OE1 = (Atom) target.getAtomNode("OE1");
                keyAtoms.add(CG);
                keyAtoms.add(CD);
                keyAtoms.add(OE1);
                break;
            }
            case MET: {
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CE = (Atom) target.getAtomNode("CE");
                Atom SD = (Atom) target.getAtomNode("SD");
                keyAtoms.add(CG);
                keyAtoms.add(SD);
                keyAtoms.add(CE);
                break;
            }
            case LYS: 
            case LYD: {
                Atom CD = (Atom) target.getAtomNode("CD");
                Atom CE = (Atom) target.getAtomNode("CE");
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom NZ = (Atom) target.getAtomNode("NZ");
                keyAtoms.add(CG);
                keyAtoms.add(CD);
                keyAtoms.add(CE);
                keyAtoms.add(NZ);
                break;
            }
            case ARG: {
                Atom CD = (Atom) target.getAtomNode("CD");
                Atom CG = (Atom) target.getAtomNode("CG");
                Atom CZ = (Atom) target.getAtomNode("CZ");
                Atom NE = (Atom) target.getAtomNode("NE");
                keyAtoms.add(CG);
                keyAtoms.add(CD);
                keyAtoms.add(NE);
                keyAtoms.add(CZ);
                break;
            }
            default:
                logger.severe("CBMC called on unsupported residue.");
        }
        // Build the chain and assign back-bonded terms.
        chain.add(N);
        chain.add(CA);
        chain.add(CB);
        chain.addAll(keyAtoms);
        for (int i = 3; i < chain.size(); i++) {
            Atom key = chain.get(i);
            Bond bond = key.getBond(chain.get(i-1));
            Angle angle = key.getAngle(chain.get(i-1), chain.get(i-2));
            Torsion torsion = key.getTorsion(chain.get(i-1),
                    chain.get(i-2), chain.get(i-3));
            BackBondedList bbl = new BackBondedList(bond, angle, torsion);
            map.put(i-3, bbl);
        }
        return map;
    }
    
    private void writeSnapshot(String suffix, boolean append) {
        if (snapshotWriter != null) {
            snapshotWriter.write(suffix, append);
        }
    }
    
    private class SnapshotWriter {
        private final MolecularAssembly mola;
        private final PDBFilter filter;
        private final boolean interleaving;
        private SnapshotWriter(MolecularAssembly mola, boolean interleaving) {
            this.mola = mola;
            this.interleaving = interleaving;
            this.filter = new PDBFilter(mola.getFile(), mola, null, null);
        }
        
        private void write(String suffix, boolean append) {
            String filename = FilenameUtils.removeExtension(mola.getFile().toString()) + "." + suffix + "-" + moveNumber;
            if (interleaving) {
                filename = mola.getFile().getAbsolutePath();
                if (!filename.contains("dyn")) {
                    filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
                }
            }
            File file = new File(filename);
            filter.writeFile(file, append);
        }
    }
    
}
