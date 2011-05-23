/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ffx.autoparm;

/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */

import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static ffx.numerics.VectorMath.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.Potential;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.bonded.Utilities;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.autoparm.PME_2;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parsers.XYZFilter;
import ffx.utilities.Keyword;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.math.linear.*;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Energy implements Potential {

    private static final Logger logger = Logger.getLogger(Energy.class.getName());
    private final Atom[] atoms;
    private final Crystal crystal;
    private final ParallelTeam parallelTeam;
    private final Bond bonds[];
    private final Angle angles[];
    private final StretchBend stretchBends[];
    private final UreyBradley ureyBradleys[];
    private final OutOfPlaneBend outOfPlaneBends[];
    private final Torsion torsions[];
    private final PiOrbitalTorsion piOrbitalTorsions[];
    private final TorsionTorsion torsionTorsions[];
    private final VanDerWaals vanderWaals;
    //private final CellCellVanDerWaals vanderWaals;
    private final ParticleMeshEwald particleMeshEwald;
    private PME_2 pme2;
    protected final int nAtoms;
    protected final int nBonds;
    protected final int nAngles;
    protected final int nStretchBends;
    protected final int nUreyBradleys;
    protected final int nOutOfPlaneBends;
    protected final int nTorsions;
    protected final int nPiOrbitalTorsions;
    protected final int nTorsionTorsions;
    protected int nVanDerWaals, nPME;
    protected final boolean bondTerm;
    protected final boolean angleTerm;
    protected final boolean stretchBendTerm;
    protected final boolean ureyBradleyTerm;
    protected final boolean outOfPlaneBendTerm;
    protected final boolean torsionTerm;
    protected final boolean piOrbitalTorsionTerm;
    protected final boolean torsionTorsionTerm;
    protected final boolean vanDerWaalsTerm;
    protected final boolean multipoleTerm;
    protected final boolean polarizationTerm;
    protected double bondEnergy, bondRMSD;
    protected double angleEnergy, angleRMSD;
    protected double stretchBendEnergy;
    protected double ureyBradleyEnergy;
    protected double outOfPlaneBendEnergy;
    protected double torsionEnergy;
    protected double piOrbitalTorsionEnergy;
    protected double torsionTorsionEnergy;
    protected double totalBondedEnergy;
    protected double vanDerWaalsEnergy;
    protected double permanentMultipoleEnergy;
    protected double polarizationEnergy;
    protected double totalElectrostaticEnergy;
    protected double totalNonBondedEnergy;
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long totalTime;
    protected double[] optimizationScaling = null;
    private static final double toSeconds = 0.000000001;
    private File structure_key;
    private File structure_xyz;
    InputStreamReader stdinput = new InputStreamReader(System.in);
    BufferedReader stdreader = new BufferedReader(stdinput);
    private MolecularAssembly molecularAssembly;
    private ArrayList<String> key = new ArrayList<String>();
    private ForceField forceField;

    public Energy(String xyz_filename) throws IOException{
        
        parallelTeam = new ParallelTeam();
        logger.info(format(" Constructing Force Field"));
        logger.info(format("\n SMP threads:                        %10d", parallelTeam.getThreadCount()));

        structure_xyz = new File(xyz_filename);
        if (!(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead())) {
            System.out.println("Couldn't find xyz file");
            System.exit(1);
        }
        int n = 1;
        String oxyzfname = null;
        String old = xyz_filename;
        while (structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()) {
            oxyzfname = xyz_filename;
            n++;
            xyz_filename = old;
            xyz_filename = xyz_filename + "_" + Integer.toString(n);
            structure_xyz = new File(xyz_filename);
        }
        structure_xyz = new File(oxyzfname);
        int index = oxyzfname.lastIndexOf(".");
        String name = oxyzfname.substring(0, index);
        String keyfname = name + ".key";
        structure_key = new File(keyfname);
        while (!(structure_key != null && structure_key.exists() && structure_key.canRead())) {
            System.out.println("Enter the Key Filename: ");
            keyfname = stdreader.readLine();
            structure_key = new File(keyfname);
            if (!(structure_key != null && structure_key.exists() && structure_key.canRead())) {
                System.out.println("Couldn't find key file");
            }
        }

        n = 1;
        String okeyfname = null;
        old = keyfname;
        while (structure_key != null && structure_key.exists() && structure_key.canRead() && structure_key.length() != 0) {
            okeyfname = keyfname;
            n++;
            keyfname = old;
            keyfname = keyfname + "_" + Integer.toString(n);
            structure_key = new File(keyfname);
        }

        structure_key = new File(okeyfname);
        
        molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure_xyz);
        CompositeConfiguration properties = Keyword.loadProperties(structure_key);
        ForceFieldFilter_2 forceFieldFilter = new ForceFieldFilter_2(properties, structure_key);
        forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);
        XYZFilter xyzFilter = new XYZFilter(structure_xyz, molecularAssembly, forceField, properties);
        xyzFilter.readFile();
        Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
        molecularAssembly.finalize(true);

        
        // Get a reference to the sorted atom array.
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;

        //ForceField forceField = molecularAssembly.getForceField();
        bondTerm = forceField.getBoolean(ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(ForceFieldBoolean.TORSIONTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(ForceFieldBoolean.TORTORTERM, true);
        vanDerWaalsTerm = forceField.getBoolean(ForceFieldBoolean.VDWTERM, true);
        multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
        polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);

        // Define the cutoff lengths.
        double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        double buff = 2.0;
        double cutOff2 = 2.0 * (max(vdwOff, ewaldOff) + buff);

        // Determine the unit cell dimensions and Spacegroup
        String spacegroup;
        double a, b, c, alpha, beta, gamma;
        boolean aperiodic;
        try {
            a = forceField.getDouble(ForceFieldDouble.A_AXIS);
            aperiodic = false;
            b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
            c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
            alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
            beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
            gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
            spacegroup = forceField.getString(ForceFieldString.SPACEGROUP, "P1");
        } catch (Exception e) {
            logger.info(" The system will be treated as aperiodic.");
            aperiodic = true;
            spacegroup = "P1";
            /**
             * Search all atom pairs to find the largest pair-wise distance.
             */
            double xr[] = new double[3];
            double maxr = 0.0;
            for (int i = 0; i < nAtoms - 1; i++) {
                double[] xi = atoms[i].getXYZ();
                for (int j = i + 1; j < nAtoms; j++) {
                    double[] xj = atoms[j].getXYZ();
                    diff(xi, xj, xr);
                    double r = r(xr);
                    if (r > maxr) {
                        maxr = r;
                    }
                }
            }
            a = 2.0 * (maxr + ewaldOff);
            b = a;
            c = a;
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }
        Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        unitCell.setAperiodic(aperiodic);

        /**
         * Do we need a ReplicatesCrystal?
         */
        int l = 1;
        int m = 1;
        n = 1;
        while (unitCell.a * l < cutOff2) {
            l++;
        }
        while (unitCell.b * m < cutOff2) {
            m++;
        }
        while (unitCell.c * n < cutOff2) {
            n++;
        }

        if (l * m * n > 1 && !aperiodic) {
            this.crystal = new ReplicatesCrystal(unitCell, l, m, n);
        } else {
            this.crystal = unitCell;
        }

        logger.info(crystal.toString());

        boolean rigidHydrogens = forceField.getBoolean(ForceFieldBoolean.RIGID_HYDROGENS, false);
        double rigidScale = forceField.getDouble(ForceFieldDouble.RIGID_SCALE, 10.0);

        if (rigidScale <= 1.0) {
            rigidScale = 1.0;
        }

        logger.info("\n Bonded Terms\n");
        if (rigidHydrogens && rigidScale > 1.0) {
            logger.info(format(" Rigid hydrogens:                    %10.2f", rigidScale));
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
            logger.info(format(" Bonds:                              %10d", nBonds));
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<ROLS> angle = molecularAssembly.getAngleList();
            nAngles = angle.size();
            angles = angle.toArray(new Angle[nAngles]);
            Arrays.sort(angles);
            logger.info(format(" Angles:                             %10d", nAngles));
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<ROLS> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = stretchBend.size();
            stretchBends = stretchBend.toArray(new StretchBend[nStretchBends]);
            Arrays.sort(stretchBends);
            logger.info(format(" Stretch-Bends:                      %10d", nStretchBends));
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<ROLS> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = ureyBradley.size();
            ureyBradleys = ureyBradley.toArray(new UreyBradley[nUreyBradleys]);
            Arrays.sort(ureyBradleys);
            logger.info(format(" Urey-Bradleys:                      %10d", nUreyBradleys));
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }

        /**
         * Set a multiplier on the force constants of bonded terms containing
         * hydrogens.
         */
        if (rigidHydrogens) {
            if (bonds != null) {
                for (Bond bond : bonds) {
                    if (bond.containsHydrogen()) {
                        bond.setRigidScale(rigidScale);
                    }
                }
            }
            if (angles != null) {
                for (Angle angle : angles) {
                    if (angle.containsHydrogen()) {
                        angle.setRigidScale(rigidScale);
                    }
                }
            }
            if (stretchBends != null) {
                for (StretchBend stretchBend : stretchBends) {
                    if (stretchBend.containsHydrogen()) {
                        stretchBend.setRigidScale(rigidScale);
                    }
                }
            }
            if (ureyBradleys != null) {
                for (UreyBradley ureyBradley : ureyBradleys) {
                    if (ureyBradley.containsHydrogen()) {
                        ureyBradley.setRigidScale(rigidScale);
                    }
                }
            }
        }

        // Collect, count, pack and sort out-of-plane bends.
        if (outOfPlaneBendTerm) {
            ArrayList<ROLS> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = outOfPlaneBend.size();
            outOfPlaneBends = outOfPlaneBend.toArray(new OutOfPlaneBend[nOutOfPlaneBends]);
            Arrays.sort(outOfPlaneBends);
            logger.info(format(" Out-of-Plane Bends:                 %10d", nOutOfPlaneBends));
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
            logger.info(format(" Torsions:                           %10d", nTorsions));
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
            logger.info(format(" Pi-Orbital Torsions:                %10d", nPiOrbitalTorsions));
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
            logger.info(format(" Torsion-Torsions:                   %10d", nTorsionTorsions));
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        if (vanDerWaalsTerm) {
            //vanderWaals = new CellCellVanDerWaals(forceField, atoms, crystal, parallelTeam);
            vanderWaals = new VanDerWaals(forceField, atoms, crystal, parallelTeam);
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            particleMeshEwald = new ParticleMeshEwald(forceField, atoms, crystal, parallelTeam,
                                                      vanderWaals.getNeighborLists());
            pme2 = new PME_2(forceField, atoms, crystal, parallelTeam,
                                                      vanderWaals.getNeighborLists(),key);
            pme2.propyze = true;
            pme2.init_prms();
        } else {
            particleMeshEwald = null;
        }

        //molecularAssembly.setPotential(this);
    }

    public double energy(boolean gradient, boolean print) {
        bondTime = 0;
        angleTime = 0;
        stretchBendTime = 0;
        ureyBradleyTime = 0;
        outOfPlaneBendTime = 0;
        torsionTime = 0;
        piOrbitalTorsionTime = 0;
        torsionTorsionTime = 0;
        vanDerWaalsTime = 0;
        electrostaticTime = 0;
        totalTime = System.nanoTime();

        // Zero out the potential energy of each bonded term.
        bondEnergy = 0.0;
        angleEnergy = 0.0;
        stretchBendEnergy = 0.0;
        ureyBradleyEnergy = 0.0;
        outOfPlaneBendEnergy = 0.0;
        torsionEnergy = 0.0;
        piOrbitalTorsionEnergy = 0.0;
        torsionTorsionEnergy = 0.0;
        totalBondedEnergy = 0.0;

        // Zero out bond and angle RMSDs.
        bondRMSD = 0.0;
        angleRMSD = 0.0;

        // Zero out the potential energy of each non-bonded term.
        vanDerWaalsEnergy = 0.0;
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        totalElectrostaticEnergy = 0.0;
        totalNonBondedEnergy = 0.0;

        // Zero out the total potential energy.
        totalEnergy = 0.0;

        // Zero out the Cartesian coordinate gradient for each atom.
        if (gradient) {
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.setXYZGradient(0.0, 0.0, 0.0);
            }
        }

        if (bondTerm) {
            bondTime = System.nanoTime();
            for (int i = 0; i < nBonds; i++) {
                Bond b = bonds[i];
                bondEnergy += b.energy(gradient);
                double value = b.getValue();
                bondRMSD += value * value;
            }
            bondRMSD = sqrt(bondRMSD / bonds.length);
            bondTime = System.nanoTime() - bondTime;
        }


        if (angleTerm) {
            angleTime = System.nanoTime();
            for (int i = 0; i < nAngles; i++) {
                Angle a = angles[i];
                angleEnergy += a.energy(gradient);
                double value = a.getValue();
                angleRMSD += value * value;
            }
            angleRMSD = sqrt(angleRMSD / angles.length);
            angleTime = System.nanoTime() - angleTime;
        }

        if (stretchBendTerm) {
            stretchBendTime = System.nanoTime();
            for (int i = 0; i < nStretchBends; i++) {
                stretchBendEnergy += stretchBends[i].energy(gradient);
            }
            stretchBendTime = System.nanoTime() - stretchBendTime;
        }

        if (ureyBradleyTerm) {
            ureyBradleyTime = System.nanoTime();
            for (int i = 0; i < nUreyBradleys; i++) {
                ureyBradleyEnergy += ureyBradleys[i].energy(gradient);
            }
            ureyBradleyTime = System.nanoTime() - ureyBradleyTime;
        }

        if (outOfPlaneBendTerm) {
            outOfPlaneBendTime = System.nanoTime();
            for (int i = 0; i < nOutOfPlaneBends; i++) {
                outOfPlaneBendEnergy += outOfPlaneBends[i].energy(gradient);
            }
            outOfPlaneBendTime = System.nanoTime() - outOfPlaneBendTime;
        }

        if (torsionTerm) {
            torsionTime = System.nanoTime();
            for (int i = 0; i < nTorsions; i++) {
                torsionEnergy += torsions[i].energy(gradient);
            }
            torsionTime = System.nanoTime() - torsionTime;
        }

        if (piOrbitalTorsionTerm) {
            piOrbitalTorsionTime = System.nanoTime();
            for (int i = 0; i < nPiOrbitalTorsions; i++) {
                piOrbitalTorsionEnergy += piOrbitalTorsions[i].energy(gradient);
            }
            piOrbitalTorsionTime = System.nanoTime() - piOrbitalTorsionTime;
            torsionTorsionTime = System.nanoTime();
        }

        if (torsionTorsionTerm) {
            for (int i = 0; i < nTorsionTorsions; i++) {
                torsionTorsionEnergy += torsionTorsions[i].energy(gradient);
            }
            torsionTorsionTime = System.nanoTime() - torsionTorsionTime;
        }

        if (vanDerWaalsTerm) {
            vanDerWaalsTime = System.nanoTime();
            vanDerWaalsEnergy = vanderWaals.energy(gradient, print);
            nVanDerWaals = this.vanderWaals.getInteractions();
            vanDerWaalsTime = System.nanoTime() - vanDerWaalsTime;
        }

        if (multipoleTerm) {
            electrostaticTime = System.nanoTime();
            totalElectrostaticEnergy = particleMeshEwald.energy(gradient, print);
            permanentMultipoleEnergy = particleMeshEwald.getPermanentEnergy();
            polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
            nPME = particleMeshEwald.getInteractions();
            electrostaticTime = System.nanoTime() - electrostaticTime;
        }
        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy;
        if (print) {
            StringBuilder sb = new StringBuilder("\n");
            if (gradient) {
                sb.append(" Computed Potential Energy and Atomic Coordinate Gradients\n");
            } else {
                sb.append(" Computed Potential Energy\n");
            }
            sb.append(this);
            logger.info(sb.toString());
        }
        system_mpoles();
        return totalEnergy;
    }

    public double getTotal() {
        return totalEnergy;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n");
        if (bondTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Bond Streching    ", bondEnergy, bonds.length,
                                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Angle Bending     ", angleEnergy, angles.length,
                                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Stretch-Bend      ", stretchBendEnergy,
                                    stretchBends.length, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Urey-Bradley      ", ureyBradleyEnergy,
                                    ureyBradleys.length, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                                    outOfPlaneBends.length, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsional Angle   ", torsionEnergy, torsions.length,
                                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                                    piOrbitalTorsions.length, piOrbitalTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsion-Torsion   ", torsionTorsionEnergy,
                                    torsionTorsions.length, torsionTorsionTime * toSeconds));
        }
        if (vanDerWaalsTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Van der Waals     ", vanDerWaalsEnergy,
                                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                                    "Atomic Multipoles ", permanentMultipoleEnergy, nPME));
        }
        if (polarizationTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Polarization      ", polarizationEnergy,
                                    nPME, electrostaticTime * toSeconds));
        }
        sb.append(String.format("\n %s %16.8f  %s %12.3f (sec)\n",
                                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));
        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format(" %s %16.8f\n", "Unit Cell         ",
                                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format(" %s %16.8f\n", "Replicates Cell   ",
                                    totalEnergy * nsymm));
        }

        return sb.toString();
    }

    public Crystal getCrystal() {
        return crystal;
    }

    public void setSoftCoreLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            vanderWaals.setLambda(lambda);
        } else {
            String message = String.format("Softcore lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    public void setElectrostaticsLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            particleMeshEwald.setLambda(lambda);
        } else {
            String message = String.format("Electrostatics lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    @Override
    public void setScaling(double scaling[]) {
        if (scaling != null) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    @Override
    public double energyAndGradient(double x[], double g[]) {
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }
        setCoordinates(x);
        double e = energy(true, false);
        getGradients(g);
        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }
        return e;
    }

    public void getGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZGradient(grad);
            double gx = grad[0];
            double gy = grad[1];
            double gz = grad[2];
            if (gx == Double.NaN || gx == Double.NEGATIVE_INFINITY || gx == Double.POSITIVE_INFINITY
                || gy == Double.NaN || gy == Double.NEGATIVE_INFINITY || gy == Double.POSITIVE_INFINITY
                || gz == Double.NaN || gz == Double.NEGATIVE_INFINITY || gz == Double.POSITIVE_INFINITY) {
                String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                                        a.toString(), gx, gy, gz);
                logger.warning(message);
            }
            g[index++] = gx;
            g[index++] = gy;
            g[index++] = gz;
        }
    }

    private void setCoordinates(double coords[]) {
        assert (coords != null);
        int index = 0;
        for (Atom a : atoms) {
            double x = coords[index++];
            double y = coords[index++];
            double z = coords[index++];
            a.moveTo(x, y, z);
        }
    }

    @Override
    public double[] getCoordinates(double x[]) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZ(xyz);
            x[index++] = xyz[0];
            x[index++] = xyz[1];
            x[index++] = xyz[2];
        }
        return x;
    }

    @Override
    public double[] getMass() {
        int n = getNumberOfVariables();
        double mass[] = new double[n];
        int i = 0;
        for (Atom a : atoms) {
            double m = a.getMass();
            mass[i++] = m;
            mass[i++] = m;
            mass[i++] = m;
        }
        return mass;
    }

    @Override
    public int getNumberOfVariables() {
        return nAtoms * 3;
    }
    
    public void system_mpoles(){
        //Find center of mass.
        double weigh = 0;
        double xyzmid[] = {0,0,0};
        double xyzcm[][] = new double[nAtoms][3];
        for(int i = 0; i < nAtoms; i++){
            weigh = weigh + atoms[i].getMass();
            for (int j = 0; j < 3; j++){
                xyzmid[j] = xyzmid[j] + atoms[i].getXYZ()[j] * atoms[i].getMass();
            }
        }
        if(weigh != 0){
            for(int j = 0; j<3; j++){
                xyzmid[j] = xyzmid[j]/weigh;
            }
        }
        
        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                xyzcm[i][j] = atoms[i].getXYZ()[j] - xyzmid[j];
            }
        }
        addInducedToGlobal();
        double netchg = 0, xdpl = 0, ydpl = 0, zdpl = 0, xxqdp = 0, xyqdp = 0, xzqdp = 0, yxqdp = 0, yyqdp = 0, yzqdp = 0, zxqdp = 0, zyqdp = 0, zzqdp = 0;
        for(int i = 0; i < nAtoms; i++){
            double charge = atoms[i].getMultipoleType().charge;
            double[] dipole = {pme2.globalMultipole[0][i][1], pme2.globalMultipole[0][i][2], pme2.globalMultipole[0][i][3]};
            //double[] dipole = atoms[i].getMultipoleType().dipole;
            netchg = netchg + charge;
            xdpl = xdpl + xyzcm[i][0] * charge + dipole[0];
            ydpl = ydpl + xyzcm[i][1] * charge + dipole[1];
            zdpl = zdpl + xyzcm[i][2] * charge + dipole[2];
            xxqdp = xxqdp + xyzcm[i][0] * xyzcm[i][0] * charge + 2 * xyzcm[i][0] * dipole[0];
            xyqdp = xyqdp + xyzcm[i][0] * xyzcm[i][1] * charge + xyzcm[i][0] * dipole[1] + xyzcm[i][1] * dipole[0];
            xzqdp = xzqdp + xyzcm[i][0] * xyzcm[i][2] * charge + xyzcm[i][0] * dipole[2] + xyzcm[i][2] * dipole[0];
            
            yxqdp = yxqdp + xyzcm[i][0] * xyzcm[i][1] * charge + xyzcm[i][0] * dipole[1] + xyzcm[i][1] * dipole[0];
            yyqdp = yyqdp + xyzcm[i][1] * xyzcm[i][1] * charge + 2 * xyzcm[i][1] * dipole[1];
            yzqdp = yzqdp + xyzcm[i][1] * xyzcm[i][2] * charge + xyzcm[i][1]*dipole[2] + xyzcm[i][2] * dipole[1];
            
            //zxqdp = zxqdp + xyzcm[i][2] * xyzcm[i][0] * charge + xyzcm[i][2] * dipole[0] + xyzcm[i][0] * dipole[2];
            zxqdp = zxqdp + xyzcm[i][0] * xyzcm[i][2] * charge + xyzcm[i][0] * dipole[2] + xyzcm[i][2] * dipole[0];
            //zyqdp = zyqdp + xyzcm[i][2] * xyzcm[i][1] * charge + xyzcm[i][2] * dipole[1] + xyzcm[i][1] * dipole[2];
            zyqdp = zyqdp + xyzcm[i][1] * xyzcm[i][2] * charge + xyzcm[i][1]*dipole[2] + xyzcm[i][2] * dipole[1];
            zzqdp = zzqdp + xyzcm[i][2] * xyzcm[i][2] * charge + 2 * xyzcm[i][2] * dipole[2];
        }
        

        
        double qave = (xxqdp + yyqdp + zzqdp)/3;
        xxqdp = 1.5 * (xxqdp - qave);
        xyqdp = 1.5 * xyqdp;
        xzqdp = 1.5 * xzqdp;
        yxqdp = 1.5 * yxqdp;
        yyqdp = 1.5 * (yyqdp - qave);
        yzqdp = 1.5 * yzqdp;
        zxqdp = 1.5 * zxqdp;
        zyqdp = 1.5 * zyqdp;
        zzqdp = 1.5 * (zzqdp - qave);
        
        for(int i = 0; i < nAtoms; i++){
            double[][] quadrupole = {{pme2.globalMultipole[0][i][4], pme2.globalMultipole[0][i][7], pme2.globalMultipole[0][i][8]},{pme2.globalMultipole[0][i][7], pme2.globalMultipole[0][i][5], pme2.globalMultipole[0][i][9]},{pme2.globalMultipole[0][i][8], pme2.globalMultipole[0][i][9], pme2.globalMultipole[0][i][6]}};
            //double[][] quadrupole = atoms[i].getMultipoleType().quadrupole;
            xxqdp = xxqdp + 3 * quadrupole[0][0];
            xyqdp = xyqdp + 3 * quadrupole[0][1];
            xzqdp = xzqdp + 3 * quadrupole[0][2];
            yxqdp = yxqdp + 3 * quadrupole[1][0];
            yyqdp = yyqdp + 3 * quadrupole[1][1];
            yzqdp = yzqdp + 3 * quadrupole[1][2];
            zxqdp = zxqdp + 3 * quadrupole[2][0];
            zyqdp = zyqdp + 3 * quadrupole[2][1];
            zzqdp = zzqdp + 3 * quadrupole[2][2];
        }
        
        xdpl = MultipoleType.DEBYE * xdpl;
        ydpl = MultipoleType.DEBYE * ydpl;
        zdpl = MultipoleType.DEBYE * zdpl;

        xxqdp = MultipoleType.DEBYE * xxqdp;
        xyqdp = MultipoleType.DEBYE * xyqdp;
        xzqdp = MultipoleType.DEBYE * xzqdp;
        yxqdp = MultipoleType.DEBYE * yxqdp;
        yyqdp = MultipoleType.DEBYE * yyqdp;
        yzqdp = MultipoleType.DEBYE * yzqdp;
        zxqdp = MultipoleType.DEBYE * zxqdp;
        zyqdp = MultipoleType.DEBYE * zyqdp;
        zzqdp = MultipoleType.DEBYE * zzqdp;
        
        double netdpl = Math.sqrt(xdpl*xdpl + ydpl*ydpl + zdpl*zdpl);
        
        RealMatrix a = new Array2DRowRealMatrix(new double[][] {{xxqdp,xyqdp,xzqdp},{yxqdp,yyqdp,yzqdp},{zxqdp,zyqdp,zzqdp}});

        EigenDecompositionImpl e = new EigenDecompositionImpl(a,1);
        a = e.getD();
        double[] netqdp = {a.getColumn(0)[0],a.getColumn(1)[1],a.getColumn(2)[2]};
        
        DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");
        String output;
        output = String.format(" Total Electric Charge:   %13s %s Electrons\n"," ",myFormatter.format(netchg));
        System.out.println(output);
        output = String.format(" Dipole Moment Magnitude: %13s %s Debyes\n"," ",myFormatter.format(netdpl));
        System.out.println(output);
        output = String.format(" Dipole X,Y,Z-Components: %13s %s %s %s\n"," ",myFormatter.format(xdpl),myFormatter.format(ydpl),myFormatter.format(zdpl));
        System.out.println(output);
        output = String.format(" Quadrupole Moment Tensor:%13s %s %s %s"," ",myFormatter.format(xxqdp),myFormatter.format(xyqdp),myFormatter.format(xzqdp));
        System.out.println(output);
        output = String.format("      (Buckinghams)       %13s %s %s %s"," ",myFormatter.format(yxqdp),myFormatter.format(yyqdp),myFormatter.format(yzqdp));
        System.out.println(output);
        output = String.format("                          %13s %s %s %s\n"," ",myFormatter.format(zxqdp),myFormatter.format(zyqdp),myFormatter.format(zzqdp));
        System.out.println(output);
        output = String.format("Principle Axes Quadrupole:%13s %s %s %s"," ",myFormatter.format(netqdp[2]),myFormatter.format(netqdp[1]),myFormatter.format(netqdp[0]));
        System.out.println(output);

    }
    
    public void addInducedToGlobal(){
        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3;j++){
                pme2.globalMultipole[0][i][j+1] = pme2.globalMultipole[0][i][j+1] + pme2.inducedDipole[0][i][j];
            }
        }
    }
    
//    public static void main(String args[]) throws IOException{
//        Energy e = new Energy("/users/gchattree/Research/Compounds/s_test4_compounds/triclosan/triclosan.xyz");
//        e.energy(false,true);
//    }
}

