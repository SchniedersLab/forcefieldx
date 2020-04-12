//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.nonbonded;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Logger;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

/**
 * Apply a restraint between groups.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainGroups {

    private static final Logger logger = Logger.getLogger(RestrainGroups.class.getName());

    /**
     * The MolecularAssembly to operate on.
     */
    final MolecularAssembly molecularAssembly;
    /**
     * Atom array.
     */
    final Atom[] atoms;
    /**
     * Number of atoms.
     */
    final int nAtoms;
    /**
     * Number of groups.
     */
    final int nGroups;
    /**
     * The atoms in each group.
     */
    final int[][] groupMembers;
    /**
     * The molecule number of each atom group (or -1 if the atoms are in different molecules
     */
    final int[] groupMolecule;
    /**
     * The mass of each group.
     */
    final double[] groupMass;
    /**
     * Number of group restraints.
     */
    final int nRestraints;
    /**
     * Index of group 1 for each restraint.
     */
    final int[] group1;
    /**
     * Index of group 2 for each restraint.
     */
    final int[] group2;
    /**
     * Force constant each restraint.
     */
    final double[] forceConstants;
    /**
     * Distance 1 for each restraint.
     */
    final double[] distance1;
    /**
     * Distance 1 for each restraint.
     */
    final double[] distance2;
    /**
     * Flag to indicate of all atoms for both groups are part of the same molecule.
     */
    final boolean[] sameMolecule;

    /**
     * Group Restraint Constructor.
     *
     * @param molecularAssembly The molecularAssembly to operate on.
     */
    public RestrainGroups(MolecularAssembly molecularAssembly) {
        this.molecularAssembly = molecularAssembly;
        this.atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;

        CompositeConfiguration compositeConfiguration = molecularAssembly.getProperties();

        String[] groups = compositeConfiguration.getStringArray("group");
        if (groups == null || groups.length < 2) {
            logger.severe(" restrain-groups requires two groups to be defined.");
        }

        String[] restrainGroups = compositeConfiguration.getStringArray("restrain-groups");
        if (restrainGroups == null || restrainGroups.length < 1) {
            logger.severe(" No restrain-groups property found.");
        }

        HashMap<Integer, ArrayList<Integer>> groupMap = new HashMap<>();
        // Loop over groups.
        for (String g : groups) {
            logger.fine(format(" Parsing group:\n %s.", g));
            // Split group lines on 1 or mores spaces.
            String[] gs = g.trim().split(" +");
            // First entry is the group number.
            int groupNumber = parseInt(gs[0]) - 1;
            ArrayList<Integer> groupAtomIDs = new ArrayList<>();
            // Loop over group ranges.
            for (int i = 1; i < gs.length; i++) {
                String[] range = gs[i].split(",");
                // Only one entry (not a range).
                if (range.length == 1) {
                    groupAtomIDs.add(parseInt(range[0]) - 1);
                } else if (range.length == 2) {
                    int start = -parseInt(range[0]) - 1;
                    int end = parseInt(range[1]) - 1;
                    if (start < 0 || start >= nAtoms || end < start || end >= nAtoms) {
                        logger.severe(format(" Property group could not be parsed:\n %s.", g));
                        continue;
                    }
                    for (int j = start; j <= end; j++) {
                        groupAtomIDs.add(j);
                    }
                } else {
                    logger.severe(format(" Property group could not be parsed:\n %s.", g));
                    continue;
                }
                groupMap.put(groupNumber, groupAtomIDs);
            }
        }

        // Pack groups into arrays.
        int[] molecule = molecularAssembly.getMoleculeNumbers();
        nGroups = groupMap.size();
        groupMembers = new int[nGroups][];
        groupMass = new double[nGroups];
        groupMolecule = new int[nGroups];
        sameMolecule = new boolean[nGroups];
        for (Integer groupNumber : groupMap.keySet()) {
            ArrayList<Integer> members = groupMap.get(groupNumber);
            if (groupNumber >= nGroups) {
                logger.severe(" Please label groups from 1 to the number of groups.");
            }
            int groupID = groupNumber;
            groupMembers[groupID] = new int[members.size()];
            int mem = 0;
            for (int member : members) {
                groupMembers[groupID][mem++] = member;
            }
            logger.fine(format(" Group %d members %s.", groupID, Arrays.toString(groupMembers[groupID])));
            int k = groupMembers[groupID][0];
            Atom atom = atoms[k];
            groupMass[groupID] = atom.getMass();
            groupMolecule[groupID] = molecule[k];
            for (int i = 1; i < groupMembers[groupID].length; i++) {
                k = groupMembers[groupID][i];
                atom = atoms[k];
                groupMass[groupID] += atom.getMass();
                if (groupMolecule[groupID] != molecule[k]) {
                    groupMolecule[groupID] = -1;
                }
            }
            groupMass[groupID] = max(groupMass[groupID], 1.0);
        }

        // Parse restrain-groups properties.
        nRestraints = restrainGroups.length;
        group1 = new int[nRestraints];
        group2 = new int[nRestraints];
        forceConstants = new double[nRestraints];
        distance1 = new double[nRestraints];
        distance2 = new double[nRestraints];
        int iRestraint = 0;
        for (String restraint : restrainGroups) {
            logger.fine(format(" Parsing restrain-groups:\n %s.", restraint));
            String[] values = restraint.trim().split(" +");
            if (values.length < 3) {
                logger.severe(format(" Property restrain-groups could not be parsed:\n %s.", restraint));
            }
            group1[iRestraint] = parseInt(values[0]) - 1;
            group2[iRestraint] = parseInt(values[1]) - 1;
            int i1 = group1[iRestraint];
            int i2 = group2[iRestraint];
            if (i1 >= nGroups || i2 >= nGroups) {
                logger.severe(format(" Property restrain-groups has invalid groups:\n %s.", restraint));
            }
            forceConstants[iRestraint] = parseDouble(values[2]);
            distance1[iRestraint] = 0.0;
            distance2[iRestraint] = 0.0;
            if (values.length > 3) {
                distance1[iRestraint] = parseDouble(values[3]);
                if (values.length > 4) {
                    distance2[iRestraint] = parseDouble(values[4]);
                }
            }
            sameMolecule[iRestraint] = false;
            if (groupMolecule[i1] > -1 && groupMolecule[i2] > -1) {
                if (groupMolecule[i1] == groupMolecule[i2]) {
                    sameMolecule[iRestraint] = true;
                }
            }
        }
    }

    /**
     * Compute energy and derivatives for group distance restraint terms.
     *
     * @param gradient If true, compute the derivative.
     * @return The potential energy.
     */
    public double energy(boolean gradient) {
        double energy = 0;
        double[] xyz = new double[3];
        double[] dr = new double[3];
        Crystal crystal = molecularAssembly.getCrystal();
        for (int i = 0; i < nRestraints; i++) {
            int i1 = group1[i];
            int i2 = group2[i];
            // Loop over the first group members.
            double xcm = 0.0;
            double ycm = 0.0;
            double zcm = 0.0;
            for (int j = 0; j < groupMembers[i1].length; j++) {
                int k = groupMembers[i1][j];
                Atom atom = atoms[k];
                double mass = atom.getMass();
                atom.getXYZ(xyz);
                xcm = xcm + xyz[0] * mass;
                ycm = ycm + xyz[1] * mass;
                zcm = zcm + xyz[2] * mass;
            }
            double mass1 = groupMass[i1];
            double xr = xcm / mass1;
            double yr = ycm / mass1;
            double zr = zcm / mass1;
            // Loop over the second group members.
            xcm = 0.0;
            ycm = 0.0;
            zcm = 0.0;
            for (int j = 0; j < groupMembers[i2].length; j++) {
                int k = groupMembers[i2][j];
                Atom atom = atoms[k];
                double mass = atom.getMass();
                atom.getXYZ(xyz);
                xcm = xcm + xyz[0] * mass;
                ycm = ycm + xyz[1] * mass;
                zcm = zcm + xyz[2] * mass;
            }
            double mass2 = groupMass[i2];
            xr = xr - xcm / mass2;
            yr = yr - ycm / mass2;
            zr = zr - zcm / mass2;
            double r2;
            if (sameMolecule[i]) {
                r2 = xr * xr + yr * yr + zr * zr;
            } else {
                dr[0] = xr;
                dr[1] = yr;
                dr[2] = zr;
                r2 = crystal.image(dr);
            }
            double r = sqrt(r2);
            double force = forceConstants[i];
            double gf1 = distance1[i];
            double gf2 = distance2[i];
            double target = r;
            if (r < gf1) target = gf1;
            if (r > gf2) target = gf2;
            double dt = r - target;
            double dt2 = dt * dt;
            double e = force * dt2;
            energy = energy + e;
            if (gradient) {
                // Compute chain rule terms needed for derivatives.
                if (r == 0.0) r = 1.0;
                double de = 2.0 * force * dt / r;
                double dedx = de * xr;
                double dedy = de * yr;
                double dedz = de * zr;
                // Increment the total energy and first derivatives for group 1.
                for (int j = 0; j < groupMembers[i1].length; j++) {
                    int k = groupMembers[i1][j];
                    Atom atom = atoms[k];
                    double mass = atom.getMass();
                    double ratio = mass / mass1;
                    atom.addToXYZGradient(dedx * ratio, dedy * ratio, dedz * ratio);
                }
                // Increment the total energy and first derivatives for group 2.
                for (int j = 0; j < groupMembers[i2].length; j++) {
                    int k = groupMembers[i2][j];
                    Atom atom = atoms[k];
                    double mass = atom.getMass();
                    double ratio = mass / mass2;
                    atom.addToXYZGradient(-dedx * ratio, -dedy * ratio, -dedz * ratio);
                }
            }
        }
        return energy;
    }

    /**
     * Get the number of group restraints.
     *
     * @return The number of group restraints.
     */
    public int getNumberOfRestraints() {
        return nRestraints;
    }

}
