/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.autoparm;

import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiTorsionType;
import ffx.potential.parameters.StretchBendType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

/**
 * Stitches PolType fragments together to form a whole
 * parent molecule force field.  
 * 
 * Stitching is done by averaging duplicate parameter values. 
 *
 * @author Rae Ann Corrigan
 */
public class Stitch {
    
    private ArrayList<File> sdfFiles = null;
    private ArrayList<File> uniqueAtomNamesTextFiles = null;
    
    private final static Logger logger = Logger.getLogger(Stitch.class.getName());

    //constructor
    public Stitch(ArrayList<File> sdfFiles, ArrayList<File> uniqueAtomNamesTextFiles) {
        this.sdfFiles = sdfFiles;
        this.uniqueAtomNamesTextFiles = uniqueAtomNamesTextFiles;
    }
    
    public ForceField combinePatches() {
        ForceField parent = new ForceField(null);
        
        int nFragments = sdfFiles.size();

        // Create array of forcefields
        ForceField[] forcefield = new ForceField[nFragments];

        // Read in parent atom names.
        String parentAtomNames[] = null;
        
        for (int i = 0; i < nFragments; i++) {
            // Loop over fragments.
            // Read in force field patch
            ForceField currentPatch = new ForceField(null);
            // Read in Rae's atom labels
            HashMap<String, Integer> fragmentMap = new HashMap<>();
            // Match atom labels to AtomType instances
            HashMap<Integer, Integer> typeMap = createMap(parentAtomNames, fragmentMap);
            // Overwrite AtomType name (i.e. PolType atom names are wrong)
            // Loop over all force field terms
            // If term includes only valid atoms (i.e. valid atom names) from the overall molecule, 
            // add it to parent force field

            // Angles
            AngleType fragAngleType = null;
            // Map
            int[] polTypeAngleClasses = fragAngleType.atomClasses;
            updateAtomClasses(polTypeAngleClasses, typeMap);
            fragAngleType.setKey(AngleType.sortKey(polTypeAngleClasses));
            String angleKey = fragAngleType.getKey();
            AngleType parentAngleType = parent.getAngleType(angleKey);
            if (parentAngleType == null) {
                parent.addForceFieldType(fragAngleType);
            } else {
                // Stitch (or average)
                AngleType averageAngleType = AngleType.average(fragAngleType, parentAngleType, polTypeAngleClasses);
                parent.addForceFieldType(averageAngleType);
            }

            // Bonds
            BondType fragBondType = null;
            // Map 
            int[] polTypeBondClasses = fragBondType.atomClasses;
            updateAtomClasses(polTypeBondClasses, typeMap);
            fragBondType.setKey(BondType.sortKey(polTypeBondClasses));
            String bondKey = fragBondType.getKey();
            BondType parentBondType = parent.getBondType(bondKey);
            if (parentBondType == null) {
                parent.addForceFieldType(fragBondType);
            } else {
                // Strich (or average)
                BondType averageBondType = BondType.average(fragBondType, parentBondType, polTypeBondClasses);
                parent.addForceFieldType(averageBondType);
            }

            // Multipoles
            MultipoleType fragMultipoleType = null;
            // Map
            int[] multipoleFrameTypes = fragMultipoleType.frameAtomTypes;
            updateAtomClasses(multipoleFrameTypes, typeMap);
            fragMultipoleType.setKey(multipoleFrameTypes);
            String multipoleKey = fragMultipoleType.getKey();
            MultipoleType parentMultipoleType = parent.getMultipoleType(multipoleKey);
            if (parentMultipoleType == null) {
                parent.addForceFieldType(fragMultipoleType);
            } else {
                // Stitch (or average)
                MultipoleType averageMultipoleType = MultipoleType.average(fragMultipoleType, parentMultipoleType, multipoleFrameTypes);
                parent.addForceFieldType(averageMultipoleType);
            }

            // Out of Plane Bends
            OutOfPlaneBendType fragOutOfPlaneBendType = null;
            // Map
            int[] polTypeOutOfPlaneBendClasses = fragOutOfPlaneBendType.atomClasses;
            updateAtomClasses(polTypeOutOfPlaneBendClasses, typeMap);
            fragOutOfPlaneBendType.setKey(polTypeOutOfPlaneBendClasses);
            String outOfPlaneBendKey = fragOutOfPlaneBendType.getKey();
            OutOfPlaneBendType parentOutOfPlaneBendType = parent.getOutOfPlaneBendType(outOfPlaneBendKey);
            if (parentOutOfPlaneBendType == null) {
                parent.addForceFieldType(fragOutOfPlaneBendType);
            } else {
                // Stitch (or average)
                OutOfPlaneBendType averageOutOfPlaneBendType = OutOfPlaneBendType.average(fragOutOfPlaneBendType, parentOutOfPlaneBendType, polTypeOutOfPlaneBendClasses);
                parent.addForceFieldType(averageOutOfPlaneBendType);
            }

            // Pi Torsions
            PiTorsionType fragPiTorsionType = null;
            // Map
            int[] polTypePiTorsionClasses = fragPiTorsionType.atomClasses;
            updateAtomClasses(polTypePiTorsionClasses, typeMap);
            fragPiTorsionType.setKey(PiTorsionType.sortKey(polTypePiTorsionClasses));
            String piTorsionKey = fragPiTorsionType.getKey();
            PiTorsionType parentPiTorsionType = parent.getPiTorsionType(piTorsionKey);
            if (parentPiTorsionType == null) {
                parent.addForceFieldType(fragPiTorsionType);
            } else {
                // Stitch (or average)
                PiTorsionType averagePiTorsionType = PiTorsionType.average(fragPiTorsionType, parentPiTorsionType, polTypePiTorsionClasses);
                parent.addForceFieldType(averagePiTorsionType);
            }

            // Stretch Bends
            StretchBendType fragStretchBendType = null;
            // Map
            int[] polTypeStretchBendClasses = fragStretchBendType.atomClasses;
            updateAtomClasses(polTypeStretchBendClasses, typeMap);
            fragStretchBendType.setKey(StretchBendType.sortKey(polTypeStretchBendClasses));
            String stretchBendKey = fragStretchBendType.getKey();
            StretchBendType parentStretchBendType = parent.getStretchBendType(stretchBendKey);
            if (parentStretchBendType == null) {
                parent.addForceFieldType(fragStretchBendType);
            } else {
                // Stitch (or average)
                StretchBendType averageStretchBendType = StretchBendType.average(fragStretchBendType, parentStretchBendType, polTypeStretchBendClasses);
                parent.addForceFieldType(averageStretchBendType);
            }

            // Torsion-Torsion
            TorsionTorsionType fragTorsionTorsionType = null;
            // Map
            int[] polTypeTorsionTorsionClasses = fragTorsionTorsionType.atomClasses;
            updateAtomClasses(polTypeTorsionTorsionClasses, typeMap);
            fragTorsionTorsionType.setKey(TorsionTorsionType.sortKey(polTypeTorsionTorsionClasses));
            String torsionTorsionKey = fragTorsionTorsionType.getKey();
            TorsionTorsionType parentTorsionTorsionType = parent.getTorsionTorsionType(torsionTorsionKey);
            if (parentTorsionTorsionType == null) {
                parent.addForceFieldType(fragTorsionTorsionType);
            } else {
                // Stitch (or average)
                TorsionTorsionType averageTorsionTorsionType = TorsionTorsionType.average(fragTorsionTorsionType, parentTorsionTorsionType, polTypeTorsionTorsionClasses);
                parent.addForceFieldType(averageTorsionTorsionType);
            }

            // Torsion
            TorsionType fragTorsionType = null;
            // Map
            int[] polTypeTorsionClasses = fragTorsionType.atomClasses;
            updateAtomClasses(polTypeTorsionClasses, typeMap);
            fragTorsionType.setKey(TorsionType.sortKey(polTypeTorsionClasses));
            String torsionKey = fragTorsionType.getKey();
            TorsionType parentTorsionType = parent.getTorsionType(torsionKey);
            if (parentTorsionType == null) {
                parent.addForceFieldType(fragTorsionType);
            } else {
                // Stitch (or average)
                TorsionType averageTorsionType = TorsionType.average(fragTorsionType, parentTorsionType, polTypeTorsionClasses);
                parent.addForceFieldType(averageTorsionType);
            }

            // Urey-Bradley
            UreyBradleyType fragUreyBradleyType = null;
            // Map
            int[] polTypeUreyBradleyClasses = fragUreyBradleyType.atomClasses;
            updateAtomClasses(polTypeUreyBradleyClasses, typeMap);
            fragUreyBradleyType.setKey(UreyBradleyType.sortKey(polTypeUreyBradleyClasses));
            String ureyBradleyKey = fragUreyBradleyType.getKey();
            UreyBradleyType parentUreyBradleyType = parent.getUreyBradleyType(ureyBradleyKey);
            if (parentUreyBradleyType == null) {
                parent.addForceFieldType(fragUreyBradleyType);
            } else {
                // Stitch (or average)
                UreyBradleyType averageUreyBradleyType = UreyBradleyType.average(fragUreyBradleyType, parentUreyBradleyType, polTypeUreyBradleyClasses);
                parent.addForceFieldType(averageUreyBradleyType);
            }

            // van der Waals
            VDWType fragVDWType = null;
            // Map
            int[] polTypeVDWClass = null;
            polTypeVDWClass[0] = fragVDWType.atomClass;
            updateAtomClasses(polTypeVDWClass, typeMap);
            fragVDWType.setKey(polTypeVDWClass);
            String vdwKey = fragVDWType.getKey();
            VDWType parentVDWType = parent.getVDWType(vdwKey);
            if (parentVDWType == null) {
                parent.addForceFieldType(fragVDWType);
            } else {
                // Stitch (or average)
                VDWType averageVDWType = VDWType.average(fragVDWType, parentVDWType, polTypeVDWClass[0]);
                parent.addForceFieldType(fragVDWType);
            }
            
        } // End "loop over fragments" for loop
        return parent;
    }
    
    private HashMap<Integer, Integer> createMap(String[] parentAtomNames, HashMap<String, Integer> fragmentMap) {
        HashMap<Integer, Integer> hashMap = new HashMap<>();
        int numParentAtoms = parentAtomNames.length;
        for (String key : fragmentMap.keySet()) {
            Integer type = fragmentMap.get(key);
            for (int i = 0; i < numParentAtoms; i++) {
                if (parentAtomNames[i].equalsIgnoreCase(key)) {
                    int newType = i + 1;
                    hashMap.put(type, newType);
                }
            }
        }
        return hashMap;
    }
    
    private void updateAtomClasses(int currentTypes[], HashMap<Integer, Integer> map) {
        if (currentTypes == null) {
            return;
        }
        for (int i = 0; i < currentTypes.length; i++) {
            currentTypes[i] = map.get(currentTypes[i]);
        }
    }
    
}
