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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;

/**
 * Stitches PolType fragments together to form a whole parent molecule force
 * field.
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

        // Create a CompositeConfiguration from the parent molecule file to pass to the ForceField consturctor
        CompositeConfiguration compositeConfiguration = new CompositeConfiguration();
        
        System.out.println("SDF files ArrayList size: "+sdfFiles.size());
        System.out.println("Text files ArrayList size: "+uniqueAtomNamesTextFiles.size());

        try {
            compositeConfiguration.addConfiguration(new PropertiesConfiguration(sdfFiles.get(0)));
        } catch (ConfigurationException e) {
            e.printStackTrace();
        }

        // This will be the final output forcefield for the full molecule
        ForceField parent = new ForceField(compositeConfiguration);

        // Minus one accounts for the parent SDF 
        // If using parent CIF for atom names, take out the minus 1
        int nFragments = sdfFiles.size() - 1;

        // Create array of forcefields
        // The forcefields in the are those of the fragments
        ForceField[] forcefield = new ForceField[nFragments];

        // Read from uniqueAtomNamesTextFiles? Or parent CIF/SDF
        // Currently read in from parent SDF
        String parentAtomNames[] = getParentAtomNames();

        // Loop over fragments.
        for (int i = 0; i < nFragments; i++) {

            // Read in force field patch
            // Need to read about FF constructor parameters
            // How to convert from SDF to FF?
            ForceField currentPatch = new ForceField(null);

            // Add new FF to array of fragment forcefields
            forcefield[i] = currentPatch;

            // Read in Rae's atom labels
            // <UniqueAtomName, Type>
            // TODO: make fragmentMap hashMap, read in atom type and name
            //HashMap<String, String> fragmentMap = new HashMap<>();
            HashMap<String, Integer> fragmentMap = new HashMap<>();

            // Match atom labels to AtomType instances
            // <type#, newType#>
            //HashMap<String, String> typeMap = createMap(parentAtomNames, fragmentMap);
            HashMap<Integer, Integer> typeMap = createMap(parentAtomNames, fragmentMap);

            // Overwrite AtomType name (i.e. PolType atom names are wrong)
            // Loop over all force field terms
            // If term includes only valid atoms (i.e. valid atom names) from the overall molecule, 
            // add it to parent force field
            // Angles
            AngleType fragAngleType = null;
            // Map
            int[] polTypeAngleClasses = fragAngleType.atomClasses;
            System.out.println("polTypeAngleClasses: ");
            for(int j = 0; j < polTypeAngleClasses.length; j++){
                System.out.println(polTypeAngleClasses[i]);
            }
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
                MultipoleType averageMultipoleType = MultipoleType.averageTypes(fragMultipoleType, parentMultipoleType, multipoleFrameTypes);
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
    /* Potential String Version
     private HashMap<String, String> createMap(String[] parentAtomNames, HashMap<String, String> fragmentMap) {
     HashMap<String, String> hashMap = new HashMap<>();
     int numParentAtoms = parentAtomNames.length;
     for (String key : fragmentMap.keySet()) {
     String type = fragmentMap.get(key);
     for (int i = 0; i < numParentAtoms; i++) {
     if (parentAtomNames[i].equalsIgnoreCase(key)) {
     String newType = type;
     hashMap.put(type, newType);
     }
     }
     }
     return hashMap;
     }*/

    private void updateAtomClasses(int currentTypes[], HashMap<Integer, Integer> map) {
        if (currentTypes == null) {
            return;
        }
        for (int i = 0; i < currentTypes.length; i++) {
            currentTypes[i] = map.get(currentTypes[i]);
        }
    }
    /* Potential String Version
     private void updateAtomClasses(String currentTypes[], HashMap<String, String> map) {
     if (currentTypes == null) {
     return;
     }
     for (int i = 0; i < currentTypes.length; i++) {
     currentTypes[i] = map.get(currentTypes[i]);
     }
     }*/

    private String[] getParentAtomNames() {

        // Read in parent file (sdfFiles.get(0))
        // Read atom names from it
        // May need to change to CIF read-in if unique atom names are necessary
        int atomCounter = 0;
        ArrayList<String> parentAtomNamesList = new ArrayList<>();

        try {
            BufferedReader read = new BufferedReader(new FileReader(sdfFiles.get(0)));
            String line;

            while ((line = read.readLine()) != null) {
                //test to see if the line read in contains unique atom name info.
                // Lines with atom names should end with 0  0  0  0  0
                if (line.contains("0  0  0  0  0")) {
                    atomCounter++;

                    String str32 = Character.toString(line.charAt(32));
                    String str33 = Character.toString(line.charAt(33));

                    String atomName = str32;

                    // If atom name is two letters (ex: Cl)
                    if (line.charAt(33) != ' ') {
                        atomName = str32.concat(str33);
                    }

                    // Add atom name to parentAtomNames array
                    parentAtomNamesList.add(atomName);

                }

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        // Convert parentAtomNamesList to parentAtomNames array to be passed back
        String[] parentAtomNames = new String[atomCounter];

        for (int i = 0; i < atomCounter; i++) {
            parentAtomNames[i] = parentAtomNamesList.get(i);
        }

        return parentAtomNames;

    }

}
