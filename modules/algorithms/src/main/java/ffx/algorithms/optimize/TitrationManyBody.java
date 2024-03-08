// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
// ******************************************************************************
package ffx.algorithms.optimize;

import static ffx.potential.bonded.RotamerLibrary.applyRotamer;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.openmm.OpenMMEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.openmm.OpenMMEnergy;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TitrationUtils;
import ffx.potential.parsers.PDBFilter;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import java.lang.String;

public class TitrationManyBody {

  private static final Logger logger = Logger.getLogger(TitrationManyBody.class.getName());

  private final ForceField forceField;
  private final List<Integer> resNumberList;
  private final double pH;
  private final String filename;
  private PDBFilter protonFilter;
  private ForceFieldEnergy potentialEnergy;

  public TitrationManyBody(String filename, ForceField forceField, List<Integer> resNumberList,
      double pH) {
    this.filename = filename;
    this.forceField = forceField;
    this.resNumberList = resNumberList;
    this.pH = pH;
  }

  public MolecularAssembly getProtonatedAssembly() {
    MolecularAssembly protonatedAssembly = new MolecularAssembly(filename);
    protonatedAssembly.setForceField(forceField);
    File structureFile = new File(filename);
    protonFilter = new PDBFilter(structureFile, protonatedAssembly, forceField,
        forceField.getProperties(), resNumberList);
    protonFilter.setRotamerTitration(true);
    protonFilter.readFile();
    protonFilter.applyAtomProperties();
    protonatedAssembly.finalize(true, forceField);
    potentialEnergy = ForceFieldEnergy.energyFactory(protonatedAssembly);
    protonatedAssembly.setFile(structureFile);

    TitrationUtils titrationUtils;
    titrationUtils = new TitrationUtils(protonatedAssembly.getForceField());
    titrationUtils.setRotamerPhBias(298.15, pH);
    for (Residue residue : protonatedAssembly.getResidueList()) {
      String resName = residue.getName();
      if (resNumberList.contains(residue.getResidueNumber())) {
        if (resName.equalsIgnoreCase("ASH") || resName.equalsIgnoreCase("GLH")
            || resName.equalsIgnoreCase("LYS") || resName.equalsIgnoreCase("HIS")
            || resName.equalsIgnoreCase("CYS")) {
          residue.setTitrationUtils(titrationUtils);
        }
      }
    }
    if (potentialEnergy instanceof OpenMMEnergy openMMEnergy) {
      boolean updateBondedTerms = forceField.getBoolean("TITRATION_UPDATE_BONDED_TERMS", true);
      openMMEnergy.getSystem().setUpdateBondedTerms(updateBondedTerms);
    }
    potentialEnergy.energy();
    return protonatedAssembly;
  }

  public MolecularAssembly[] getProtonatedAssemblies() {
    logger.info("Getting protonated assemblies");
    MolecularAssembly molecularAssembly = getProtonatedAssembly();
    List<Character> altLocs = protonFilter.getAltLocs();
    int locs = 1;
    if(altLocs!=null){
      locs = altLocs.size();
      for (int i = 0; i < locs; i++) {
        if (altLocs.get(i) == null) {
          altLocs.remove(altLocs.get(i));
        }
      }
    }
    MolecularAssembly[] molecularAssemblies = new MolecularAssembly[locs];
    molecularAssemblies[0] = molecularAssembly;
    for (int i = 0; i < altLocs.size(); i++) {
      if (i != 0) {
        logger.info(filename);
        MolecularAssembly newAssembly = new MolecularAssembly(filename);
        newAssembly.setForceField(forceField);
        File structureFile = new File(filename);
        protonFilter = new PDBFilter(structureFile, newAssembly, forceField,
            forceField.getProperties(), resNumberList);
        logger.info(newAssembly.getResidueList().toString());
        protonFilter.setRotamerTitration(true);
        protonFilter.setAltID(newAssembly, altLocs.get(i));
        protonFilter.readFile();
        logger.info(newAssembly.getResidueList().get(0).getAtomList().get(0).getAltLoc().toString());
        protonFilter.applyAtomProperties();
        newAssembly.finalize(true, forceField);
        potentialEnergy = ForceFieldEnergy.energyFactory(newAssembly);

        TitrationUtils titrationUtils;
        titrationUtils = new TitrationUtils(molecularAssembly.getForceField());
        titrationUtils.setRotamerPhBias(298.15, pH);
        for (Residue residue : molecularAssembly.getResidueList()) {
          String resName = residue.getName();
          if (resNumberList.contains(residue.getResidueNumber())) {
            if (resName.equalsIgnoreCase("ASH") || resName.equalsIgnoreCase("GLH")
                || resName.equalsIgnoreCase("LYS") || resName.equalsIgnoreCase("HIS")
                || resName.equalsIgnoreCase("CYS")) {
              residue.setTitrationUtils(titrationUtils);
            }
          }
        }
        potentialEnergy.energy();
        molecularAssemblies[i] = newAssembly;
      }
    }
    return molecularAssemblies;
  }


  public boolean excludeExcessAtoms(Set<Atom> excludeAtoms, int[] optimalRotamers,
      List<Residue> residueList) {
    boolean isTitrating = false;
    int i = 0;
    for (Residue residue : residueList) {
      Rotamer rotamer = residue.getRotamers()[optimalRotamers[i++]];
      applyRotamer(residue, rotamer);
      if (rotamer.isTitrating) {
        isTitrating = true;
        AminoAcidUtils.AminoAcid3 aa3 = rotamer.aminoAcid3;
        residue.setName(aa3.name());
        switch (aa3) {
          case HID, GLU -> {
            // No HE2
            Atom HE2 = residue.getAtomByName("HE2", true);
            excludeAtoms.add(HE2);
          }
          case HIE -> {
            // No HD1
            Atom HD1 = residue.getAtomByName("HD1", true);
            excludeAtoms.add(HD1);
          }
          case ASP -> {
            // No HD2
            Atom HD2 = residue.getAtomByName("HD2", true);
            excludeAtoms.add(HD2);
          }
          case LYD -> {
            // No HZ3
            Atom HZ3 = residue.getAtomByName("HZ3", true);
            excludeAtoms.add(HZ3);
          }
          case CYD -> {
            // No HG
            Atom HG = residue.getAtomByName("HG", true);
            excludeAtoms.add(HG);
          }
          default -> {
          }
          // Do nothing.
        }
      }
    }

    return isTitrating;
  }


}
