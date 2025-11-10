//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.algorithms.commands;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.algorithms.optimize.TitrationManyBody;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.RotamerLibrary.NucleicSugarPucker;
import ffx.potential.parameters.TitrationUtils;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXBinding;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.getExtension;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * The SaveRotamers script saves out rotamers.
 * <br>
 * Usage:
 * <br>
 * ffxc SaveRotamers [options] &lt;filename&gt;
 */
@Command(description = " Save out rotamers.", name = "SaveRotamers")
public class SaveRotamers extends AlgorithmsCommand {

  @Option(names = {"--chain", "-c"}, paramLabel = " ",
      description = "Single character chain name.")
  private char c = ' ';

  @Option(names = {"--library", "-l"}, paramLabel = "1", defaultValue = "1",
      description = "Available rotamer libraries are (1) Ponder and Richards or (2) Richardson.")
  private int library;

  @Option(names = {"--resid", "-r"}, paramLabel = "1", defaultValue = "1",
      description = "Residue number.")
  private int resID;

  @Option(names = {"--independent", "-i"}, paramLabel = "false", defaultValue = "false",
      description = "Independent draws nucleic acid rotamers independently of chain context.")
  private boolean independent;

  @Option(names = {"--start", "-s"}, paramLabel = "0", defaultValue = "0",
      description = "First rotamer to draw (indexed from rotamer 0).")
  private int start;

  @Option(names = {"--finish", "-f"}, paramLabel = "-1", defaultValue = "-1",
      description = "Last rotamer to draw (indexed from rotamer 0).")
  private int finish;

  @Option(names = {"--all", "-x"}, paramLabel = "-1", defaultValue = "-1",
      description = "Draw all rotamers beginning from the passed index (overrides other options).")
  private int all;

  @Option(names = {"--upstreamPucker", "-u"}, paramLabel = "false", defaultValue = "false",
      description = "Adjusts the pucker of the 5' residue to match the rotamer.")
  private boolean upstreamPucker;

  @Option(names = {"--tR", "--titrateResidue"}, paramLabel = "false", defaultValue = "false",
      description = "Titrate residues.")
  private boolean titrateResidue;

  @Parameters(arity = "1", paramLabel = "file",
      description = "The atomic coordinate file in XYZ or PDB format.")
  private String filename = null;

  private TitrationManyBody titrationManyBody;
  private TitrationUtils titrationUtils;

  public SaveRotamers() {
    super();
  }

  public SaveRotamers(FFXBinding binding) {
    super(binding);
  }

  public SaveRotamers(String[] args) {
    super(args);
  }

  @Override
  public SaveRotamers run() {
    if (!init()) {
      return this;
    }

    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    filename = activeAssembly.getFile().getAbsolutePath();
    CompositeConfiguration properties = activeAssembly.getProperties();

    if (titrateResidue) {
      List<Residue> residues = activeAssembly.getResidueList();
      List<Integer> resNumberList = new ArrayList<>();
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber());
      }
      titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
          resNumberList, 7.0);
      activeAssembly = titrationManyBody.getProtonatedAssembly();
    }
    RotamerLibrary rLib = new RotamerLibrary(
        RotamerLibrary.ProteinLibrary.intToProteinLibrary(library), true);
    String chain = String.valueOf(c);

    boolean saveAllRotamers = false;
    int allStart = 0;
    if (all > -1) {
      allStart = all;
      saveAllRotamers = true;
    }

    logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".");

    RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains());
    Polymer polymer = activeAssembly.getChain(chain);
    if (polymer == null) {
      logger.info(" Polymer + " + chain + " does not exist.");
      return this;
    }
    Residue residue = polymer.getResidue(resID);
    if (residue == null) {
      logger.info(" Residue + " + resID + " does not exist.");
      return this;
    }

    residue.setRotamers(rLib);
    Rotamer[] rotamers = residue.getRotamers();
    if (rotamers == null) {
      logger.severe(" There are no rotamers for residue + " + residue.toString());
    }
    int nrotamers = rotamers.length;

    boolean isDeoxy = false;
    Residue prevResidue = null;
    if (upstreamPucker) {
      try {
        if (residue.getResidueType() == Residue.ResidueType.NA) {
          prevResidue = residue.getPreviousResidue();
          if (prevResidue == null) {
            upstreamPucker = false;
          } else {
            Atom HOs = (Atom) prevResidue.getAtomNode("HO'");
            if (HOs == null) {
              isDeoxy = true;
            }
          }
        } else {
          upstreamPucker = false;
        }
      } catch (Exception e) {
        upstreamPucker = false;
      }
    }

    String ext = getExtension(filename);
    filename = removeExtension(filename);
    Set<Atom> excludeAtoms = new HashSet<>();
    List<Residue> removeAtomList = new ArrayList<>();
    if (saveAllRotamers) {
      if (allStart >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.");
      } else {
        for (int i = allStart; i < nrotamers; i++) {
          RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
          logger.info(rotamers[i].getName());
          if (upstreamPucker) {
            double prevDelta = rotamers[i].chi1;
            NucleicSugarPucker sugarPucker = NucleicSugarPucker.checkPucker(prevDelta, isDeoxy);
            RotamerLibrary.applySugarPucker(prevResidue, sugarPucker, isDeoxy, true);
          }
          if (ext.toUpperCase().contains("XYZ")) {
            logger.info(format(" Saving rotamer %d", i));
            algorithmFunctions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"));
          } else {
            if (titrateResidue) {
              excludeAtoms.clear();
              removeAtomList.clear();
              removeAtomList.add(residue);

              int[] currentRot = new int[]{i};
              titrationManyBody.excludeExcessAtoms(excludeAtoms, currentRot, removeAtomList);
              File modelFile = new File(filename + ".pdb");
              PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
                  properties);
              if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
                logger.info(format(" Save failed for %s", activeAssembly));
              }
            } else {
              logger.info(format(" Saving rotamer %d", i));
              algorithmFunctions.saveAsPDB(activeAssembly, new File(filename + ".pdb"));
            }
          }
        }
      }
    } else {
      if (start >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.");
      } else {
        if (finish >= nrotamers) {
          finish = nrotamers - 1;
        } else if (finish < start) {
          logger.info(
              " Specified finish point is before the start point drawing only rotamer " + start);
          finish = start;
        }
        for (int i = start; i <= finish; i++) {
          RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
          if (upstreamPucker) {
            double prevDelta = rotamers[i].chi1;
            NucleicSugarPucker sugarPucker = NucleicSugarPucker.checkPucker(prevDelta, isDeoxy);
            RotamerLibrary.applySugarPucker(prevResidue, sugarPucker, isDeoxy, true);
          }
          if (ext.toUpperCase().contains("XYZ")) {
            logger.info(format(" Saving rotamer %d", i));
            algorithmFunctions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"));
          } else {
            logger.info(format(" Saving rotamer %d", i));
            algorithmFunctions.saveAsPDB(activeAssembly, new File(filename + ".pdb"));
          }
        }
      }
    }

    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    return Collections.emptyList();
  }
}
