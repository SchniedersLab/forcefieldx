// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.cli;

import static ffx.utilities.StringUtils.parseAtomRanges;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.logging.Logger;
import picocli.CommandLine.Option;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;

/**
 * Represents command line options for scripts that support atom selections.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AtomSelectionOptions {

  private static final Logger logger = Logger.getLogger(AtomSelectionOptions.class.getName());

  /** --aa or --activeAtoms Ranges of active atoms [NONE, ALL, Range(s): 1-3,6-N]. */
  @Option(
      names = {"--aa", "--active"},
      paramLabel = "<selection>",
      defaultValue = "",
      description = "Ranges of active atoms [NONE, ALL, Range(s): 1-3,6-N].")
  public String activeAtoms;

  /** --ia or --inactiveAtoms Ranges of inactive atoms [NONE, ALL, Range(s): 1-3,6-N]. */
  @Option(
      names = {"--ia", "--inactive"},
      paramLabel = "<selection>",
      defaultValue = "",
      description = "Ranges of inactive atoms [NONE, ALL, Range(s): 1-3,6-N].")
  public String inactiveAtoms;

  /**
   * Set active atoms for a MolecularAssembly.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  public void setActiveAtoms(MolecularAssembly molecularAssembly) {
    // First, evaluate the inactive atom selection.
    setInactive(molecularAssembly);

    // Second, evaluate the active atom selection, which takes precedence over the inactive flag.
    setActive(molecularAssembly);
  }

  private void setInactive(MolecularAssembly assembly) {
    actOnAtoms(assembly, inactiveAtoms, (Atom a, Boolean b) -> a.setActive(!b),
            "inactive", " Inactive atoms");
  }

  private void setActive(MolecularAssembly assembly) {
    actOnAtoms(assembly, activeAtoms, Atom::setActive, "active", " Active atoms");
  }

  public static void actOnAtoms(@Nonnull MolecularAssembly assembly, @Nullable String selection,
                                @Nonnull BiConsumer<Atom, Boolean> action, @Nonnull String description,
                                @Nullable String keyType) {
    if (selection == null || selection.equalsIgnoreCase("")) {
      // Empty or null string -- no changes.
      return;
    }

    Atom[] atoms = assembly.getAtomArray();

    // No atoms selected.
    if (selection.equalsIgnoreCase("NONE")) {
      for (Atom atom : atoms) {
        action.accept(atom, false);
      }
      logger.info(" No atoms are " + description + ".\n");
      return;
    }

    // All atoms selected
    if (selection.equalsIgnoreCase("ALL")) {
      for (Atom atom : atoms) {
        action.accept(atom, true);
      }
      logger.info(" All atoms are " + description + ".\n");
      return;
    }

    // A range(s) of atoms are active.
    int nAtoms = atoms.length;
    for (Atom atom : atoms) {
      action.accept(atom, false);
    }
    if (keyType != null) {
      List<Integer> atomRanges = parseAtomRanges(keyType, selection, nAtoms);
      for (int i : atomRanges) {
        action.accept(atoms[i], true);
      }
    }
    logger.info(" " + description + " atoms set to: " + selection + "\n");
  }
}
