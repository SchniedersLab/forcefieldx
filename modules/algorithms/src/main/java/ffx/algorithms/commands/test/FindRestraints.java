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
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsScript;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import groovy.lang.Binding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.String.format;

/**
 * The FindRestraints script identifies guest molecule atoms that should be restrained based on
 * their proximity to specific atoms in a host molecule.
 * <br>
 * Usage:
 * <br>
 * ffxc test.FindRestraints [options] &lt;filename&gt;
 */
@Command(description = " Find guest atoms to restrain near host molecule.", name = "test.FindRestraints")
public class FindRestraints extends AlgorithmsScript {

  /**
   * --hostName Molecule name of the host in the file.
   */
  @Option(names = {"--hostName"}, paramLabel = "BCD", defaultValue = "BCD",
      description = "Host molecule name in the file.")
  private String hostName;

  /**
   * --guestName Molecule name of the guest in the file.
   */
  @Option(names = {"--guestName"}, paramLabel = "LIG", defaultValue = "LIG",
      description = "Ligand molecule name in the file.")
  private String guestName;

  /**
   * --distanceCutoff Cutoff to use when selecting guest atoms near host COM.
   */
  @Option(names = {"--distanceCutoff"}, paramLabel = "5", defaultValue = "5",
      description = "Cutoff to use when selecting guest atoms near host COM")
  private double distanceCutoff;

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private List<String> filenames;

  /**
   * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
   * Originally it is declared in the run method
   */
  public Potential potential;

  /**
   * Get the potential object.
   * @return The potential object.
   */
  public Potential getPotentialObject() {
    return potential;
  }

  /**
   * FindRestraints Constructor.
   */
  public FindRestraints() {
    super();
  }

  /**
   * FindRestraints Constructor.
   * @param binding The Groovy Binding to use.
   */
  public FindRestraints(Binding binding) {
    super(binding);
  }

  /**
   * FindRestraints constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public FindRestraints(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public FindRestraints run() {
    if (!init()) {
      return this;
    }

    activeAssembly = getActiveAssembly(filenames.get(0));
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    Molecule[] molArr = activeAssembly.getMoleculeArray();

    List<Atom> restrainHostList = new ArrayList<>();
    List<Atom> restrainList = new ArrayList<>();
    double[] COM = new double[3];
    double[] subCOM = new double[3];
    int[] restrainHostIndices = new int[]{11, 16, 17, 20, 23, 26, 31, 32, 39, 40, 51, 63, 64, 70, 71, 82, 94, 95, 101, 102, 113, 125, 126, 132, 133, 144, 156, 157, 163, 164, 175, 187, 188, 191, 198};
    for (Molecule molecule : molArr) {
      logger.info(format(" Molecule name: " + molecule.getName()));
      if (molecule.getName().contains(hostName)) {
        Atom[] host_atoms = molecule.getAtomList().toArray(new Atom[0]);
        COM = getCOM(host_atoms);
        logger.info(format(" Center of mass of host molecule: " + Arrays.toString(COM)));
        for (Atom atom : host_atoms) {
          if (contains(restrainHostIndices, atom.getIndex())) {
            restrainHostList.add(atom);
            logger.info(format(" Atom: " + atom));
          }
        }
        Atom[] subAtoms = restrainHostList.toArray(new Atom[0]);
        subCOM = getCOM(subAtoms);
        logger.info(format(" Center of mass of subsection host atoms: " + Arrays.toString(subCOM)));
        double comdist = Math.sqrt(Math.pow(subCOM[0] - COM[0], 2) +
            Math.pow(subCOM[1] - COM[1], 2) +
            Math.pow(subCOM[2] - COM[2], 2));
        logger.info(format(" Distance between COMs: " + comdist));
      } else if (molecule.getName().contains(guestName)) {
        Atom[] guest_atoms = molecule.getAtomList().toArray(new Atom[0]);
        for (Atom atom : guest_atoms) {
          double dist = Math.sqrt(Math.pow(atom.getXYZ().get()[0] - subCOM[0], 2) +
              Math.pow(atom.getXYZ().get()[1] - subCOM[1], 2) +
              Math.pow(atom.getXYZ().get()[2] - subCOM[2], 2));
          //logger.info(format("Atom: " + atom));
          //logger.info(format("XYZ: " + atom.getXYZ().get()));
          //logger.info(format("Distance from host COM: " + dist));

          if (dist < distanceCutoff && atom.isHeavy()) {
            restrainList.add(atom);
          }
        }
      }
    }
    logger.info(format(" Number of atoms to restrain: " + restrainList.size()));
    int[] restrainIndices = restrainList.stream()
        .mapToInt(Atom::getIndex)
        .toArray();
    logger.info(format(" Restrain list indices: " + Arrays.toString(restrainIndices)));

    return this;
  }

  /**
   * Gets the center of mass of a set of atoms.
   *
   * @param atoms Array of atoms
   * @return x,y,z coordinates of center of mass
   */
  private static double[] getCOM(Atom[] atoms) {
    // Get center of mass of moleculeOneAtoms
    double[] COM = new double[3];
    double totalMass = 0.0;
    for (Atom s : atoms) {
      double[] pos = s.getXYZ().get();
      COM[0] += pos[0] * s.getMass();
      COM[1] += pos[1] * s.getMass();
      COM[2] += pos[2] * s.getMass();
      totalMass += s.getMass();
    }
    totalMass = 1 / totalMass;
    COM[0] *= totalMass;
    COM[1] *= totalMass;
    COM[2] *= totalMass;

    return COM;
  }

  /**
   * Checks if an int array contains a specific value.
   *
   * @param array Array to search
   * @param value Value to find
   * @return true if array contains value
   */
  private static boolean contains(int[] array, int value) {
    return Arrays.stream(array).anyMatch(element -> element == value);
  }
}
