//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.numerics.Potential;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import java.util.*;
import static java.lang.Math.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import static java.lang.String.format;

/**
 * The FindRestraints script identifies guest molecule atoms that should be restrained based on
 * their proximity to specific atoms in a host molecule.
 * <br>
 * Usage:
 * <br>
 * ffxc test.FindRestraints [options] &lt;filename&gt;
 */
@Command(description = "Find guest atoms for COM or Boresch restraints.",
        name = "test.FindRestraints")
public class FindRestraints extends AlgorithmsCommand {

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

  // ------------------------
  // Boresch Options
  // ------------------------

  @Option(names = {"--boresch"}, defaultValue = "false",
      description = "Use Boresch anchor selection instead of COM-based selection.")
  private boolean boresch;

  // From GHOAT paper, H1, H2, H3 are explicitly defined by the user
  @Option(names = {"--H1"}, paramLabel = "int",
      description = "Host anchor H1 atom index.")
  private Integer H1Index;

  @Option(names = {"--H2"}, paramLabel = "int",
      description = "Host anchor H2 atom index.")
  private Integer H2Index;

  @Option(names = {"--H3"}, paramLabel = "int",
      description = "Host anchor H3 atom index.")
  private Integer H3Index;


  @Option(names = {"--l1Range"}, defaultValue = "3.5",
      description = "Cylinder diameter/height for G1 search.")
  private double l1Range;

  @Option(names = {"--minAdis"}, defaultValue = "1.5",
      description = "Minimum distance between guest anchors.")
  private double minAdis;

  @Option(names = {"--maxAdis"}, defaultValue = "5.0",
      description = "Maximum distance between guest anchors.")
  private double maxAdis;

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
   * @param binding The Binding to use.
   */
  public FindRestraints(FFXBinding binding) {
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
    Molecule host = null;
    Molecule guest = null;

    for (Molecule mol : activeAssembly.getMoleculeArray()) {
      if (mol.getName().contains(hostName)) host = mol;
      if (mol.getName().contains(guestName)) guest = mol;
    }

    if (host == null || guest == null) {
      logger.severe("Host or Guest molecule not found.");
      return this;
    }

    if (boresch) {
      runBoreschMode(host, guest);
    } else {
      runCOMMode();
    }

    return this;
  }

  // ================================================================
  // ==================== COM MODE =================================
  // ================================================================
  private void runCOMMode() {   
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

  }

  // ================================================================
  // ==================== BORESCH MODE ==============================
  // ================================================================

    /**
     * Calculate rotation matrix to align host to Z-axis.
     *
     * The returned matrix has rows = {x̂, ŷ, ẑ} so that r' = R * r.
     */
    private double[][] getAlignmentRotationMatrix(Atom[] hostAtoms, double[] com) {
      // Compute inertia tensor
      double[][] I = new double[3][3];

      for (Atom atom : hostAtoms) {
        double m = atom.getMass();
        double[] r = atom.getXYZ().get();
        // Use coordinates relative to COM
        double rx = r[0] - com[0];
        double ry = r[1] - com[1];
        double rz = r[2] - com[2];

        I[0][0] += m * (ry * ry + rz * rz);
        I[1][1] += m * (rx * rx + rz * rz);
        I[2][2] += m * (rx * rx + ry * ry);

        I[0][1] -= m * rx * ry;
        I[0][2] -= m * rx * rz;
        I[1][2] -= m * ry * rz;
      }

      I[1][0] = I[0][1];
      I[2][0] = I[0][2];
      I[2][1] = I[1][2];

      // Diagonalize inertia tensor
      EigenDecomposition ed = new EigenDecomposition(new Array2DRowRealMatrix(I));
      double[] evals = ed.getRealEigenvalues();
      double[][] evecs = ed.getV().getData(); // columns are eigenvectors

      int zIndex = 0;
      int xIndex = 0;
      for (int i = 1; i < 3; i++) {
        if (evals[i] > evals[zIndex]) {
          zIndex = i;
        }
        if (evals[i] < evals[xIndex]) {
          xIndex = i;
        }
      }
      int yIndex = 3 - zIndex - xIndex;

      double[] zHat = normalize(new double[]{evecs[0][zIndex], evecs[1][zIndex], evecs[2][zIndex]});
      double[] xHat = normalize(new double[]{evecs[0][xIndex], evecs[1][xIndex], evecs[2][xIndex]});
      double[] yHat = normalize(cross(zHat, xHat));
      xHat = normalize(cross(yHat, zHat));

      return new double[][]{xHat, yHat, zHat};
    }

  /**
   * Apply translation and rotation transformation to atoms
   */
    private void applyTransformation(Atom[] atoms, double[] com, double[][] rotationMatrix) {
      for (Atom atom : atoms) {
        double[] r = atom.getXYZ().get();

        // Translate to origin
        double tx = r[0] - com[0];
        double ty = r[1] - com[1];
        double tz = r[2] - com[2];

        // Rotate: r' = R * r
        double[] rotated = new double[3];
        for (int i = 0; i < 3; i++) {
          rotated[i] = rotationMatrix[i][0] * tx
              + rotationMatrix[i][1] * ty
              + rotationMatrix[i][2] * tz;
        }

        atom.setXYZ(rotated);
      }
    }

    private static double[][] captureCoordinates(Atom[] atoms) {
      double[][] coords = new double[atoms.length][3];
      for (int i = 0; i < atoms.length; i++) {
        double[] r = atoms[i].getXYZ().get();
        coords[i][0] = r[0];
        coords[i][1] = r[1];
        coords[i][2] = r[2];
      }
      return coords;
    }

    private static void restoreCoordinates(Atom[] atoms, double[][] coords) {
      for (int i = 0; i < atoms.length; i++) {
        atoms[i].setXYZ(coords[i]);
      }
    }

  private void runBoreschMode(Molecule host, Molecule guest) {

    if (H1Index == null || H2Index == null || H3Index == null) {
      logger.severe("Boresch mode requires --H1 --H2 --H3 indices.");
      return;
    }

    Atom H1 = findAtomByIndex(host, H1Index);
    Atom H2 = findAtomByIndex(host, H2Index);
    Atom H3 = findAtomByIndex(host, H3Index);

    if (H1 == null || H2 == null || H3 == null) {
      logger.severe("Could not locate host anchor atoms.");
      return;
    }

    // Calculate principal axis from HOST non-hydrogen atoms ONLY
    Atom[] hostAtoms = host.getAtomList().toArray(new Atom[0]);
    Atom[] guestAtoms = guest.getAtomList().toArray(new Atom[0]);
    
    // Use COM of non-hydrogen host atoms as the origin
    double[] com = getCOMNonHydrogen(hostAtoms);
    logger.info(format("Host COM (non-H): (%.3f, %.3f, %.3f)", com[0], com[1], com[2]));
    
    double[][] rotationMatrix = getAlignmentRotationMatrix(hostAtoms, com);

    // Temporarily align to select anchors in a consistent frame
    double[][] hostCoords = captureCoordinates(hostAtoms);
    double[][] guestCoords = captureCoordinates(guestAtoms);
    applyTransformation(hostAtoms, com, rotationMatrix);
    applyTransformation(guestAtoms, com, rotationMatrix);

    Atom G1 = selectG1(guest);
    if (G1 == null) return;
    
    Atom G2 = selectG2(guest, G1);
    if (G2 == null) return;
    
    Atom G3 = selectG3(guest, G1, G2);
    if (G3 == null) return;

    // Restore original coordinates for geometry evaluation
    restoreCoordinates(hostAtoms, hostCoords);
    restoreCoordinates(guestAtoms, guestCoords);

    // Log selected atoms
    logger.info(format("\n=== Selected Boresch Anchors ==="));
    logger.info(format("Host anchors:"));
    logger.info(format("  H1: %s (index %d, residue %d)", H1.getName(), H1.getIndex(), H1.getResidueNumber()));
    logger.info(format("  H2: %s (index %d, residue %d)", H2.getName(), H2.getIndex(), H2.getResidueNumber()));
    logger.info(format("  H3: %s (index %d, residue %d)", H3.getName(), H3.getIndex(), H3.getResidueNumber()));
    logger.info(format("Guest anchors:"));
    logger.info(format("  G1: %s (index %d, residue %d)", G1.getName(), G1.getIndex(), G1.getResidueNumber()));
    logger.info(format("  G2: %s (index %d, residue %d)", G2.getName(), G2.getIndex(), G2.getResidueNumber()));
    logger.info(format("  G3: %s (index %d, residue %d)", G3.getName(), G3.getIndex(), G3.getResidueNumber()));

    // Note: G1, G2, G3 are geometric anchor atoms and do not require bonding constraints

    // Calculate all geometric parameters
    double r = distance(H1.getXYZ().get(), G1.getXYZ().get());
    double thetaA = angle(H2, H1, G1);
    double thetaB = angle(H1, G1, G2);
    double phiA = dihedral(H3, H2, H1, G1);
    double phiB = dihedral(H2, H1, G1, G2);
    double phiC = dihedral(H1, G1, G2, G3);

    logger.info(format("\n=== Boresch Restraint Parameters ==="));
    logger.info(format("Distance H1-G1 (r): %.3f Å", r));
    logger.info(format("Angle H2-H1-G1 (θA): %.2f°", thetaA));
    logger.info(format("Angle H1-G1-G2 (θB): %.2f°", thetaB));
    logger.info(format("Dihedral H3-H2-H1-G1 (φA): %.2f°", phiA));
    logger.info(format("Dihedral H2-H1-G1-G2 (φB): %.2f°", phiB));
    logger.info(format("Dihedral H1-G1-G2-G3 (φC): %.2f°", phiC));
    
    // Validate geometry (warn about near-degenerate angles)
    if (thetaA < 20.0 || thetaA > 160.0) {
        logger.warning(format("Angle θA = %.2f° is near-degenerate.", thetaA));
    }
    if (thetaB < 20.0 || thetaB > 160.0) {
        logger.warning(format("Angle θB = %.2f° is near-degenerate.", thetaB));
    }
  }



  /**
   * Calculate dihedral angle between four atoms
   */
  private static double dihedral(Atom a, Atom b, Atom c, Atom d) {
      double[] v1 = subtract(a.getXYZ().get(), b.getXYZ().get());
      double[] v2 = subtract(c.getXYZ().get(), b.getXYZ().get());
      double[] v3 = subtract(d.getXYZ().get(), c.getXYZ().get());
      
      double[] n1 = cross(v1, v2);
      double[] n2 = cross(v2, v3);
      
      double[] m1 = cross(n1, v2);
      
      double x = dot(n1, n2);
      double y = dot(m1, n2);
      
      return toDegrees(atan2(y, x));
  }

  private static double[] cross(double[] a, double[] b) {
      return new double[]{
          a[1]*b[2] - a[2]*b[1],
          a[2]*b[0] - a[0]*b[2],
          a[0]*b[1] - a[1]*b[0]
      };
  }

  private static double dot(double[] a, double[] b) {
      return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }


  // ================================================================
  // ==================== GUEST SELECTION ===========================
  // ================================================================

  private Atom selectG1(Molecule guest) {
      double radius = l1Range / 2.0;
      double radius2 = radius * radius;
      double halfHeight = l1Range / 2.0;

      // MISSING: You need to declare these variables
      Atom G1 = null;
      double bestR2 = Double.MAX_VALUE;

      // MISSING: You need to iterate through guest atoms
      for (Atom atom : guest.getAtomList()) {
          if (!atom.isHeavy()) continue;

          double[] xyz = atom.getXYZ().get();
          double r2 = xyz[0]*xyz[0] + xyz[1]*xyz[1];

          if (Math.abs(xyz[2]) <= halfHeight && r2 <= radius2) {
              if (r2 < bestR2) {
                  bestR2 = r2;
                  G1 = atom;
              }
          }
      }
      
      if (G1 == null) {
          logger.severe("anch1 error: No valid G1 anchor found.");
          return null;  
      }

      return G1;  
  }

  private Atom selectG2(Molecule guest, Atom G1) {

    Atom best = null;
    double bestProj = Double.MAX_VALUE;

    double[] g1 = G1.getXYZ().get();

    for (Atom atom : guest.getAtomList()) {
      if (atom == G1 || !atom.isHeavy()) continue;

      double d = distance(g1, atom.getXYZ().get());
      if (d < minAdis || d > maxAdis) continue;

      double[] diff = subtract(atom.getXYZ().get(), g1);
      double proj2 = diff[0] * diff[0] + diff[1] * diff[1];

      if (proj2 < bestProj) {
        best = atom;
        bestProj = proj2;
      }
    }
    if (best == null) {
        logger.severe("anch2 error: Could not identify valid G2 anchor.");
    }

    return best;
  }

  private Atom selectG3(Molecule guest, Atom G1, Atom G2) {
      Atom best = null;
      double bestAngleDiff = Double.MAX_VALUE;

      List<Atom> candidates = new ArrayList<>();

      for (Atom atom : guest.getAtomList()) {
          if (atom == G1 || atom == G2 || !atom.isHeavy()) continue;

          // Check distance from G2
          double d = distance(G2.getXYZ().get(), atom.getXYZ().get());
          if (d < minAdis || d > maxAdis) continue;

          // Check distance from G1
          double d_g1 = distance(G1.getXYZ().get(), atom.getXYZ().get());
          if (d_g1 < minAdis) continue;

          candidates.add(atom);
      }
      
      if (candidates.isEmpty()) {
          logger.warning("No valid G3 candidates found");
          return null;
      }
      
      // Find atom with angle closest to 90°
      for (Atom atom : candidates) {
          double ang = angle(G1, G2, atom);
          
          // Reject degenerate angles
          if (ang < 20.0 || ang > 160.0) continue;
          
          double diff = Math.abs(ang - 90.0);
          
          if (diff < bestAngleDiff) {
              bestAngleDiff = diff;
              best = atom;
          }
      }
      
      if (best == null) {
          logger.warning("No G3 candidate with acceptable angle found");
      } else {
          logger.info(String.format("Selected G3: %s (index %d)", 
              best.getName(), best.getIndex()));
      }
      
      return best;
  }

  // ================================================================
  // ==================== GEOMETRY HELPERS ==========================
  // ================================================================

  private static double distance(double[] a, double[] b) {
    return sqrt(pow(a[0] - b[0], 2)
            + pow(a[1] - b[1], 2)
            + pow(a[2] - b[2], 2));
  }

  private static double[] subtract(double[] a, double[] b) {
    return new double[]{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  private static double[] normalize(double[] v) {
    double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm == 0.0) {
      return new double[]{0.0, 0.0, 0.0};
    }
    return new double[]{v[0] / norm, v[1] / norm, v[2] / norm};
  }

  private static double angle(Atom a, Atom b, Atom c) {
    double[] v1 = subtract(a.getXYZ().get(), b.getXYZ().get());
    double[] v2 = subtract(c.getXYZ().get(), b.getXYZ().get());

    double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double mag1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double mag2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

    return toDegrees(acos(dot / (mag1 * mag2)));
  }

  private static Atom findAtomByIndex(Molecule mol, int index) {
    for (Atom atom : mol.getAtomList()) {
      if (atom.getIndex() == index) return atom;
    }
    return null;
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
   * Gets the center of mass of non-hydrogen atoms in a set.
   *
   * @param atoms Array of atoms
   * @return x,y,z coordinates of center of mass (non-hydrogen only)
   */
  private static double[] getCOMNonHydrogen(Atom[] atoms) {
    double[] COM = new double[3];
    double totalMass = 0.0;
    for (Atom s : atoms) {
      if (!s.isHeavy()) continue;  // Skip hydrogen atoms
      double[] pos = s.getXYZ().get();
      COM[0] += pos[0] * s.getMass();
      COM[1] += pos[1] * s.getMass();
      COM[2] += pos[2] * s.getMass();
      totalMass += s.getMass();
    }
    if (totalMass == 0.0) {
      return new double[]{0.0, 0.0, 0.0};
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
