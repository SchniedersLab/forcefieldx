// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.bonded;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.DoubleMath.length;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * Restrain the position of atoms to their initial coordinates.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RestrainPosition extends BondedTerm implements LambdaInterface {

  private static final Logger logger = Logger.getLogger(RestrainPosition.class.getName());

  private final double[][] equilibriumCoordinates;

  /**
   * Force constant variable stores K/2 in Kcal/mol/A. E = K/2 * dx^2.
   */
  private final double forceConstant;

  /**
   * Flat bottom radius in Angstroms.
   */
  private final double flatBottom;

  private final double[] a1 = new double[3];
  private final double[] dx = new double[3];
  private final double lambdaExp = 1.0;
  private final double[] lambdaGradient;
  private final int nAtoms;
  private final boolean lambdaTerm;
  private double lambda = 1.0;
  private double lambdaPow = pow(lambda, lambdaExp);
  private double dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
  private double d2LambdaPow = 0;
  private double dEdL = 0.0;
  private double d2EdL2 = 0.0;

  /**
   * Restrain atoms to a position in the global coordinate frame.
   *
   * @param atoms                  The atoms to restrain.
   * @param equilibriumCoordinates The equilibrium coordinates.
   * @param forceConst             The force constant in kcal/mol/A^2.
   * @param flatBottom             The flat bottom radius in Angstroms.
   * @param lambdaTerm             If true, apply lambda to this restraint.
   */
  public RestrainPosition(Atom[] atoms, double[][] equilibriumCoordinates,
                          double forceConst, double flatBottom, boolean lambdaTerm) {
    this.atoms = atoms;
    this.equilibriumCoordinates = equilibriumCoordinates;
    this.forceConstant = forceConst;
    this.flatBottom = flatBottom;
    nAtoms = atoms.length;
    this.lambdaTerm = lambdaTerm;
    if (lambdaTerm) {
      lambdaGradient = new double[nAtoms * 3];
    } else {
      lambdaGradient = null;
      lambda = 1.0;
      lambdaPow = 1.0;
      dLambdaPow = 0.0;
      d2LambdaPow = 0.0;
    }

    if (logger.isLoggable(Level.FINE)) {
      logger.info(" RestrainPosition: " + lambdaTerm);
      for (Atom atom : atoms) {
        logger.fine(atom.toString());
      }
    }
  }

  /**
   * Returns a copy of the atoms array.
   *
   * @return Copy of the atom array.
   */
  public Atom[] getAtoms() {
    Atom[] retArray = new Atom[nAtoms];
    arraycopy(atoms, 0, retArray, 0, nAtoms);
    return retArray;
  }

  /**
   * Returns the force constant in kcal/mol/Angstrom^2.
   *
   * @return a double.
   */
  public double getForceConstant() {
    return forceConstant;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getLambda() {
    return lambda;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setLambda(double lambda) {
    if (lambdaTerm) {
      this.lambda = lambda;
      double lambdaWindow = 1.0;
      if (this.lambda <= lambdaWindow) {
        double dldgl = 1.0 / lambdaWindow;
        double l = dldgl * this.lambda;
        double l2 = l * l;
        double l3 = l2 * l;
        double l4 = l2 * l2;
        double l5 = l4 * l;
        double c3 = 10.0;
        double c4 = -15.0;
        double c5 = 6.0;
        double threeC3 = 3.0 * c3;
        double sixC3 = 6.0 * c3;
        double fourC4 = 4.0 * c4;
        double twelveC4 = 12.0 * c4;
        double fiveC5 = 5.0 * c5;
        double twentyC5 = 20.0 * c5;
        lambdaPow = c3 * l3 + c4 * l4 + c5 * l5;
        dLambdaPow = (threeC3 * l2 + fourC4 * l3 + fiveC5 * l4) * dldgl;
        d2LambdaPow = (sixC3 * l + twelveC4 * l2 + twentyC5 * l3) * dldgl * dldgl;
      } else {
        lambdaPow = 1.0;
        dLambdaPow = 0.0;
        d2LambdaPow = 0.0;
      }
    } else {
      this.lambda = 1.0;
      lambdaPow = 1.0;
      dLambdaPow = 0.0;
      d2LambdaPow = 0.0;
    }
  }

  /**
   * getNumAtoms.
   *
   * @return a int.
   */
  public int getNumAtoms() {
    return nAtoms;
  }

  /**
   * Returns the original coordinates of this restraint, indexed by atoms then x,y,z. This is the
   * opposite order of the internal storage.
   *
   * @return Original coordinates [atoms][xyz]
   */
  public double[][] getEquilibriumCoordinates() {
    double[][] equilibriumCoords = new double[nAtoms][3];
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        // Mild subtlety here: it is stored internally as [xyz][atoms], but returned as [atoms][xyz]
        equilibriumCoords[i][j] = equilibriumCoordinates[j][i];
      }
    }
    return equilibriumCoords;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    if (lambdaTerm) {
      return d2EdL2;
    } else {
      return 0.0;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    if (lambdaTerm) {
      return dEdL;
    } else {
      return 0.0;
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] gradient) {
    if (lambdaTerm) {
      int n3 = nAtoms * 3;
      for (int i = 0; i < n3; i++) {
        gradient[i] += lambdaGradient[i];
      }
    }
  }


  /**
   * Calculates energy and gradients for this coordinate restraint.
   *
   * @param gradient Calculate gradients
   * @return Energy in the coordinate restraint
   */
  @Override
  public double energy(boolean gradient, int threadID, AtomicDoubleArray3D grad, AtomicDoubleArray3D lambdaGrad) {
    double e = 0.0;
    dEdL = 0.0;
    d2EdL2 = 0.0;
    double fx2 = forceConstant * 2.0;
    for (int i = 0; i < nAtoms; i++) {
      // Current atomic coordinates.
      Atom atom = atoms[i];
      // Compute the separation distance.
      atom.getXYZ(a1);
      dx[0] = a1[0] - equilibriumCoordinates[0][i];
      dx[1] = a1[1] - equilibriumCoordinates[1][i];
      dx[2] = a1[2] - equilibriumCoordinates[2][i];
      double r = length(dx);
      double dr = r;
      if (flatBottom > 0.0) {
        dr = Math.max(0.0, r - flatBottom);
      }
      value = dr;
      double dr2 = dr * dr;
      e += dr2;
      if (gradient || lambdaTerm) {
        double scale = fx2 * dr;
        if (r > 0.0) {
          scale /= r;
        }
        final double dedx = dx[0] * scale;
        final double dedy = dx[1] * scale;
        final double dedz = dx[2] * scale;
        final int index = atom.getIndex() - 1;
        if (gradient) {
          grad.add(threadID, index, lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
        }
        if (lambdaTerm) {
          lambdaGrad.add(threadID, index, dLambdaPow * dedx, dLambdaPow * dedy, dLambdaPow * dedz);
        }
      }
    }

    if (lambdaTerm) {
      dEdL = dLambdaPow * forceConstant * e;
      d2EdL2 = d2LambdaPow * forceConstant * e;
    }

    energy = forceConstant * e * lambdaPow;

    // logger.info(format(" RestrainPosition: %12.6f %12.6f %12.6f", energy, dEdL, d2EdL2));

    return energy;
  }

  public static RestrainPosition[] parseRestrainPositions(MolecularAssembly molecularAssembly) {
    List<RestrainPosition> restrainPositionList = new ArrayList<>();
    ForceField forceField = molecularAssembly.getForceField();
    CompositeConfiguration properties = forceField.getProperties();
    if (properties.containsKey("restrain-position")) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      String[] lines = properties.getStringArray("restrain-position");
      for (String line : lines) {
        RestrainPosition restrain = parseRestrainPosition(line, atoms, false);
        if (restrain != null) {
          restrainPositionList.add(restrain);
        }
      }
    }
    if (properties.containsKey("restrain-position-lambda")) {
      Atom[] atoms = molecularAssembly.getAtomArray();
      String[] lines = properties.getStringArray("restrain-position-lambda");
      for (String line : lines) {
        RestrainPosition restrain = parseRestrainPosition(line, atoms, true);
        if (restrain != null) {
          restrainPositionList.add(restrain);
        }
      }
    }

    if (restrainPositionList.isEmpty()) {
      return null;
    }

    return restrainPositionList.toArray(new RestrainPosition[0]);
  }

  /**
   * Parse a Restrain-Position line and return a RestrainPosition instance.
   *
   * @param line      The restraint line.
   * @param atoms     The atoms in the system.
   * @param useLambda If true, apply lambda to this restraint.
   * @return A RestrainPosition object.
   */
  public static RestrainPosition parseRestrainPosition(String line, Atom[] atoms, boolean useLambda) {
    String[] tokens = line.split("\\s+");
    if (tokens.length < 1) {
      throw new IllegalArgumentException("Invalid restraint line: " + line);
    }

    int nAtoms = atoms.length;
    Atom[] atomArray;
    double[][] coordinates;
    double forceConstant = 50.0;
    double flatBottom = 0.0;
    // First token is the atom index.
    int start = parseInt(tokens[0]);
    if (start < 0) {
      // Range of atoms.
      start = -start - 1;
      int end = parseInt(tokens[1]) - 1;
      if (start >= nAtoms || end < start || end >= nAtoms) {
        logger.severe(format(" Property restrain-position could not be parsed:\n %s.", line));
        return null;
      }
      int n = end - start + 1;
      atomArray = new Atom[n];
      coordinates = new double[3][n];
      for (int j = start, index = 0; j <= end; j++, index++) {
        atomArray[index] = atoms[j];
        coordinates[0][index] = atoms[j].getX();
        coordinates[1][index] = atoms[j].getY();
        coordinates[2][index] = atoms[j].getZ();
      }
      if (tokens.length > 2) {
        forceConstant = parseDouble(tokens[2]);
      }
      if (tokens.length > 3) {
        flatBottom = parseDouble(tokens[3]);
      }
      logger.fine(format(" Restrain-Position of Atoms %d to %d (K=%8.4f, D=%8.4f)",
          start + 1, end + 1, forceConstant, flatBottom));
    } else {
      // Single atom.
      Atom atom = atoms[start - 1];
      atomArray = new Atom[1];
      atomArray[0] = atom;
      double x = atom.getX();
      double y = atom.getY();
      double z = atom.getZ();
      if (tokens.length > 1) {
        x = parseDouble(tokens[1]);
      }
      if (tokens.length > 2) {
        y = parseDouble(tokens[2]);
      }
      if (tokens.length > 3) {
        z = parseDouble(tokens[3]);
      }
      if (tokens.length > 4) {
        forceConstant = parseDouble(tokens[4]);
      }
      if (tokens.length > 5) {
        flatBottom = parseDouble(tokens[5]);
      }
      logger.fine(format(
          " Restrain-Position of %s to (%12.6f, %12.6f, %12.6f) with K=%8.4f and D=%8.4f",
          atom, x, y, z, forceConstant, flatBottom));
      coordinates = new double[3][1];
      coordinates[0][0] = x;
      coordinates[1][0] = y;
      coordinates[2][0] = z;
    }
    return new RestrainPosition(atomArray, coordinates, forceConstant, flatBottom, useLambda);
  }
}
