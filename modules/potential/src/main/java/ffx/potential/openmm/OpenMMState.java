package ffx.potential.openmm;

import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_get;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

/**
 * Retrieve state information from an OpenMM Simulation.
 */
public class OpenMMState {
  /**
   * Logger.
   */
  private static final Logger logger = Logger.getLogger(OpenMMState.class.getName());

  /**
   * Potential energy (kcal/mol).
   */
  public final double potentialEnergy;
  /**
   * Kinetic energy (kcal/mol).
   */
  public final double kineticEnergy;
  /**
   * Total energy (kcal/mol).
   */
  public final double totalEnergy;
  /**
   * Pointer to an OpenMM state.
   */
  private final PointerByReference state;
  /**
   * Mask of information to retrieve.
   */
  private final int mask;
  /**
   * Array of atoms.
   */
  private final Atom[] atoms;
  /**
   * Degrees of freedom.
   */
  private final int n;
  /**
   * Number of atoms.
   */
  private final int nAtoms;

  /**
   * Construct an OpenMM State with the requested information.
   *
   * @param openMMContext OpenMM Context.
   * @param mask          Mask of information to retrieve.
   * @param atoms         Array of atoms.
   */
  public OpenMMState(OpenMMContext openMMContext, int mask, Atom[] atoms, int dof) {
    this.mask = mask;
    this.atoms = atoms;
    this.n = dof;
    nAtoms = atoms.length;

    // Retrieve the state.
    state = openMMContext.getState(mask);

    if (stateContains(OpenMM_State_Energy)) {
      potentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
      kineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
      totalEnergy = potentialEnergy + kineticEnergy;
    } else {
      potentialEnergy = 0.0;
      kineticEnergy = 0.0;
      totalEnergy = 0.0;
    }
  }

  /**
   * Check to see if the state contains the requested information.
   *
   * @param flag Information to check for.
   * @return boolean indicating whether the state contains the requested information.
   */
  private boolean stateContains(int flag) {
    return (mask & flag) == flag;
  }

  public void free() {
    if (state != null) {
      logger.fine(" Free OpenMM State.");
      OpenMM_State_destroy(state);
      logger.fine(" Free OpenMM State completed.");
    }
  }

  /**
   * The force array contains the OpenMM force information for all atoms. The returned array a
   * contains accelerations for active atoms only.
   *
   * @param a Acceleration components for only active atomic coordinates.
   * @return The acceleration for each active atomic coordinate.
   */
  public double[] getAccelerations(@Nullable double[] a) {
    if (!stateContains(OpenMM_State_Forces)) {
      return a;
    }

    if (a == null || a.length != n) {
      a = new double[n];
    }

    PointerByReference forcePointer = OpenMM_State_getForces(state);

    for (int i = 0, index = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        double mass = atom.getMass();
        OpenMM_Vec3 acc = OpenMM_Vec3Array_get(forcePointer, i);
        double xx = acc.x * OpenMM_AngstromsPerNm / mass;
        double yy = acc.y * OpenMM_AngstromsPerNm / mass;
        double zz = acc.z * OpenMM_AngstromsPerNm / mass;
        a[index++] = xx;
        a[index++] = yy;
        a[index++] = zz;
        atom.setAcceleration(xx, yy, zz);
      }
    }
    return a;
  }

  /**
   * The force array contains the OpenMM force information for all atoms. The returned array g
   * contains components for active atoms only.
   *
   * @param g Gradient components for only active atomic coordinates.
   * @return g The gradient includes only active atoms
   */
  public double[] getGradient(@Nullable double[] g) {
    if (!stateContains(OpenMM_State_Forces)) {
      return g;
    }

    if (g == null || g.length != n) {
      g = new double[n];
    }

    PointerByReference forcePointer = OpenMM_State_getForces(state);
    for (int i = 0, index = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        OpenMM_Vec3 force = OpenMM_Vec3Array_get(forcePointer, i);
        double xx = -force.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        double yy = -force.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        double zz = -force.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
        if (isNaN(xx) || isInfinite(xx) || isNaN(yy) || isInfinite(yy) || isNaN(zz) || isInfinite(
            zz)) {
          StringBuilder sb = new StringBuilder(
              format(" The gradient of atom %s is (%8.3f,%8.3f,%8.3f).", atom, xx, yy, zz));
          double[] vals = new double[3];
          atom.getVelocity(vals);
          sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          atom.getAcceleration(vals);
          sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          atom.getPreviousAcceleration(vals);
          sb.append(
              format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          throw new EnergyException(sb.toString());
        }
        g[index++] = xx;
        g[index++] = yy;
        g[index++] = zz;
        atom.setXYZGradient(xx, yy, zz);
      }
    }
    return g;
  }

  /**
   * Read the periodic lattice vectors from a state.
   *
   * <p>The crystal instance will be updated, and passed to the ForceFieldEnergy instance.
   */
  public double[][] getPeriodicBoxVectors() {
    if (!stateContains(OpenMM_State_Positions)) {
      return null;
    }

    OpenMM_Vec3 a = new OpenMM_Vec3();
    OpenMM_Vec3 b = new OpenMM_Vec3();
    OpenMM_Vec3 c = new OpenMM_Vec3();
    OpenMM_State_getPeriodicBoxVectors(state, a, b, c);
    double[][] latticeVectors = new double[3][3];
    latticeVectors[0][0] = a.x * OpenMM_AngstromsPerNm;
    latticeVectors[0][1] = a.y * OpenMM_AngstromsPerNm;
    latticeVectors[0][2] = a.z * OpenMM_AngstromsPerNm;
    latticeVectors[1][0] = b.x * OpenMM_AngstromsPerNm;
    latticeVectors[1][1] = b.y * OpenMM_AngstromsPerNm;
    latticeVectors[1][2] = b.z * OpenMM_AngstromsPerNm;
    latticeVectors[2][0] = c.x * OpenMM_AngstromsPerNm;
    latticeVectors[2][1] = c.y * OpenMM_AngstromsPerNm;
    latticeVectors[2][2] = c.z * OpenMM_AngstromsPerNm;
    return latticeVectors;
  }

  /**
   * The positions array contains the OpenMM atomic position information for all atoms. The
   * returned array x contains coordinates only for active atoms.
   *
   * @param x Atomic coordinates only for active atoms.
   * @return x The atomic coordinates for only active atoms.
   */
  public double[] getPositions(@Nullable double[] x) {
    if (!stateContains(OpenMM_State_Positions)) {
      return x;
    }

    if (x == null || x.length != n) {
      x = new double[n];
    }

    PointerByReference positionsPointer = OpenMM_State_getPositions(state);
    for (int i = 0, index = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        OpenMM_Vec3 pos = OpenMM_Vec3Array_get(positionsPointer, i);
        double xx = pos.x * OpenMM_AngstromsPerNm;
        double yy = pos.y * OpenMM_AngstromsPerNm;
        double zz = pos.z * OpenMM_AngstromsPerNm;
        x[index++] = xx;
        x[index++] = yy;
        x[index++] = zz;
        atom.moveTo(xx, yy, zz);
      }
    }
    return x;
  }

  /**
   * The positions array contains the OpenMM atomic position information for all atoms. The
   * returned array x contains coordinates for active atoms only.
   *
   * @param v Velocity only for active atomic coordinates.
   * @return v The velocity for each active atomic coordinate.
   */
  public double[] getVelocities(@Nullable double[] v) {
    if (!stateContains(OpenMM_State_Velocities)) {
      return v;
    }

    if (v == null || v.length != n) {
      v = new double[n];
    }

    PointerByReference velocitiesPointer = OpenMM_State_getVelocities(state);
    for (int i = 0, index = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      if (atom.isActive()) {
        OpenMM_Vec3 vel = OpenMM_Vec3Array_get(velocitiesPointer, i);
        double xx = vel.x * OpenMM_AngstromsPerNm;
        double yy = vel.y * OpenMM_AngstromsPerNm;
        double zz = vel.z * OpenMM_AngstromsPerNm;
        v[index++] = xx;
        v[index++] = yy;
        v[index++] = zz;
        atom.setVelocity(xx, yy, zz);
      }
    }
    return v;
  }
}
