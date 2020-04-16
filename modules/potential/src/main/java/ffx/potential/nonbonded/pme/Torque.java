package ffx.potential.nonbonded.pme;

import static ffx.numerics.math.DoubleMath.X;
import static ffx.numerics.math.DoubleMath.add;
import static ffx.numerics.math.DoubleMath.dot;
import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.DoubleMath.normalize;
import static ffx.numerics.math.DoubleMath.scale;
import static ffx.numerics.math.DoubleMath.sub;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Torque {

  private static final Logger logger = Logger.getLogger(Torque.class.getName());

  private final double[] vecZ = new double[3];
  private final double[] vecX = new double[3];
  private final double[] vecY = new double[3];
  private final double[] sumXY = new double[3];
  private final double[] sumZXY = new double[3];
  private final double[] r = new double[3];
  private final double[] del = new double[3];
  private final double[] eps = new double[3];
  private final double[] xZXY = new double[3];
  private final double[] xXZ = new double[3];
  private final double[] xYZ = new double[3];
  private final double[] xYX = new double[3];
  private final double[] xXYZ = new double[3];
  private final double[] xxZXYZ = new double[3];
  private final double[] xxZXYX = new double[3];
  private final double[] xxZXYY = new double[3];
  private final double[] t1 = new double[3];
  private final double[] t2 = new double[3];
  private final double[] localOrigin = new double[3];

  private int[][] axisAtom;
  private MultipoleFrameDefinition[] frame;
  private double[][][] coordinates;

  public void init(int[][] axisAtom, MultipoleFrameDefinition[] frame, double[][][] coordinates) {
    this.axisAtom = axisAtom;
    this.frame = frame;
    this.coordinates = coordinates;
  }

  public void torque(int i, int iSymm, double[] trq, int[] frameIndex, double[][] g) {
    final int[] ax = axisAtom[i];
    // Ions, for example, have no torque.
    if (frame[i] == MultipoleType.MultipoleFrameDefinition.NONE) {
      return;
    }

    // Get the local frame type and the frame-defining atoms.
    int nAxisAtoms = ax.length;
    final int iz = ax[0];
    int ix = -1;
    int iy = -1;
    if (nAxisAtoms > 1) {
      ix = ax[1];
      if (nAxisAtoms > 2) {
        iy = ax[2];
      }
    }
    frameIndex[0] = iz;
    frameIndex[1] = ix;
    frameIndex[2] = iy;
    frameIndex[3] = i;
    double[] x = coordinates[iSymm][0];
    double[] y = coordinates[iSymm][1];
    double[] z = coordinates[iSymm][2];
    localOrigin[0] = x[i];
    localOrigin[1] = y[i];
    localOrigin[2] = z[i];
    vecZ[0] = x[iz];
    vecZ[1] = y[iz];
    vecZ[2] = z[iz];
    sub(vecZ, localOrigin, vecZ);
    double rZ = length(vecZ);
    if (frame[i] != MultipoleType.MultipoleFrameDefinition.ZONLY) {
      vecX[0] = x[ix];
      vecX[1] = y[ix];
      vecX[2] = z[ix];
      sub(vecX, localOrigin, vecX);
    } else {
      vecX[0] = 1.0;
      vecX[1] = 0.0;
      vecX[2] = 0.0;
      double dot = vecZ[0] / rZ;
      if (abs(dot) > 0.866) {
        vecX[0] = 0.0;
        vecX[1] = 1.0;
      }
    }
    double rX = length(vecX);
    if (frame[i] == MultipoleType.MultipoleFrameDefinition.ZTHENBISECTOR
        || frame[i] == MultipoleType.MultipoleFrameDefinition.THREEFOLD) {
      vecY[0] = x[iy];
      vecY[1] = y[iy];
      vecY[2] = z[iy];
      sub(vecY, localOrigin, vecY);
    } else {
      X(vecZ, vecX, vecY);
    }
    double rY = length(vecY);
    scale(vecZ, 1.0 / rZ, vecZ);
    scale(vecX, 1.0 / rX, vecX);
    scale(vecY, 1.0 / rY, vecY);
    // Find the perpendicular and angle for each pair of axes.
    X(vecX, vecZ, xXZ);
    X(vecY, vecZ, xYZ);
    X(vecY, vecX, xYX);
    normalize(xXZ, xXZ);
    normalize(xYZ, xYZ);
    normalize(xYX, xYX);
    // Compute the sine of the angle between the rotation axes.
    double cosZX = dot(vecZ, vecX);
    double sinZX = sqrt(1.0 - cosZX * cosZX);
    /*
     * Negative of dot product of torque with unit vectors gives
     * result of infinitesimal rotation along these vectors.
     */
    double dPhidZ = -(trq[0] * vecZ[0] + trq[1] * vecZ[1] + trq[2] * vecZ[2]);
    double dPhidX = -(trq[0] * vecX[0] + trq[1] * vecX[1] + trq[2] * vecX[2]);
    double dPhidY = -(trq[0] * vecY[0] + trq[1] * vecY[1] + trq[2] * vecY[2]);
    switch (frame[i]) {
      case ZONLY:
        for (int j = 0; j < 3; j++) {
          double dZ = xXZ[j] * dPhidX / (rZ * sinZX) + xYZ[j] * dPhidY / rZ;
          g[0][j] = dZ; // Atom Z
          g[3][j] = -dZ; // Atom I
        }
        break;
      case ZTHENX:
        for (int j = 0; j < 3; j++) {
          double dZ = xXZ[j] * dPhidX / (rZ * sinZX) + xYZ[j] * dPhidY / rZ;
          double dX = -xXZ[j] * dPhidZ / (rX * sinZX);
          g[0][j] = dZ; // Atom Z
          g[1][j] = dX; // Atom X
          g[3][j] = -(dZ + dX); // Atom I
        }
        break;
      case BISECTOR:
        for (int j = 0; j < 3; j++) {
          double dZ = xXZ[j] * dPhidX / (rZ * sinZX) + 0.5 * xYZ[j] * dPhidY / rZ;
          double dX = -xXZ[j] * dPhidZ / (rX * sinZX) + 0.5 * xYX[j] * dPhidY / rX;
          g[0][j] = dZ; // Atom Z
          g[1][j] = dX; // Atom X
          g[3][j] = -(dZ + dX); // Atom I
        }
        break;
      case ZTHENBISECTOR:
        // Build some additional axes needed for the Z-then-Bisector method
        add(vecX, vecY, sumXY);
        X(vecZ, sumXY, xZXY);
        normalize(sumXY, sumXY);
        normalize(xZXY, xZXY);
        // Find the perpendicular and angle for each pair of axes.
        X(sumXY, vecZ, xXYZ);
        X(xZXY, vecZ, xxZXYZ);
        X(xZXY, vecX, xxZXYX);
        X(xZXY, vecY, xxZXYY);
        normalize(xXYZ, xXYZ);
        normalize(xxZXYZ, xxZXYZ);
        normalize(xxZXYX, xxZXYX);
        normalize(xxZXYY, xxZXYY);
        // Compute the sine of the angle between the rotation axes
        double cosZXY = dot(vecZ, sumXY);
        double sinZXY = sqrt(1.0 - cosZXY * cosZXY);
        double cosXZXY = dot(vecX, xZXY);
        double sinXZXY = sqrt(1.0 - cosXZXY * cosXZXY);
        double cosYZXY = dot(vecY, xZXY);
        double sinYZXY = sqrt(1.0 - cosYZXY * cosYZXY);
        // Compute the projection of v and w onto the ru-plane
        scale(xZXY, -cosXZXY, t1);
        scale(xZXY, -cosYZXY, t2);
        add(vecX, t1, t1);
        add(vecY, t2, t2);
        normalize(t1, t1);
        normalize(t2, t2);
        double ut1cos = dot(vecZ, t1);
        double ut1sin = sqrt(1.0 - ut1cos * ut1cos);
        double ut2cos = dot(vecZ, t2);
        double ut2sin = sqrt(1.0 - ut2cos * ut2cos);
        double dPhidR = -(trq[0] * sumXY[0] + trq[1] * sumXY[1] + trq[2] * sumXY[2]);
        double dPhidZXY = -(trq[0] * xZXY[0] + trq[1] * xZXY[1] + trq[2] * xZXY[2]);
        for (int j = 0; j < 3; j++) {
          double dZ = xXYZ[j] * dPhidR / (rZ * sinZXY) + xxZXYZ[j] * dPhidZXY / rZ;
          double dX = (sinXZXY * xZXY[j] - cosXZXY * t1[j]) * dPhidZ / (rX * (ut1sin + ut2sin));
          double dY = (sinYZXY * xZXY[j] - cosYZXY * t2[j]) * dPhidZ / (rY * (ut1sin + ut2sin));
          g[0][j] = dZ; // Atom Z
          g[1][j] = dX; // Atom X
          g[2][j] = dY; // Atom Y
          g[3][j] = -(dZ + dX + dY); // Atom I
        }
        break;
      case THREEFOLD:
        add(vecZ, vecX, sumZXY);
        add(vecY, sumZXY, sumZXY);
        double rZXY = length(sumZXY);
        scale(sumZXY, 1.0 / rZXY, sumZXY);
        double cosYSum = dot(vecY, sumZXY);
        double cosZSum = dot(vecZ, sumZXY);
        double cosXSum = dot(vecX, sumZXY);
        add(vecX, vecY, r);
        normalize(r, r);
        double cosRZ = dot(r, vecZ);
        double sinRZ = sqrt(1.0 - cosRZ * cosRZ);
        dPhidR = -trq[0] * r[0] - trq[1] * r[1] - trq[2] * r[2];
        X(r, vecZ, del);
        normalize(del, del);
        double dPhidDel = -trq[0] * del[0] - trq[1] * del[1] - trq[2] * del[2];
        X(del, vecZ, eps);
        for (int j = 0; j < 3; j++) {
          double dZ = del[j] * dPhidR / (rZ * sinRZ) + eps[j] * dPhidDel * cosZSum / (rZ * rZXY);
          g[0][j] = dZ; // Atom Z
          g[3][j] -= dZ; // Atom I
        }
        add(vecZ, vecY, r);
        normalize(r, r);
        double cosRX = dot(r, vecX);
        double sinRX = sqrt(1.0 - cosRX * cosRX);
        dPhidR = -trq[0] * r[0] - trq[1] * r[1] - trq[2] * r[2];
        X(r, vecX, del);
        normalize(del, del);
        dPhidDel = -trq[0] * del[0] - trq[1] * del[1] - trq[2] * del[2];
        X(del, vecX, eps);
        for (int j = 0; j < 3; j++) {
          double dX = del[j] * dPhidR / (rX * sinRX) + eps[j] * dPhidDel * cosXSum / (rX * rZXY);
          g[1][j] = dX; // Atom X
          g[3][j] -= dX; // Atom I
        }
        add(vecZ, vecX, r);
        normalize(r, r);
        double cosRY = dot(r, vecY);
        double sinRY = sqrt(1.0 - cosRY * cosRY);
        dPhidR = -trq[0] * r[0] - trq[1] * r[1] - trq[2] * r[2];
        X(r, vecY, del);
        normalize(del, del);
        dPhidDel = -trq[0] * del[0] - trq[1] * del[1] - trq[2] * del[2];
        X(del, vecY, eps);
        for (int j = 0; j < 3; j++) {
          double dY = del[j] * dPhidR / (rY * sinRY) + eps[j] * dPhidDel * cosYSum / (rY * rZXY);
          g[2][j] = dY; // Atom Y
          g[3][j] = -dY; // Atom I
        }
        break;
      default:
        String message = "Fatal exception: Unsupported frame definition: " + frame[i] + "\n";
        logger.log(Level.SEVERE, message);
    }
  }
}
