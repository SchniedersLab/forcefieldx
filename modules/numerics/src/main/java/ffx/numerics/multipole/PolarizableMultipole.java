package ffx.numerics.multipole;

public class PolarizableMultipole {

  private static final double oneThird = 1.0 / 3.0;
  private static final double twoThirds = 2.0 / 3.0;

  /**
   * Partial charge.
   */
  protected double q;
  /**
   * Dipole x-component.
   */
  protected double dx;
  /**
   * Dipole y-component.
   */
  protected double dy;
  /**
   * Dipole z-component.
   */
  protected double dz;
  /**
   * Quadrupole xx-component multiplied by 1/3.
   */
  protected double qxx;
  /**
   * Quadrupole yy-component multiplied by 1/3.
   */
  protected double qyy;
  /**
   * Quadrupole zz-component multiplied by 1/3.
   */
  protected double qzz;
  /**
   * Quadrupole xy-component multiplied by 2/3.
   */
  protected double qxy;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected double qxz;
  /**
   * Quadrupole xz-component multiplied by 2/3.
   */
  protected double qyz;
  /**
   * Induced dipole x-component.
   */
  protected double ux;
  /**
   * Induced dipole y-component.
   */
  protected double uy;
  /**
   * Induced dipole z-component.
   */
  protected double uz;
  /**
   * Induced dipole chain rule x-component.
   */
  protected double px;
  /**
   * Induced dipole chain rule y-component.
   */
  protected double py;
  /**
   * Induced dipole chain rule z-component.
   */
  protected double pz;
  /**
   * Averaged induced dipole + induced dipole chain-rule x-component: sx = 0.5 * (ux + px).
   */
  protected double sx;
  /**
   * Averaged induced dipole + induced dipole chain-rule y-component: sy = 0.5 * (uy + py).
   */
  protected double sy;
  /**
   * Averaged induced dipole + induced dipole chain-rule z-component: sz = 0.5 * (uz + pz).
   */
  protected double sz;

  /**
   * PolarizableMultipole constructor.
   *
   * @param Q Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public PolarizableMultipole(double[] Q, double[] u, double[] uCR) {
    setMultipole(Q);
    setInducedDipole(u, uCR);
  }

  /**
   * Set the permanent multipole.
   *
   * @param Q Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   * @param u Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public void set(double[] Q, double[] u, double[] uCR) {
    setMultipole(Q);
    setInducedDipole(u, uCR);
  }

  /**
   * Set the permanent multipole.
   *
   * @param Q Multipole Q[q, dx, dy, dz, qxx, qyy, qzz, qxy, qxz, qyz]
   */
  public void setMultipole(double[] Q) {
    q = Q[0];
    dx = Q[1];
    dy = Q[2];
    dz = Q[3];
    qxx = Q[4] * oneThird;
    qyy = Q[5] * oneThird;
    qzz = Q[6] * oneThird;
    qxy = Q[7] * twoThirds;
    qxz = Q[8] * twoThirds;
    qyz = Q[9] * twoThirds;
  }

  /**
   * Set the induced dipole.
   *
   * @param u Induced dipole u[ux, uy, uz]
   * @param uCR Induced dipole chain-rule uCR[ux, uy, uz]
   */
  public void setInducedDipole(double[] u, double[] uCR) {
    ux = u[0];
    uy = u[1];
    uz = u[2];
    px = uCR[0];
    py = uCR[1];
    pz = uCR[2];
    sx = 0.5 * (ux + px);
    sy = 0.5 * (uy + py);
    sz = 0.5 * (uz + pz);
  }

  /**
   * Clear the induced dipoles.
   */
  public void clearInducedDipoles() {
    ux = 0.0;
    uy = 0.0;
    uz = 0.0;
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    sx = 0.0;
    sy = 0.0;
    sz = 0.0;
  }

  /**
   * Scale the averaged induced dipole.
   *
   * @param scaleInduction a double.
   * @param scaleEnergy a double.
   */
  public final void applyMasks(double scaleInduction, double scaleEnergy) {
    // [Ux, Uy, Uz] resulted from induction masking rules, and we now apply the energy mask.
    // [Px, Py, Pz] resulted from energy masking rules, and we now apply the induction mask.
    sx = 0.5 * (ux * scaleEnergy + px * scaleInduction);
    sy = 0.5 * (uy * scaleEnergy + py * scaleInduction);
    sz = 0.5 * (uz * scaleEnergy + pz * scaleInduction);
  }

}
