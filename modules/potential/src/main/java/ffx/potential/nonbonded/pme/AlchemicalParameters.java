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
package ffx.potential.nonbonded.pme;

import edu.rit.pj.IntegerSchedule;
import edu.rit.util.Range;
import ffx.crystal.Crystal;
import ffx.potential.parameters.ForceField;

import java.util.logging.Logger;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;

public class AlchemicalParameters {

  private static final Logger logger = Logger.getLogger(AlchemicalParameters.class.getName());

  /**
   * The polarization model to use.
   */
  private final Polarization polarization;
  /**
   * The alchemical mode to use.
   */
  public final AlchemicalMode mode;
  /**
   * Constant α in: r' = sqrt(r^2 + α*(1 - L)^2)
   */
  public double permLambdaAlpha = 1.0;
  /**
   * Power on L in front of the pairwise multipole potential.
   */
  public double permLambdaExponent = 3.0;
  /**
   * Begin turning on permanent multipoles at Lambda = 0.4;
   */
  public double permLambdaStart = 0.4;
  /**
   * Finish turning on permanent multipoles at Lambda = 1.0;
   */
  public double permLambdaEnd = 1.0;
  /**
   * Start turning on polarization later in the Lambda path to prevent SCF convergence problems when
   * atoms nearly overlap.
   */
  public double polLambdaStart = 0.75;
  public double polLambdaEnd = 1.0;
  /**
   * Power on L in front of the polarization energy.
   */
  public double polLambdaExponent = 3.0;
  /**
   * Intramolecular electrostatics for the ligand in vapor is included by default.
   */
  public boolean doLigandVaporElec = true;
  /**
   * Intramolecular electrostatics for the ligand in done in GK implicit solvent.
   */
  public boolean doLigandGKElec = false;
  /**
   * Condensed phase SCF without the ligand present is included by default. For DualTopologyEnergy
   * calculations it can be turned off.
   */
  public boolean doNoLigandCondensedSCF = true;
  /**
   * lAlpha = α*(1 - L)^2
   */
  public double lAlpha = 0.0;
  public double dlAlpha = 0.0;
  public double d2lAlpha = 0.0;
  public double dEdLSign = 1.0;
  /**
   * lPowPerm = L^permanentLambdaExponent
   */
  public double lPowPerm = 1.0;
  public double dlPowPerm = 0.0;
  public double d2lPowPerm = 0.0;
  public boolean doPermanentRealSpace = true;
  public double permanentScale = 1.0;
  /**
   * lPowPol = L^polarizationLambdaExponent
   */
  public double lPowPol = 1.0;
  public double dlPowPol = 0.0;
  public double d2lPowPol = 0.0;
  public boolean doPolarization = true;
  /**
   * When computing the polarization energy at L there are 3 pieces.
   *
   * <p>1.) Upol(1) = The polarization energy computed normally (ie. system with ligand).
   *
   * <p>2.) Uenv = The polarization energy of the system without the ligand.
   *
   * <p>3.) Uligand = The polarization energy of the ligand by itself.
   *
   * <p>Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
   *
   * <p>Set polarizationScale to L for part 1. Set polarizationScale to (1-L) for parts 2 and 3.
   */
  public double polarizationScale = 1.0;
  /**
   * The polarization Lambda value goes from 0.0 .. 1.0 as the global lambda value varies between
   * polLambdaStart .. polLambadEnd.
   */
  public double polLambda = 1.0;
  /**
   * The permanent Lambda value goes from 0.0 .. 1.0 as the global lambda value varies between
   * permLambdaStart .. permLambdaEnd.
   */
  public double permLambda = 1.0;
  /**
   * Boundary conditions for the vapor end of the alchemical path.
   */
  public Crystal vaporCrystal = null;
  public int[][][] vaporLists = null;
  public Range[] vacuumRanges = null;
  public IntegerSchedule vaporPermanentSchedule = null;
  public IntegerSchedule vaporEwaldSchedule = null;

  public AlchemicalParameters(ForceField forceField, boolean lambdaTerm,
                              boolean nnTerm, Polarization polarization) {
    this.polarization = polarization;

    AlchemicalMode tempMode;
    try {
      tempMode = AlchemicalMode.valueOf(forceField.getString("ALCHEMICAL_MODE", "OST").toUpperCase());
    } catch (IllegalArgumentException e) {
      logger.info(" Invalid value for alchemical-mode; reverting to OST");
      tempMode = AlchemicalMode.OST;
    }
    mode = tempMode;

    if (lambdaTerm) {
      // Values of PERMANENT_LAMBDA_ALPHA below 2 can lead to unstable  trajectories.
      permLambdaAlpha = forceField.getDouble("PERMANENT_LAMBDA_ALPHA", 2.0);
      if (permLambdaAlpha < 0.0 || permLambdaAlpha > 3.0) {
        logger.warning("Invalid value for permanent-lambda-alpha (<0.0 || >3.0); reverting to 2.0");
        permLambdaAlpha = 2.0;
      }

      /*
       A PERMANENT_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at the
       beginning of the permanent schedule. Choosing a power of 3 or
       greater ensures a smooth dU/dL and d2U/dL2 over the schedule.

       A value of 0.0 is also admissible for when ExtendedSystem is
       scaling multipoles rather than softcoring them.
      */
      permLambdaExponent = forceField.getDouble("PERMANENT_LAMBDA_EXPONENT", 3.0);
      if (permLambdaExponent < 0.0) {
        logger.warning("Invalid value for permanent-lambda-exponent (<0.0); reverting to 3.0");
        permLambdaExponent = 3.0;
      }

      /*
       A POLARIZATION_LAMBDA_EXPONENT of 2 gives a non-zero d2U/dL2 at
       the beginning of the polarization schedule. Choosing a power of 3
       or greater ensures a smooth dU/dL and d2U/dL2 over the schedule.

       A value of 0.0 is also admissible: when polarization is not being
       softcored but instead scaled, as by ExtendedSystem.
      */
      polLambdaExponent = forceField.getDouble("POLARIZATION_LAMBDA_EXPONENT", 3.0);
      if (polLambdaExponent < 0.0) {
        logger.warning("Invalid value for polarization-lambda-exponent (<0.0); reverting to 3.0");
        polLambdaExponent = 3.0;
      }

      // Values of PERMANENT_LAMBDA_START below 0.5 can lead to unstable trajectories.
      permLambdaStart = forceField.getDouble("PERMANENT_LAMBDA_START", 0.4);
      if (permLambdaStart < 0.0 || permLambdaStart > 1.0) {
        logger.warning("Invalid value for perm-lambda-start (<0.0 || >1.0); reverting to 0.4");
        permLambdaStart = 0.4;
      }

      // Values of PERMANENT_LAMBDA_END must be greater than permLambdaStart and <= 1.0.
      permLambdaEnd = forceField.getDouble("PERMANENT_LAMBDA_END", 1.0);
      if (permLambdaEnd < permLambdaStart || permLambdaEnd > 1.0) {
        logger.warning("Invalid value for perm-lambda-end (<start || >1.0); reverting to 1.0");
        permLambdaEnd = 1.0;
      }

      /*
         The POLARIZATION_LAMBDA_START defines the point in the lambda
         schedule when the condensed phase polarization of the ligand
         begins to be turned on. If the condensed phase polarization
         is considered near lambda=0, then SCF convergence is slow,
         even with Thole damping. In addition, 2 (instead of 1)
         condensed phase SCF calculations are necessary from the
         beginning of the window to lambda=1.
        */
      polLambdaStart = forceField.getDouble("POLARIZATION_LAMBDA_START", 0.75);
      if (polLambdaStart < 0.0 || polLambdaStart > 1.0) {
        logger.warning("Invalid value for polarization-lambda-start; reverting to 0.75");
        polLambdaStart = 0.75;
      }

      /*
         The POLARIZATION_LAMBDA_END defines the point in the lambda
         schedule when the condensed phase polarization of ligand has
         been completely turned on. Values other than 1.0 have not been tested.
        */
      polLambdaEnd = forceField.getDouble("POLARIZATION_LAMBDA_END", 1.0);
      if (polLambdaEnd < polLambdaStart || polLambdaEnd > 1.0) {
        logger.warning(
            "Invalid value for polarization-lambda-end (<start || >1.0); reverting to 1.0");
        polLambdaEnd = 1.0;
      }

      // The LAMBDA_VAPOR_ELEC defines if intra-molecular electrostatics of the ligand in vapor
      // will be considered.
      if (mode == AlchemicalMode.OST) {
        doLigandVaporElec = forceField.getBoolean("LIGAND_VAPOR_ELEC", true);
        doLigandGKElec = forceField.getBoolean("LIGAND_GK_ELEC", false);
        doNoLigandCondensedSCF = forceField.getBoolean("NO_LIGAND_CONDENSED_SCF", true);
      } else {
        doLigandVaporElec = false;
        doLigandGKElec = false;
        doNoLigandCondensedSCF = false;
      }
    } else if (nnTerm) {
      permLambdaAlpha = 0.0;
      permLambdaExponent = 1.0;
      polLambdaExponent = 1.0;
      permLambdaStart = 0.0;
      permLambdaEnd = 1.0;
      polLambdaStart = 0.0;
      polLambdaEnd = 1.0;
      // The LAMBDA_VAPOR_ELEC defines if intramolecular electrostatics of the neural network
      // atoms in vapor will be removed.
      if (mode == AlchemicalMode.OST) {
        doLigandVaporElec = forceField.getBoolean("LIGAND_VAPOR_ELEC", true);
      } else {
        doLigandVaporElec = false;
      }
      doLigandGKElec = false;
      doNoLigandCondensedSCF = false;
    }
  }

  public String toString() {
    StringBuilder sb = new StringBuilder("   Alchemical Parameters\n");
    sb.append(
        format(
            "    Permanent Multipole Range:      %5.3f-%5.3f\n", permLambdaStart, permLambdaEnd));
    sb.append(format("    Permanent Multipole Softcore Alpha:   %5.3f\n", permLambdaAlpha));
    sb.append(format("    Permanent Multipole Lambda Exponent:  %5.3f\n", permLambdaExponent));
    if (polarization != Polarization.NONE) {
      sb.append(format("    Polarization Lambda Exponent:         %5.3f\n", polLambdaExponent));
      sb.append(
          format(
              "    Polarization Range:             %5.3f-%5.3f\n", polLambdaStart, polLambdaEnd));
      sb.append(format("    Condensed SCF Without Ligand:         %B\n", doNoLigandCondensedSCF));
    }
    if (!doLigandGKElec) {
      sb.append(format("    Vapor Electrostatics:                 %B\n", doLigandVaporElec));
    } else {
      sb.append(format("    GK Electrostatics at L=0:             %B\n", doLigandGKElec));
    }
    return sb.toString();
  }

  /*
   * f = sqrt(r^2 + lAlpha)
   *
   * df/dL = -alpha * (1.0 - lambda) / f
   *
   * g = 1 / sqrt(r^2 + lAlpha)
   *
   * dg/dL = alpha * (1.0 - lambda) / (r^2 + lAlpha)^(3/2)
   *
   * define dlAlpha = alpha * 1.0 - lambda)
   *
   * then df/dL = -dlAlpha / f and dg/dL = dlAlpha * g^3
   *
   * Multipoles are turned on from permLambdaStart .. permLambdaEnd.
   *
   * @param lambda
   */
  public void update(double lambda) {
    lPowPerm = 1.0;
    permLambda = 1.0;
    dlPowPerm = 0.0;
    d2lPowPerm = 0.0;
    lAlpha = 0.0;
    dlAlpha = 0.0;
    d2lAlpha = 0.0;
    if (lambda < permLambdaStart) {
      lPowPerm = 0.0;
      permLambda = 0.0;
    } else if (lambda <= permLambdaEnd) {
      double permWindow = permLambdaEnd - permLambdaStart;
      double permLambdaScale = 1.0 / permWindow;
      permLambda = permLambdaScale * (lambda - permLambdaStart);
      if (mode == AlchemicalMode.OST) {
        lAlpha = permLambdaAlpha * (1.0 - permLambda) * (1.0 - permLambda);
        dlAlpha = permLambdaAlpha * (1.0 - permLambda);
        d2lAlpha = -permLambdaAlpha;

        lPowPerm = pow(permLambda, permLambdaExponent);
        dlPowPerm = permLambdaExponent * pow(permLambda, permLambdaExponent - 1.0);
        d2lPowPerm = 0.0;
        if (permLambdaExponent >= 2.0) {
          d2lPowPerm =
              permLambdaExponent
                  * (permLambdaExponent - 1.0)
                  * pow(permLambda, permLambdaExponent - 2.0);
        }

        dlAlpha *= permLambdaScale;
        d2lAlpha *= (permLambdaScale * permLambdaScale);
        dlPowPerm *= permLambdaScale;
        d2lPowPerm *= (permLambdaScale * permLambdaScale);
      }
    }

    // Polarization is turned on from polarizationLambdaStart .. polarizationLambdaEnd.
    lPowPol = 1.0;
    polLambda = 1.0;
    dlPowPol = 0.0;
    d2lPowPol = 0.0;
    if (lambda < polLambdaStart) {
      lPowPol = 0.0;
      polLambda = 0.0;
    } else if (lambda <= polLambdaEnd) {
      double polWindow = polLambdaEnd - polLambdaStart;
      double polLambdaScale = 1.0 / polWindow;
      polLambda = polLambdaScale * (lambda - polLambdaStart);
      if (mode == AlchemicalMode.OST) {
        if (polLambdaExponent > 0.0) {
          lPowPol = pow(polLambda, polLambdaExponent);
          if (polLambdaExponent >= 1.0) {
            dlPowPol = polLambdaExponent * pow(polLambda, polLambdaExponent - 1.0);
            if (polLambdaExponent >= 2.0) {
              d2lPowPol =
                  polLambdaExponent
                      * (polLambdaExponent - 1.0)
                      * pow(polLambda, polLambdaExponent - 2.0);
            }
          }
        }
        // Add the chain rule term due to shrinking the lambda range for the polarization energy.
        dlPowPol *= polLambdaScale;
        d2lPowPol *= (polLambdaScale * polLambdaScale);
      }
    }
  }

  /**
   * For OST mode, we are calculating analytic dU/dL, d2U/dL2 and d2U/dL/dX for the permanent and
   * polarization energy terms.
   * <p>
   * For SCALE mode, the permanent multipoles and polarizabilities are scaled by lambda for
   * alchemical atoms.
   */
  public enum AlchemicalMode {
    OST, SCALE
  }

}
