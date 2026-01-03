// ******************************************************************************
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
// ******************************************************************************
package ffx.xray.scatter;

import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * The XRayScatteringParameters class is a record used to store and access X-ray scattering parameters
 * for specific atoms or elements. It provides methods to retrieve such parameters based
 * on their atomic properties or a unique key. Instances of this class consist of an atom's
 * name, a description, and its associated X-ray form factor data.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public record XRayScatteringParameters(String name,
                                       int atomicNumber,
                                       int charge,
                                       int numberOfGaussians,
                                       double[][] formFactor) {

  private static final Logger logger = Logger.getLogger(XRayScatteringParameters.class.getName());

  private static final HashMap<String, XRayScatteringParameters> formFactorsSuCoppens = new HashMap<>();
  private static final HashMap<String, XRayScatteringParameters> formFactorsCCTBX = new HashMap<>();

  /**
   * Returns a string representation of the X-Ray scattering parameters.
   * The string includes the atom name, description, the first form factor element,
   * and the array of additional form factor parameters.
   *
   * @return A formatted string representation of the object.
   */
  @Override
  public String toString() {
    return format(" %s %d Charge: %d Gaussian Amplitudes (%d): %s",
        name, atomicNumber, charge, numberOfGaussians, Arrays.toString(formFactor[1]));
  }

  /**
   * Retrieves the scattering parameters based on the provided atomic parameters and preference for using 3 Gaussian parameters.
   *
   * @param charge       the formal charge of the atom.
   * @param atomicNumber the atomic number of the element.
   * @param use3G        a boolean flag indicating whether to prefer 3 Gaussian parameters if available.
   * @return The X-ray scattering parameters, or null if no matching data is found.
   */
  public static XRayScatteringParameters getFormFactor(int atomicNumber, int charge, boolean use3G) {
    XRayScatteringParameters parameters = null;
    if (use3G) {
      parameters = getFormFactorCCTBX(atomicNumber, charge);
    }
    // Su and Coppens parameters if use3G is false or CCTBX parameters were not found.
    if (parameters == null) {
      parameters = getFormFactorSuCoppens(atomicNumber, charge);
    }
    // No scattering parameters were found.
    if (parameters == null) {
      String message = format(" Parameters not found for %d with charge %d", atomicNumber, charge);
      logger.severe(message);
    }

    return parameters;
  }

  /**
   * Get the Su and Coppens scattering parameters.
   *
   * @param atomicNumber The atomic number.
   * @param charge       The formal charge of the atom.
   * @return The scattering parameters.
   */
  public static XRayScatteringParameters getFormFactorSuCoppens(int atomicNumber, int charge) {
    String key = Integer.toString(atomicNumber);
    String keyWithCharge = key + "_" + charge;
    if (formFactorsSuCoppens.containsKey(keyWithCharge)) {
      return formFactorsSuCoppens.get(keyWithCharge);
    } else return formFactorsSuCoppens.getOrDefault(key, null);
  }

  /**
   * Get the CCTBX 3 Gaussian scattering parameters.
   *
   * @param atomicNumber The atomic number.
   * @param charge       The formal charge of the atom.
   * @return The scattering parameters.
   */
  public static XRayScatteringParameters getFormFactorCCTBX(int atomicNumber, int charge) {
    String key = Integer.toString(atomicNumber);
    String keyWithCharge = key + "_" + charge;
    if (formFactorsCCTBX.containsKey(keyWithCharge)) {
      return formFactorsCCTBX.get(keyWithCharge);
    } else return formFactorsCCTBX.getOrDefault(key, null);
  }

  static {
    // Load Su & Coppens scattering parameters.
    String[] atoms = XrayParametersSuCoppens.atoms;
    String[] atomsi = XrayParametersSuCoppens.atomsi;
    double[][][] ffactors = XrayParametersSuCoppens.ffactors;
    for (int i = 0; i < atoms.length; i++) {
      String environment = atomsi[i];
      String[] descriptions = environment.split("_");
      int atomicNumber = Integer.parseInt(descriptions[0]);
      int charge = 0;
      if (descriptions.length > 1) {
        charge = Integer.parseInt(descriptions[1]);
      }
      int numberOfGaussians = ffactors[i][1].length;
      XRayScatteringParameters factor = new XRayScatteringParameters(atoms[i],
          atomicNumber, charge, numberOfGaussians, ffactors[i]);
      formFactorsSuCoppens.put(atomsi[i], factor);
    }

    // Load CCTBX 3 Gaussian scattering parameters.
    atoms = XrayParametersCCTBX.atoms;
    atomsi = XrayParametersCCTBX.atomsi;
    ffactors = XrayParametersCCTBX.ffactors;
    for (int i = 0; i < atoms.length; i++) {
      String environment = atomsi[i];
      String[] descriptions = environment.split("_");
      int atomicNumber = Integer.parseInt(descriptions[0]);
      int charge = 0;
      if (descriptions.length > 1) {
        charge = Integer.parseInt(descriptions[1]);
      }
      int numberOfGaussians = ffactors[i][1].length;
      XRayScatteringParameters factor = new XRayScatteringParameters(atoms[i],
          atomicNumber, charge, numberOfGaussians, ffactors[i]);
      formFactorsCCTBX.put(atomsi[i], factor);
    }

  }

}
