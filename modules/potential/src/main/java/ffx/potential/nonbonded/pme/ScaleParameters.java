// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import static ffx.potential.parameters.ForceField.ELEC_FORM.PAM;
import static ffx.utilities.PropertyGroup.ElectrostaticsFunctionalForm;
import static java.lang.String.format;

import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ELEC_FORM;
import ffx.utilities.FFXProperty;

import java.util.logging.Logger;

/**
 * Scale factors and masking rules for electrostatics.
 */
public class ScaleParameters {

  private static final Logger logger = Logger.getLogger(ScaleParameters.class.getName());

  private static final double DEFAULT_MPOLE_12_SCALE = 0.0;
  private static final double DEFAULT_MPOLE_13_SCALE = 0.0;
  private static final double DEFAULT_MPOLE_14_SCALE = 1.0;
  private static final double DEFAULT_MPOLE_15_SCALE = 1.0;
  private static final double DEFAULT_CHG_12_SCALE = 0.0;
  private static final double DEFAULT_CHG_13_SCALE = 0.0;
  private static final double DEFAULT_CHG_14_SCALE = 2.0;
  private static final double DEFAULT_CHG_15_SCALE = 1.0;

  /**
   * The interaction energy between 1-2 multipoles is scaled by m12scale.
   */
  @FFXProperty(name = "mpole-12-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to permanent atomic multipole
          electrostatic interactions between 1-2 connected atoms, i.e., atoms that are directly bonded.
          The default value of 0.0 is used, if the mpole-12-scale property is not given
          in either the parameter file or the property file.
          """)
  @FFXProperty(name = "chg-12-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to charge-charge electrostatic
          interactions between 1-2 connected atoms, i.e., atoms that are directly bonded.
          The default value of 0.0 is used, if the chg-12-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double m12scale;

  /**
   * The interaction energy between 1-3 multipoles is scaled by m13scale.
   */
  @FFXProperty(name = "mpole-13-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to permanent atomic multipole
          electrostatic interactions between 1-3 connected atoms, i.e., atoms separated by two covalent bonds.
          The default value of 0.0 is used, if the mpole-13-scale keyword is not given
          in either the parameter file or the property file.
          """)
  @FFXProperty(name = "chg-13-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to charge-charge electrostatic
          interactions between 1-3 connected atoms, i.e., atoms separated by two covalent bonds.
          The default value of 0.0 is used, if the chg-13-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double m13scale;

  /**
   * The interaction energy between 1-4 multipoles is scaled by m14scale.
   */
  @FFXProperty(name = "mpole-14-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to permanent atomic multipole
          electrostatic interactions between 1-4 connected atoms, i.e., atoms separated by three covalent bonds.
          The default value of 1.0 is used, if the mpole-14-scale keyword is not given
          in either the parameter file or the property file.
          """)
  @FFXProperty(name = "chg-14-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to charge-charge electrostatic
          interactions between 1-4 connected atoms, i.e., atoms separated by three covalent bonds.
          The default value of 1.0 is used, if the chg-14-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double m14scale;

  /**
   * The interaction energy between 1-5 multipoles is scaled by m15scale.
   */
  @FFXProperty(name = "mpole-15-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to permanent atomic multipole
          electrostatic interactions between 1-5 connected atoms, i.e., atoms separated by four covalent bonds.
          The default value of 1.0 is used, if the mpole-15-scale keyword is not given
          in either the parameter file or the property file.
          """)
  @FFXProperty(name = "chg-15-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to charge-charge electrostatic
          interactions between 1-5 connected atoms, i.e., atoms separated by four covalent bonds.
          The default value of 1.0 is used, if the chg-15-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double m15scale;

  private static final double DEFAULT_DIRECT_11_SCALE = 0.0;
  private static final double DEFAULT_DIRECT_12_SCALE = 1.0;
  private static final double DEFAULT_DIRECT_13_SCALE = 1.0;
  private static final double DEFAULT_DIRECT_14_SCALE = 1.0;

  /**
   * DIRECT-11-SCALE factor.
   */
  @FFXProperty(name = "direct-11-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to the permanent (direct) field
          due to atoms within a polarization group during an induced dipole calculation,
          i.e., atoms that are in the same polarization group as the atom being polarized.
          The default value of 0.0 is used, if the direct-11-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double d11scale;

  /**
   * The DIRECT_12_SCALE factor is assumed to be 1.0. If this assumption is violated by a keyword,
   * FFX will exit.
   */
  @FFXProperty(name = "direct-12-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to the permanent (direct) field
          due to atoms in 1-2 polarization groups during an induced dipole calculation,
          i.e., atoms that are in polarization groups directly connected to the group containing the atom being polarized.
          The default value of 0.0 is used, if the direct-12-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double d12scale;

  /**
   * The DIRECT_13_SCALE factor is assumed to be 1.0. If this assumption is violated by a keyword,
   * FFX will exit.
   */
  @FFXProperty(name = "direct-13-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to the permanent (direct) field
          due to atoms in 1-3 polarization groups during an induced dipole calculation,
          i.e., atoms that are in polarization groups separated by one group from the group containing the atom being polarized.
          The default value of 0.0 is used, if the direct-13-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double d13scale;

  /**
   * The DIRECT_14_SCALE factor is assumed to be 1.0. If this assumption is violated by a keyword,
   * FFX will exit.
   */
  @FFXProperty(name = "direct-14-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to the permanent (direct) field
          due to atoms in 1-4 polarization groups during an induced dipole calculation,
          i.e., atoms that are in polarization groups separated by two groups from the group containing the atom being polarized.
          The default value of 1.0 is used, if the direct-14-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double d14scale;

  private static final double DEFAULT_POLAR_12_SCALE = 0.0;
  private static final double DEFAULT_POLAR_13_SCALE = 0.0;
  private static final double DEFAULT_POLAR_14_SCALE = 1.0;
  private static final double DEFAULT_POLAR_15_SCALE = 1.0;

  /**
   * The interaction energy between a permanent multipole and polarizable site that are 1-2 is scaled
   * by p12scale.
   */
  @FFXProperty(name = "polar-12-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-2 connected atoms located in different polarization groups.
          The default value of 0.0 is used, if the polar-12-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double p12scale;

  /**
   * The interaction energy between a permanent multipole and polarizable site that are 1-3 is scaled
   * by p13scale.
   */
  @FFXProperty(name = "polar-13-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """ 
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-3 connected atoms located in different polarization groups.
          The default value of 0.0 is used, if the polar-13-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double p13scale;

  /**
   * The interaction energy between a permanent multipole and polarizable site that are 1-4 is scaled
   * by p14scale.
   */
  @FFXProperty(name = "polar-14-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-4 connected atoms located in different polarization groups.
          The default value of 1.0 is used, if the polar-14-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double p14scale;

  /**
   * The interaction energy between a permanent multipole and polarizable site that are 1-5 is scaled
   * by p15scale. Only 1.0 is supported.
   */
  @FFXProperty(name = "polar-15-scale", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-5 connected atoms located in different polarization groups.
          The default value of 1.0 is used, if the polar-15-scale keyword is not given
          in either the parameter file or the property file.
          """)
  public final double p15scale;

  private static final double DEFAULT_POLAR_12_INTRA = 0.0;

  private static final double DEFAULT_POLAR_13_INTRA = 0.0;

  private static final double DEFAULT_POLAR_14_INTRA = 0.5;

  private static final double DEFAULT_POLAR_15_INTRA = 1.0;

  /**
   * An intra-12-scale factor other than 0.0 is not supported and will cause FFX to exit.
   */
  @FFXProperty(name = "polar-12-intra", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-2 connected atoms located in the same polarization group.
          The default value of 0.0 is used, if the polar-12-intra keyword is not given
          in either the parameter file or the property file.
          """)
  public final double intra12Scale;
  /**
   * An intra-13-scale factor other than 0.0 is not supported and will cause FFX to exit.
   */
  @FFXProperty(name = "polar-13-intra", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-3 connected atoms located in the same polarization group.
          The default value of 0.0 is used, if the polar-13-intra keyword is not given
          in either the parameter file or the property file.
          """)
  public final double intra13Scale;
  /**
   * Provides a multiplicative scale factor that is applied to polarization interactions between 1-4
   * connected atoms located in the same polarization group.
   */
  @FFXProperty(name = "polar-14-intra", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "0.5",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-4 connected atoms located in the same polarization group.
          The default value of 0.5 is used, if the polar-14-intra keyword is not given
          in either the parameter file or the property file.
          """)
  public final double intra14Scale;
  /**
   * An intra-15-scale factor other than 1.0 is not supported and will cause FFX to exit.
   */
  @FFXProperty(name = "polar-15-intra", propertyGroup = ElectrostaticsFunctionalForm, defaultValue = "1.0",
      description = """
          Provides a multiplicative scale factor that is applied to polarization interactions
          between 1-5 connected atoms located in the same polarization group.
          The default value of 1.0 is used, if the polar-15-intra keyword is not given
          in either the parameter file or the property file.
          """)
  public final double intra15Scale;

  public ScaleParameters(ELEC_FORM elecForm, ForceField forceField) {
    if (elecForm == PAM) {
      double m12 = forceField.getDouble("MPOLE_12_SCALE", DEFAULT_MPOLE_12_SCALE);
      if (m12 > 1.0) {
        m12 = 1.0 / m12;
      }
      m12scale = m12;
      double m13 = forceField.getDouble("MPOLE_13_SCALE", DEFAULT_MPOLE_13_SCALE);
      if (m13 > 1.0) {
        m13 = 1.0 / m13;
      }
      m13scale = m13;
      double m14 = forceField.getDouble("MPOLE_14_SCALE", DEFAULT_MPOLE_14_SCALE);
      if (m14 > 1.0) {
        m14 = 1.0 / m14;
      }
      m14scale = m14;
      double m15 = forceField.getDouble("MPOLE_15_SCALE", DEFAULT_MPOLE_15_SCALE);
      if (m15 > 1.0) {
        m15 = 1.0 / m15;
      }
      m15scale = m15;
    } else {
      double m12 = forceField.getDouble("CHG_12_SCALE", DEFAULT_CHG_12_SCALE);
      if (m12 > 1.0) {
        m12 = 1.0 / m12;
      }
      m12scale = m12;
      double m13 = forceField.getDouble("CHG_13_SCALE", DEFAULT_CHG_13_SCALE);
      if (m13 > 1.0) {
        m13 = 1.0 / m13;
      }
      m13scale = m13;
      double m14 = forceField.getDouble("CHG_14_SCALE", DEFAULT_CHG_14_SCALE);
      if (m14 > 1.0) {
        m14 = 1.0 / m14;
      }
      m14scale = m14;
      double m15 = forceField.getDouble("CHG_15_SCALE", DEFAULT_CHG_15_SCALE);
      if (m15 > 1.0) {
        m15 = 1.0 / m15;
      }
      m15scale = m15;
    }

    // Multiplicative scale factors applied to polarization interactions between connected
    // atoms located in the same polarization group.
    intra12Scale = forceField.getDouble("POLAR_12_INTRA", DEFAULT_POLAR_12_INTRA);
    if (intra12Scale != 0.0) {
      logger.severe(format(" Unsupported polar-12-intra parameter: %8.6f", intra12Scale));
    }
    intra13Scale = forceField.getDouble("POLAR_13_INTRA", DEFAULT_POLAR_13_INTRA);
    if (intra13Scale != 0.0) {
      logger.severe(format(" Unsupported polar-13-intra parameter: %8.6f", intra13Scale));
    }
    intra14Scale = forceField.getDouble("POLAR_14_INTRA", DEFAULT_POLAR_14_INTRA);
    intra15Scale = forceField.getDouble("POLAR_15_INTRA", DEFAULT_POLAR_15_INTRA);
    if (intra15Scale != 1.0) {
      logger.severe(format(" Unsupported polar-15-intra parameter: %8.6f", intra15Scale));
    }

    // Polarization groups masking rules.
    d11scale = forceField.getDouble("DIRECT_11_SCALE", DEFAULT_DIRECT_11_SCALE);
    d12scale = forceField.getDouble("DIRECT_12_SCALE", DEFAULT_DIRECT_12_SCALE);
    d13scale = forceField.getDouble("DIRECT_13_SCALE", DEFAULT_DIRECT_13_SCALE);
    d14scale = forceField.getDouble("DIRECT_14_SCALE", DEFAULT_DIRECT_14_SCALE);
    if (d12scale != 1.0) {
      logger.severe(format(" Unsupported direct-12-scale parameter: %8.6f", d12scale));
    }
    if (d13scale != 1.0) {
      logger.severe(format(" Unsupported direct-13-scale parameter: %8.6f", d13scale));
    }
    if (d14scale != 1.0) {
      logger.severe(format(" Unsupported direct-14-scale parameter: %8.6f", d14scale));
    }

    // Polarization energy masking rules.
    p12scale = forceField.getDouble("POLAR_12_SCALE", DEFAULT_POLAR_12_SCALE);
    p13scale = forceField.getDouble("POLAR_13_SCALE", DEFAULT_POLAR_13_SCALE);
    p14scale = forceField.getDouble("POLAR_14_SCALE", DEFAULT_POLAR_14_SCALE);
    p15scale = forceField.getDouble("POLAR_15_SCALE", DEFAULT_POLAR_15_SCALE);
    if (p15scale != 1.0) {
      logger.severe(format(" Unsupported polar-15-scale parameter: %8.6f", p15scale));
    }
  }
}
