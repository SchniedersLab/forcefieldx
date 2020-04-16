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

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.BaseScript;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import org.apache.log4j.PropertyConfigurator;

/**
 * Base class for scripts in the Potentials package, providing some key functions.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class PotentialScript extends BaseScript {

  /** An instance of PotentialFunctions passed into the current context. */
  public PotentialsFunctions potentialFunctions;

  /**
   * An active MolecularAssembly passed into the current context or loaded by the Script from a file
   * argument.
   */
  public MolecularAssembly activeAssembly;

  /**
   * Reclaims resources associated with all Potential objects associated with this script.
   *
   * @return If all Potentials had resources reclaimed.
   */
  public boolean destroyPotentials() {
    boolean allSucceeded = true;
    for (Potential potent : getPotentials()) {
      if (potent != null) {
        allSucceeded = allSucceeded && potent.destroy();
      }
    }
    return allSucceeded;
  }

  /**
   * Returns a List of all Potential objects associated with this script. Should be written to
   * tolerate nulls, as many tests run help() and exit without instantiating their Potentials.
   *
   * @return All Potentials. Sometimes empty, never null.
   */
  public List<Potential> getPotentials() {
    List<Potential> plist = new ArrayList<>();
    if (activeAssembly != null && activeAssembly.getPotentialEnergy() != null) {
      plist.add(activeAssembly.getPotentialEnergy());
    }
    return plist;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Execute the BaseScript init method, then load potential functions.
   */
  @Override
  public boolean init() {
    if (!super.init()) {
      return false;
    }

    if (context.hasVariable("functions")) {
      // FFX is running.
      potentialFunctions = (PotentialsFunctions) context.getVariable("functions");
    } else {
      // Potential package is running.
      potentialFunctions = new PotentialsUtils();
      // Turn off log4j.
      Properties properties = new Properties();
      properties.setProperty("log4j.threshold", "OFF");
      properties.setProperty("log4j2.level", "OFF");
      properties.setProperty("org.apache.logging.log4j.level", "OFF");
      PropertyConfigurator.configure(properties);
    }

    activeAssembly = null;
    if (context.hasVariable("active")) {
      activeAssembly = (MolecularAssembly) context.getVariable("active");
    }

    return true;
  }
}
