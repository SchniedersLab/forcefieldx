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
package ffx.potential.groovy;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import ffx.potential.parameters.BioType;
import ffx.potential.utils.PotentialTest;
import java.util.List;
import org.junit.Test;

/** Test the Biotype script. */
public class BiotypeTest extends PotentialTest {

  @Test
  public void testBiotype() {
    // Set-up the input arguments for the Biotype script.
    String[] args = {"src/main/java/ffx/potential/structures/acetanilide.xyz"};
    binding.setVariable("args", args);

    // Create and evaluate the script.
    Biotype biotype = new Biotype(binding).run();
    potentialScript = biotype;

    // Check the Biotype results.
    List<BioType> bioTypes = biotype.getBioTypes();
    assertNotNull(bioTypes);
    assertEquals(19, bioTypes.size());
    BioType bioType = bioTypes.get(0);
    assertTrue(
        " Check the value of the first Biotype.",
        bioType.toString().trim().equalsIgnoreCase(
            "biotype      1  C     \"ace                    \"    405  C     C     N"));

    // Check that the bioTypes variable is available via the Binding.
    assertEquals(bioTypes, binding.getVariable("bioTypes"));
  }

  @Test
  public void testBiotypeHelp() {
    // Set-up the input arguments for the Biotype script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Create and evaluate the script.
    Biotype biotype = new Biotype(binding).run();
    potentialScript = biotype;
  }
}
