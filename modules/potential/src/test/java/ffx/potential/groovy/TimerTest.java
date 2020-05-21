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
package ffx.potential.groovy;

import ffx.potential.utils.PotentialTest;
import org.junit.Test;

/** Test the Energy script. */
public class TimerTest extends PotentialTest {

  @Test
  public void testTimer() {
    // Set-up the input arguments for the Timer script.
    String[] args = {"-n", "2", "-v", "-g", "src/main/java/ffx/potential/structures/ubiquitin.xyz"};
    binding.setVariable("args", args);

    // Construct and evaluate the Timer script.
    Timer timer = new Timer(binding).run();
    potentialScript = timer;
  }

  @Test
  public void testTimerHelp() {
    // Set-up the input arguments for the Timer script.
    String[] args = {"-h"};
    binding.setVariable("args", args);

    // Construct and evaluate the Timer script.
    Timer timer = new Timer(binding).run();
    potentialScript = timer;
  }
}
