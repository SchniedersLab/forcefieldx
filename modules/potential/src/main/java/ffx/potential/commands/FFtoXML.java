//******************************************************************************
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
//******************************************************************************
package ffx.potential.commands;

import ffx.potential.cli.PotentialScript;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.OpenMMXmlFilter;
import ffx.utilities.Keyword;
import groovy.lang.Binding;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;

/**
 * The FFtoXML command saves a force field as an XML file usable by OpenMM.
 *
 * Usage:
 *   ffxc FFtoXML &lt;filename&gt;
 */
@Command(name = "FFtoXML", description = " Write a force field as an XML file.")
public class FFtoXML extends PotentialScript {

  /** The final argument is a PRM/KEY file. */
  @Parameters(arity = "1", paramLabel = "file",
      description = "The force field file.")
  private String filename = null;

  public FFtoXML() { super(); }
  public FFtoXML(Binding binding) { super(binding); }
  public FFtoXML(String[] args) { super(args); }

  @Override
  public FFtoXML run() {
    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    // Create PRM/KEY file's Force Field
    CompositeConfiguration props = Keyword.loadProperties(null);
    props.setProperty("parameters", filename);
    ForceFieldFilter forceFieldFilter = new ForceFieldFilter(props);
    ForceField forceField = forceFieldFilter.parse();

    // Use the current base directory, or update if necessary based on the given filename.
    String saveName = getBaseDirString(filename) + removeExtension(getName(filename));

    OpenMMXmlFilter xmlFilter = new OpenMMXmlFilter(forceField, saveName);
    try {
      xmlFilter.toXML();
    } catch (Exception e) {
      logger.severe(" Error writing XML: " + e.getMessage());
    }

    return this;
  }
}
