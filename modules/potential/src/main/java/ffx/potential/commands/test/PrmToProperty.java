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
package ffx.potential.commands.test;

import ffx.potential.cli.PotentialScript;
import ffx.potential.parameters.ChargeType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.VDWType;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.utilities.Keyword;
import org.apache.commons.configuration2.CompositeConfiguration;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;
import groovy.lang.Binding;

import java.util.List;
import java.util.Map;

/**
 * The PrmToProperty script converts a TINKER *.prm file to Java properties.
 * Usage: ffxc test.PrmToProperty <filename>
 */
@Command(description = "PrmToProperty converts a TINKER *.prm file to Java properties.",
        name = "test.PrmToProperty")
public class PrmToProperty extends PotentialScript {

    /**
     * -t or --tinker Remove line continuation characters from multi-line force field types (i.e. Tinker prm format).
     */
    @Option(names = {"-t", "--tinker"}, paramLabel = "false", defaultValue = "false",
            description = "Remove line continuation characters from multi-line force field types (i.e. Tinker prm format).")
    private boolean tinker = false;

    /**
     * -c or --useCharges Output charge keywords instead of multipoles.
     */
    @Option(names = {"-c", "--useCharges"}, paramLabel = "false", defaultValue = "false",
            description = "Output charge keywords instead of multipole keywords.")
    private boolean useCharges = false;

    /**
     * -z or --zeroVDW Add a zero VDW class with zero parameters.
     */
    @Option(names = {"-z", "--zeroVDW"}, paramLabel = "false", defaultValue = "false",
            description = "Add a zero VDW class with zero parameters.")
    private boolean zeroVDW = false;

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = "TINKER *.prm file(s).")
    private List<String> filenames = null;

    /**
     * PrmToProperty Constructor.
     */
    public PrmToProperty() {
        super();
    }

    /**
     * PrmToProperty Constructor.
     * @param binding Groovy Binding to use.
     */
    public PrmToProperty(Binding binding) {
        super(binding);
    }

    /**
     * PrmToProperty constructor that sets the command line arguments.
     * @param args Command line arguments.
     */
    public PrmToProperty(String[] args) {
        super(args);
    }

    /**
     * Execute the script.
     */
    @Override
    public PrmToProperty run() {

        if (!init()) {
            return this;
        }

        if (filenames == null || filenames.size() < 1) {
            logger.info(helpString());
            return this;
        }

        // Read in the command line file.
        String prmName = filenames.get(0);
        CompositeConfiguration properties = Keyword.loadProperties(null);
        properties.setProperty("parameters", prmName);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);

        ForceField forceField = forceFieldFilter.parse();

        for (int i = 1; i < filenames.size(); i++) {
            String filename = filenames.get(i);
            properties = Keyword.loadProperties(null);
            properties.setProperty("parameters", filename);
            forceFieldFilter = new ForceFieldFilter(properties);
            ForceField forceField2 = forceFieldFilter.parse();
            forceField.append(forceField2);
        }

        if (useCharges) {
            // Convert multipole types to charge types.
            @SuppressWarnings("unchecked")
            Map<String, MultipoleType> multipoleTypes = (Map<String, MultipoleType>) (Map<?, ?>)
                    forceField.getTypes(ForceField.ForceFieldType.MULTIPOLE);
            for (String key : multipoleTypes.keySet()) {
                MultipoleType multipoleType = multipoleTypes.get(key);
                ChargeType chargeType = new ChargeType(multipoleType.frameAtomTypes[0], multipoleType.getCharge());
                forceField.addForceFieldType(chargeType);
            }
            // Clear multipole types.
            forceField.clearForceFieldType(ForceField.ForceFieldType.MULTIPOLE);
        }

        if (zeroVDW) {
            // Add a zero VDW class with zero parameters.
            forceField.addForceFieldType(new VDWType(0, 0.0, 0.0, -1.0));
        }

        if (forceField != null) {
            StringBuffer sb = forceField.toStringBuffer();
            if (tinker) {
                logger.info(sb.toString().replace('\\', ' '));
            } else {
                logger.info(sb.toString());
            }
        }

        return this;
    }
}
