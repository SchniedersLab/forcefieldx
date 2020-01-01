//******************************************************************************
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
//******************************************************************************
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.PropertiesConfiguration;
import org.apache.commons.configuration2.builder.FileBasedConfigurationBuilder;
import org.apache.commons.configuration2.builder.fluent.Parameters;
import org.apache.commons.configuration2.ex.ConfigurationException;

import ffx.potential.parameters.AngleTorsionType;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BaseType;
import ffx.potential.parameters.BioType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldName;
import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.RelativeSolvationType;
import ffx.potential.parameters.SoluteType;
import ffx.potential.parameters.StretchBendType;
import ffx.potential.parameters.StretchTorsionType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;
import ffx.utilities.Keyword;
import static ffx.potential.parameters.AngleType.AngleFunction;
import static ffx.potential.parameters.BondType.BondFunction;

/**
 * The ForceFieldFilter Class is used to parse and store molecular mechanics
 * data from keyword/property and parameter (*.PRM) files.
 * <p>
 * Alternatively, an Apache Commons "Configuration" instance can be parsed.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ForceFieldFilter {

    private static final Logger logger = Logger.getLogger(ForceFieldFilter.class.getName());

    /**
     * Default internal force field.
     */
    private static final ForceFieldName DEFAULT_FORCE_FIELD = ForceFieldName.AMOEBA_BIO_2018;
    /**
     * An external force field parameter file.
     */
    private final File forceFieldFile;
    /**
     * The ForceField instance that will be returned.
     */
    private ForceField forceField;
    /**
     * A CompositeConfiguration instance to search for force field types.
     */
    private final CompositeConfiguration properties;

    /**
     * <p>
     * Constructor for ForceFieldFilter.</p>
     *
     * @param properties a
     *                   {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    public ForceFieldFilter(CompositeConfiguration properties) {
        this.properties = properties;
        if (properties.containsKey("parameters")) {
            String fileName = properties.getString("parameters");
            if (properties.containsKey("propertyFile")) {
                String propertyName = properties.getString("propertyFile");
                logger.info(" Property File: " + propertyName);
                File propertyFile = new File(propertyName);
                forceFieldFile = parseParameterLocation(fileName, propertyFile);
            } else {
                forceFieldFile = parseParameterLocation(fileName, null);
            }
        } else {
            forceFieldFile = null;
        }
        forceField = new ForceField(properties);
    }

    /**
     * <p>
     * parse</p>
     *
     * @return a {@link ffx.potential.parameters.ForceField} object.
     */
    public ForceField parse() {
        try {
            // Try to parse an external (ie. not in the FFX jar) parameter file.
            if (forceFieldFile != null) {
                File fileToOpen = forceFieldFile;
                if (!fileToOpen.exists()) {
                    fileToOpen = new File(forceFieldFile.getAbsolutePath() + ".prm");
                    if (!fileToOpen.exists()) {
                        logger.log(Level.INFO, " {0} does not exist.", forceFieldFile);
                        return null;
                    }
                }
                if (!fileToOpen.canRead()) {
                    logger.log(Level.INFO, " {0} can not be read.", fileToOpen);
                    return null;
                }
                parse(new FileInputStream(fileToOpen));
            } else {
                // Parse an internal parameter file and add it to the composite configuration.
                String defaultFFstring = DEFAULT_FORCE_FIELD.toString().toUpperCase().replaceAll("_", "-");
                String forceFieldString = properties.getString("forcefield", defaultFFstring);
                ForceFieldName ff;
                try {
                    ff = ForceField.ForceFieldName.valueOf(forceFieldString.toUpperCase().replace('-', '_'));
                } catch (Exception e) {
                    logger.info(format(" The forcefield property %s was not recognized.\n", forceFieldString));
                    ff = DEFAULT_FORCE_FIELD;
                }
                URL url = ForceField.getForceFieldURL(ff);
                if (url != null) {
                    forceField.forceFieldURL = url;
                    try {
                        FileBasedConfigurationBuilder<PropertiesConfiguration> builder =
                                new FileBasedConfigurationBuilder<>(PropertiesConfiguration.class)
                                        .configure(
                                                new Parameters().properties().setURL(url).setThrowExceptionOnMissing(true)
                                                        //.setListDelimiterHandler(new DefaultListDelimiterHandler(','))
                                                        .setIncludesAllowed(false));
                        PropertiesConfiguration forcefieldConfiguration = builder.getConfiguration();
                        String name = ForceField.toPropertyForm(ff.toString());
                        forcefieldConfiguration.setHeader("Internal force field (" + name + ").");
                        properties.addConfiguration(forcefieldConfiguration);
                    } catch (ConfigurationException e) {
                        e.printStackTrace();
                        logger.warning(e.toString());
                    }
                }
            }

            // Overwrite parameters of the forceFieldFile with those from the CompositeConfiguration.
            if (properties != null) {
                parse(properties);
            }
        } catch (FileNotFoundException e) {
            String message = "Exception parsing force field.";
            logger.log(Level.WARNING, message, e);
        }

        String forceFieldName = forceField.getString("forcefield", "none");
        if (forceFieldName != null) {
            forceFieldName = forceFieldName.toUpperCase();
            if (forceFieldName.contains("AMBER") || forceFieldName.contains("OPLS") ||
                forceFieldName.contains("CHARMM")) {
                forceField.setBondFunction(BondFunction.HARMONIC);
                forceField.setAngleFunction(AngleFunction.HARMONIC);
            }
        }

        return forceField;
    }

    private void parse(CompositeConfiguration properties) {
        try {
            int numConfigs = properties.getNumberOfConfigurations();

            /*
              Loop over the configurations starting with lowest precedence.
              This way higher precedence entries will overwrite lower
              precedence entries within the ForceField instance.
             */
            for (int n = numConfigs - 1; n >= 0; n--) {
                Configuration config = properties.getConfiguration(n);
                if (config instanceof PropertiesConfiguration) {
                    PropertiesConfiguration propertiesConfiguration = (PropertiesConfiguration) config;
                    if (propertiesConfiguration.getHeader() != null) {
                        logger.info(" Parsing: " + propertiesConfiguration.getHeader());
                    }
                }

                Iterator i = config.getKeys();
                while (i.hasNext()) {
                    String key = (String) i.next();

                    // If the key is not recognized as a force field keyword, continue to the next key.
                    if (!ForceField.isForceFieldKeyword(key)) {
                        continue;
                    }

                    String[] list = config.getStringArray(key);
                    for (String s : list) {
                        // Add back the key to the input line.
                        s = key + " " + s;

                        // Split the line on the pound symbol to remove comments.
                        String input = s.split("#+")[0];
                        String[] tokens = input.trim().split(" +");

                        // Parse force field types.
                        ForceFieldType type;
                        try {
                            type = ForceFieldType.valueOf(key.toUpperCase());
                        } catch (Exception e) {
                            break;
                        }

                        BaseType baseType = null;
                        switch (type) {
                            case ATOM:
                                baseType = AtomType.parse(input, tokens);
                                break;
                            case ANGTORS:
                                baseType = AngleTorsionType.parse(input, tokens);
                                break;
                            case ANGLE:
                                baseType = AngleType.parse(input, tokens);
                                break;
                            case ANGLEP:
                                baseType = AngleType.parseInPlane(input, tokens);
                                break;
                            case BIOTYPE:
                                baseType = BioType.parse(input, tokens);
                                break;
                            case BOND:
                                baseType = BondType.parse(input, tokens);
                                break;
                            case CHARGE:
                                baseType = MultipoleType.parseChargeType(input, tokens);
                                break;
                            case MULTIPOLE:
                                baseType = MultipoleType.parse(input, tokens);
                                break;
                            case OPBEND:
                                baseType = OutOfPlaneBendType.parse(input, tokens);
                                break;
                            case STRBND:
                                baseType = StretchBendType.parse(input, tokens);
                                break;
                            case PITORS:
                                baseType = PiTorsionType.parse(input, tokens);
                                break;
                            case IMPTORS:
                                baseType = ImproperTorsionType.parse(input, tokens);
                                break;
                            case TORSION:
                                baseType = TorsionType.parse(input, tokens);
                                break;
                            case IMPROPER:
                                baseType = TorsionType.parseImproper(input, tokens);
                                break;
                            case STRTORS:
                                baseType = StretchTorsionType.parse(input, tokens);
                                break;
                            case TORTORS:
                                baseType = TorsionTorsionType.parse(input, tokens);
                                break;
                            case UREYBRAD:
                                baseType = UreyBradleyType.parse(input, tokens);
                                break;
                            case VDW:
                                baseType = VDWType.parse(input, tokens);
                                break;
                            case VDW14:
                                baseType = VDWType.parseVDW14(input, tokens);
                                break;
                            case POLARIZE:
                                baseType = PolarizeType.parse(input, tokens);
                                break;
                            case RELATIVESOLV:
                                baseType = RelativeSolvationType.parse(input, tokens);
                                break;
                            case SOLUTE:
                                baseType = SoluteType.parse(input, tokens);
                                break;
                            default:
                                logger.log(Level.WARNING, "ForceField type recognized, but not stored:{0}", type);
                        }
                        if (baseType != null) {
                            forceField.addForceFieldType(baseType);
                        }
                    }
                }
            }
        } catch (Exception e) {
            String message = "Exception parsing force field.";
            logger.log(Level.WARNING, message, e);
        }

        logger.info("");
    }

    private void parse(InputStream stream) {
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(stream));
            while (br.ready()) {
                String input = br.readLine();
                parse(input, br);
            }
        } catch (IOException e) {
            String message = "Error parsing force field parameters.\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private void parse(String input, BufferedReader br) {

        // Split the line on the pound symbol to remove comments.
        String[] inputs = input.split("#");

        if (inputs.length < 1) {
            return;
        }

        input = inputs[0].split("#")[0];

        // Split the line into tokens between instances of 1 or more spaces.
        String[] tokens = input.trim().split(" +");

        // Check for the case of no tokens or a Keyword.
        // if (parseKeyword(tokens)) {
        //     return;
        // }

        if (tokens[0].equalsIgnoreCase("")) {
            return;
        }

        try {
            ForceFieldType type = ForceFieldType.valueOf(tokens[0].toUpperCase());
            BaseType baseType = null;
            switch (type) {
                case ATOM:
                    baseType = AtomType.parse(input, tokens);
                    break;
                case ANGLE:
                    baseType = AngleType.parse(input, tokens);
                    break;
                case ANGLEP:
                    baseType = AngleType.parseInPlane(input, tokens);
                    break;
                case ANGTORS:
                    baseType = AngleTorsionType.parse(input, tokens);
                    break;
                case BIOTYPE:
                    baseType = BioType.parse(input, tokens);
                    break;
                case BOND:
                    baseType = BondType.parse(input, tokens);
                    break;
                case CHARGE:
                    baseType = MultipoleType.parseChargeType(input, tokens);
                    break;
                case MULTIPOLE:
                    baseType = MultipoleType.parse(input, tokens, br);
                    break;
                case OPBEND:
                    baseType = OutOfPlaneBendType.parse(input, tokens);
                    break;
                case STRBND:
                    baseType = StretchBendType.parse(input, tokens);
                    break;
                case PITORS:
                    baseType = PiTorsionType.parse(input, tokens);
                    break;
                case IMPTORS:
                    baseType = ImproperTorsionType.parse(input, tokens);
                    break;
                case STRTORS:
                    baseType = StretchTorsionType.parse(input, tokens);
                    break;
                case TORSION:
                    baseType = TorsionType.parse(input, tokens);
                    break;
                case IMPROPER:
                    baseType = TorsionType.parseImproper(input, tokens);
                    break;
                case TORTORS:
                    baseType = TorsionTorsionType.parse(input, tokens, br);
                    break;
                case UREYBRAD:
                    baseType = UreyBradleyType.parse(input, tokens);
                    break;
                case VDW:
                    baseType = VDWType.parse(input, tokens);
                    break;
                case VDW14:
                    baseType = VDWType.parseVDW14(input, tokens);
                    break;
                case POLARIZE:
                    baseType = PolarizeType.parse(input, tokens);
                    break;
                case RELATIVESOLV:
                    baseType = RelativeSolvationType.parse(input, tokens);
                    break;
                case SOLUTE:
                    baseType = SoluteType.parse(input, tokens);
                    break;
                default:
                    logger.log(Level.WARNING, "ForceField type recognized, but not stored:{0}", type);
            }
            if (baseType != null) {
                forceField.addForceFieldType(baseType);
            }
            return;
        } catch (Exception e) {
            // Note -- this serves to skip blank lines in *.patch files but also hide an actual bug's exception
            // String message = "Exception parsing force field parametesr.\n";
            // logger.log(Level.WARNING, message, e);
        }

        // Otherwise -- add this entry as a property.
        try {
            String key = tokens[0];
            String value = input.replaceFirst(tokens[0], "").trim();
            forceField.addProperty(key, value);
        } catch (Exception e) {
            logger.info(" Ignored line: " + input);
        }

    }

    /**
     * <p>
     * parseParameterLocation</p>
     *
     * @param parameterLocation a {@link java.lang.String} object.
     * @param keyFile           a {@link java.io.File} object.
     * @return a {@link java.io.File} object.
     */
    private static File parseParameterLocation(String parameterLocation, File keyFile) {
        File parameterFile = null;
        if (parameterLocation != null && !parameterLocation.equalsIgnoreCase("NONE")) {
            // Remove quotes
            parameterLocation = parameterLocation.replaceAll("\"", "");
            // Append the suffix if necessary
            parameterFile = new File(parameterLocation);
            // If the location is not absolute, check if it is relative
            // to the key file location.
            if (!parameterFile.exists() && keyFile != null) {
                parameterFile = new File(keyFile.getParent() + File.separator + parameterLocation);
            }
        }
        return parameterFile;
    }

    /**
     * Parse a Force Field parameter file and echo the results with slashes.
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        if (args == null || args.length < 1) {
            System.out.println("Usage: ForceFieldFilter <file.prm>");
            System.exit(-1);
        }
        CompositeConfiguration properties = Keyword.loadProperties(null);
        properties.setProperty("parameters", args[0]);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();
        if (forceField != null) {
            forceField.print();
        }
    }
}
