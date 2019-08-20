//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
import java.util.Arrays;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.PropertiesConfiguration;
import org.apache.commons.configuration2.builder.FileBasedConfigurationBuilder;
import org.apache.commons.configuration2.builder.fluent.Parameters;
import org.apache.commons.configuration2.ex.ConfigurationException;
import static org.apache.commons.math3.util.FastMath.abs;

import ffx.potential.parameters.AngleTorsionType;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BioType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import ffx.potential.parameters.ForceField.ForceFieldName;
import ffx.potential.parameters.ForceField.ForceFieldString;
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
import static ffx.potential.parameters.ForceField.toEnumForm;

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
    private static final ForceFieldName DEFAULT_FORCE_FIELD = ForceFieldName.AMOEBA_BIO_2018;

    private ForceField forceField;

    private final CompositeConfiguration properties;

    private final File forceFieldFile;

    private boolean convertRadiusToDiameter = false;
    private boolean convertSigmaToRMin = false;
    private double torsionScale = -1.0;
    private double improperTorsionScale = -1.0;

    /**
     * <p>Setter for the field <code>convertRadiusToDiameter</code>.</p>
     *
     * @param convertRadiusToDiameter a boolean.
     */
    public void setConvertRadiusToDiameter(boolean convertRadiusToDiameter) {
        this.convertRadiusToDiameter = convertRadiusToDiameter;
    }

    /**
     * <p>Setter for the field <code>convertSigmaToRMin</code>.</p>
     *
     * @param convertSigmaToRMin a boolean.
     */
    public void setConvertSigmaToRMin(boolean convertSigmaToRMin) {
        this.convertSigmaToRMin = convertSigmaToRMin;
    }

    /**
     * <p>Setter for the field <code>torsionScale</code>.</p>
     *
     * @param torsionScale a double.
     */
    public void setTorsionScale(double torsionScale) {
        this.torsionScale = torsionScale;
    }

    /**
     * <p>Setter for the field <code>improperTorsionScale</code>.</p>
     *
     * @param improperTorsionScale a double.
     */
    public void setImproperTorsionScale(double improperTorsionScale) {
        this.improperTorsionScale = improperTorsionScale;
    }

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

        return forceField;
    }

    private void parse(CompositeConfiguration properties) {
        try {
            int numConfigs = properties.getNumberOfConfigurations();

            if (numConfigs > 0) {
                logger.info(" Parsing properties from: ");
            }

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
                        logger.info("  " + propertiesConfiguration.getHeader());
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

                        switch (type) {
                            case ATOM:
                                parseAtomType(input, tokens);
                                break;
                            case ANGTORS:
                                parseAngleTorsion(input, tokens);
                                break;
                            case ANGLE:
                                parseAngleType(input, tokens);
                                break;
                            case ANGLEP:
                                parseAngleTypeInPlane(input, tokens);
                                break;
                            case BIOTYPE:
                                parseBioType(input, tokens);
                                break;
                            case BOND:
                                parseBondType(input, tokens);
                                break;
                            case CHARGE:
                                parseChargeType(input, tokens);
                                break;
                            case MULTIPOLE:
                                parseMultipoleType(input, tokens);
                                break;
                            case OPBEND:
                                parseOutOfPlaneBendType(input, tokens);
                                break;
                            case STRBND:
                                parseStretchBendType(input, tokens);
                                break;
                            case PITORS:
                                parsePiTorsionType(input, tokens);
                                break;
                            case IMPTORS:
                                parseImpTorsType(input, tokens);
                                break;
                            case TORSION:
                                parseTorsion(input, tokens);
                                break;
                            case IMPROPER:
                                parseImproperTorsion(input, tokens);
                                break;
                            case STRTORS:
                                parseStretchTorsion(input, tokens);
                                break;
                            case TORTORS:
                                parseTorsionTorsionType(input, tokens);
                                break;
                            case UREYBRAD:
                                parseUreyBradleyType(input, tokens);
                                break;
                            case VDW:
                                parseVDWType(input, tokens);
                                break;
                            case VDW14:
                                parseVDW14Type(input, tokens);
                                break;
                            case POLARIZE:
                                parsePolarizeType(input, tokens);
                                break;
                            case RELATIVESOLV:
                                parseRelativeSolvation(input, tokens);
                                break;
                            case SOLUTE:
                                parseSolute(input, tokens);
                                break;
                            default:
                                logger.log(Level.WARNING, "ForceField type recognized, but not stored:{0}", type);
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
        if (parseKeyword(tokens)) {
            return;
        }

        try {
            ForceFieldType type = ForceFieldType.valueOf(tokens[0].toUpperCase());
            switch (type) {
                case ATOM:
                    parseAtomType(input, tokens);
                    break;
                case ANGLE:
                    parseAngleType(input, tokens);
                    break;
                case ANGLEP:
                    parseAngleTypeInPlane(input, tokens);
                    break;
                case ANGTORS:
                    parseAngleTorsion(input, tokens);
                    break;
                case BIOTYPE:
                    parseBioType(input, tokens);
                    break;
                case BOND:
                    parseBondType(input, tokens);
                    break;
                case CHARGE:
                    parseChargeType(input, tokens);
                    break;
                case MULTIPOLE:
                    parseMultipoleType(input, tokens, br);
                    break;
                case OPBEND:
                    parseOutOfPlaneBendType(input, tokens);
                    break;
                case STRBND:
                    parseStretchBendType(input, tokens);
                    break;
                case PITORS:
                    parsePiTorsionType(input, tokens);
                    break;
                case IMPTORS:
                    parseImpTorsType(input, tokens);
                    break;
                case STRTORS:
                    parseStretchTorsion(input, tokens);
                    break;
                case TORSION:
                    parseTorsion(input, tokens);
                    break;
                case IMPROPER:
                    parseImproperTorsion(input, tokens);
                    break;
                case TORTORS:
                    parseTorsionTorsionType(input, tokens, br);
                    break;
                case UREYBRAD:
                    parseUreyBradleyType(input, tokens);
                    break;
                case VDW:
                    parseVDWType(input, tokens);
                    break;
                case VDW14:
                    parseVDW14Type(input, tokens);
                    break;
                case POLARIZE:
                    parsePolarizeType(input, tokens);
                    break;
                case RELATIVESOLV:
                    parseRelativeSolvation(input, tokens);
                    break;
                case SOLUTE:
                    parseSolute(input, tokens);
                    break;
                default:
                    logger.log(Level.WARNING, "ForceField type recognized, but not stored:{0}", type);
            }
        } catch (Exception e) {
            // Note -- this serves to skip blank lines in *.patch files but also hide an actual bug's exception
            // String message = "Exception parsing force field parametesr.\n";
            // logger.log(Level.WARNING, message, e);
        }
    }

    private boolean parseKeyword(String[] tokens) {
        String keyword = toEnumForm(tokens[0]);
        try {
            // Parse Keywords with a String value.
            ForceFieldString ffString = ForceFieldString.valueOf(keyword);
            int len = tokens.length;
            if (len == 1) {
                forceField.addForceFieldString(ffString, null);
            } else if (len == 2) {
                forceField.addForceFieldString(ffString, tokens[1]);
            } else {
                StringBuilder stringBuilder = new StringBuilder(tokens[1]);
                for (int i = 2; i < len; i++) {
                    stringBuilder.append(" ");
                    stringBuilder.append(tokens[i]);
                }
                forceField.addForceFieldString(ffString, stringBuilder.toString());
            }
        } catch (Exception e) {
            try {
                // Parse Keywords with a Double value.
                ForceFieldDouble ffDouble = ForceFieldDouble.valueOf(keyword);
                double value = parseDouble(tokens[1]);
                forceField.addForceFieldDouble(ffDouble, value);
            } catch (Exception e2) {
                try {
                    // Parse Keywords with an Integer value.
                    ForceFieldInteger ffInteger = ForceFieldInteger.valueOf(keyword);
                    int value = parseInt(tokens[1]);
                    forceField.addForceFieldInteger(ffInteger, value);
                } catch (Exception e3) {
                    try {
                        // Parse Keywords with a Boolean value.
                        ForceFieldBoolean ffBoolean = ForceFieldBoolean.valueOf(keyword);
                        boolean value = true;
                        if (tokens.length > 1 && tokens[0].toUpperCase().endsWith("TERM")) {

                            // Handle the token "ONLY" specially to shut off all other terms.
                            if (tokens[1].equalsIgnoreCase("ONLY")) {
                                for (ForceFieldBoolean term : ForceFieldBoolean.values()) {
                                    if (term.toString().toUpperCase().endsWith("TERM")) {
                                        forceField.addForceFieldBoolean(term, false);
                                    }
                                }
                            } else if (tokens[1].equalsIgnoreCase("NONE")) {
                                // Legacy support for the "NONE" token.
                                value = false;
                            } else {
                                value = Boolean.parseBoolean(tokens[1]);
                            }
                        }
                        forceField.addForceFieldBoolean(ffBoolean, value);
                        forceField.log(keyword);
                    } catch (Exception e4) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    private void parseRelativeSolvation(String input, String[] tokens) {
        if (tokens.length < 3) {
            logger.log(Level.WARNING, "Invalid RELATIVE_SOLVATION type:\n{0}", input);
            return;
        }
        String resName = tokens[1];
        try {
            double relSolvValue = parseDouble(tokens[2]);
            RelativeSolvationType rtype = new RelativeSolvationType(resName, relSolvValue);
            forceField.addForceFieldType(rtype);
        } catch (NumberFormatException ex) {
            String message = "Exception parsing RELATIVE_SOLVATION type:\n" + input + "\n";
            logger.log(Level.SEVERE, message, ex);
        }

    }

    private AngleType parseAngleType(String input, String[] tokens) {
        if (tokens.length < 6) {
            logger.log(Level.WARNING, "Invalid ANGLE type:\n{0}", input);
        } else {
            int[] atomClasses = new int[3];
            double forceConstant = 0.0;
            int angles = 0;
            double[] bondAngle = null;
            try {
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                forceConstant = parseDouble(tokens[4]);
                angles = tokens.length - 5;
                bondAngle = new double[angles];
                for (int i = 0; i < angles; i++) {
                    bondAngle[i] = parseDouble(tokens[5 + i]);
                }
            } catch (NumberFormatException e) {
                String message = "Exception parsing ANGLE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
            double[] newBondAngle = new double[angles];
            arraycopy(bondAngle, 0, newBondAngle, 0, angles);

            String forceFieldName = forceField.toString().toUpperCase();
            AngleType angleType;
            if (forceFieldName.contains("OPLS") || forceFieldName.contains("CHARMM")
                    || forceFieldName.contains("AMBER")) {
                angleType = new AngleType(atomClasses, forceConstant,
                        newBondAngle, AngleType.AngleFunction.HARMONIC);
            } else {
                angleType = new AngleType(atomClasses, forceConstant,
                        newBondAngle, AngleType.AngleFunction.SEXTIC);
            }
            forceField.addForceFieldType(angleType);
            return angleType;
        }
        return null;
    }

    private AngleType parseAngleTypeInPlane(String input, String[] tokens) {
        if (tokens.length < 6) {
            logger.log(Level.WARNING, "Invalid ANGLEP type:\n{0}", input);
        } else {
            int[] atomClasses = new int[3];
            double forceConstant = 0.0;
            int angles = 0;
            double[] bondAngle = null;
            try {
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                forceConstant = parseDouble(tokens[4]);
                angles = tokens.length - 5;
                bondAngle = new double[angles];
                for (int i = 0; i < angles; i++) {
                    bondAngle[i] = parseDouble(tokens[5 + i]);
                }
            } catch (NumberFormatException e) {
                String message = "Exception parsing ANGLEP type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
            double[] newBondAngle = new double[angles];
            arraycopy(bondAngle, 0, newBondAngle, 0, angles);

            String forceFieldName = forceField.toString().toUpperCase();
            AngleType angleType;
            if (forceFieldName.contains("OPLS") || forceFieldName.contains("CHARMM")
                    || forceFieldName.contains("AMBER")) {
                angleType = new AngleType(atomClasses, forceConstant,
                        newBondAngle, AngleType.AngleFunction.HARMONIC, AngleType.AngleMode.IN_PLANE);
            } else {
                angleType = new AngleType(atomClasses, forceConstant,
                        newBondAngle, AngleType.AngleFunction.SEXTIC, AngleType.AngleMode.IN_PLANE);
            }
            forceField.addForceFieldType(angleType);
            return angleType;
        }
        return null;
    }

    private AtomType parseAtomType(String input, String[] tokens) {
        if (tokens.length < 7) {
            logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
        } else {
            try {
                int index = 1;
                // Atom Type
                int type = parseInt(tokens[index++]);
                // Atom Class
                int atomClass;
                // The following try/catch is a nasty hack to check for one of the
                // the following two cases:
                //
                // NUMBER TYPE CLASS IDENTIFIER ... (example is OPLSAA)
                // vs.
                // NUMBER TYPE IDENTIFIER ... (example is OPLSUA)
                //
                // If there is no atom class, a harmless exception will be caught
                // and the atomClass field will remain equal to null.
                try {
                    atomClass = parseInt(tokens[index]);
                    // If the parseInt succeeds, this force field has atom classes.
                    index++;
                } catch (NumberFormatException e) {
                    // Some force fields do not use atom classes.
                    atomClass = -1;
                }
                // Name
                String name = tokens[index].intern();
                // The "environment" string may contain spaces,
                // and is therefore surrounded in quotes located at "first" and
                // "last".
                int first = input.indexOf("\"");
                int last = input.lastIndexOf("\"");
                if (first >= last) {
                    logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
                    return null;
                }
                // Environment
                String environment = input.substring(first, last + 1).intern();
                // Shrink the tokens array to only include entries
                // after the environment field.
                tokens = input.substring(last + 1).trim().split(" +");
                index = 0;
                // Atomic Number
                int atomicNumber = parseInt(tokens[index++]);
                // Atomic Mass
                double mass = parseDouble(tokens[index++]);
                // Hybridization
                int hybridization = parseInt(tokens[index]);
                AtomType atomType = new AtomType(type, atomClass, name,
                        environment, atomicNumber, mass, hybridization);
                forceField.addForceFieldType(atomType);

                return atomType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing CHARGE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);

            }
        }
        return null;
    }

    private BioType parseBioType(String input, String[] tokens) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid BIOTYPE type:\n{0}", input);
        } else {
            try {
                int index = parseInt(tokens[1]);
                String atomName = tokens[2];
                // The "residue" string may contain spaces,
                // and is therefore surrounded in quotes located at "first" and
                // "last".
                int first = input.indexOf("\"");
                int last = input.lastIndexOf("\"");
                if (first >= last) {
                    logger.log(Level.WARNING, "Invalid BIOTYPE type:\n{0}", input);
                    return null;
                }
                // Environment
                String moleculeName = input.substring(first, last + 1).intern();
                // Shrink the tokens array to only include entries
                // after the environment field.
                tokens = input.substring(last + 1).trim().split(" +");
                int atomType = parseInt(tokens[0]);
                int bondCount = tokens.length - 1;
                String[] bonds = null;
                if (bondCount > 0) {
                    bonds = new String[bondCount];
                    arraycopy(tokens, 1, bonds, 0, bondCount);
                }
                BioType bioType = new BioType(index, atomName, moleculeName, atomType, bonds);
                forceField.addForceFieldType(bioType);
                return bioType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing BIOTYPE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private BondType parseBondType(String input, String[] tokens) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid BOND type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[2];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                double forceConstant = parseDouble(tokens[3]);
                double distance = parseDouble(tokens[4]);
                String forceFieldName = forceField.toString().toUpperCase();
                BondType bondType;
                if (forceFieldName.contains("OPLS") || forceFieldName.contains("CHARMM")
                        || forceFieldName.contains("AMBER")) {
                    bondType = new BondType(atomClasses, forceConstant,
                            distance, BondType.BondFunction.HARMONIC);
                } else {
                    bondType = new BondType(atomClasses, forceConstant,
                            distance, BondType.BondFunction.QUARTIC);
                }
                forceField.addForceFieldType(bondType);
                return bondType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing BOND type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * Map charge parameters to a Multipole instance.
     *
     * @param input  Input string.
     * @param tokens Input string tokens.
     */
    private MultipoleType parseChargeType(String input, String[] tokens) {
        if (tokens.length < 3) {
            logger.log(Level.WARNING, "Invalid CHARGE type:\n{0}", input);
        } else {
            try {
                int[] atomTypes = new int[]{parseInt(tokens[1]), 0, 0};
                double partialCharge = parseDouble(tokens[2]);
                double[] dipole = new double[3];
                double[][] quadrupole = new double[3][3];
                MultipoleType.MultipoleFrameDefinition frameDefinition
                        = MultipoleType.MultipoleFrameDefinition.ZTHENX;
                MultipoleType multipoleType = new MultipoleType(partialCharge, dipole,
                        quadrupole, atomTypes, frameDefinition, true);
                forceField.addForceFieldType(multipoleType);
                return multipoleType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing CHARGE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private MultipoleType parseMultipoleType(String input, String[] tokens, BufferedReader br) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
        } else {
            try {
                int numTypes = tokens.length - 2;
                int[] atomTypes = new int[numTypes];
                for (int i = 0; i < numTypes; i++) {
                    atomTypes[i] = parseInt(tokens[i + 1]);
                }
                MultipoleType.MultipoleFrameDefinition frameDefinition
                        = MultipoleType.MultipoleFrameDefinition.ZTHENX;
                if (atomTypes.length == 3 && (atomTypes[1] < 0 || atomTypes[2] < 0)) {
                    frameDefinition = MultipoleType.MultipoleFrameDefinition.BISECTOR;
                } else if (atomTypes.length == 4 && atomTypes[2] < 0 && atomTypes[3] < 0) {
                    if (atomTypes[1] < 0) {
                        frameDefinition = MultipoleType.MultipoleFrameDefinition.TRISECTOR;
                    } else {
                        frameDefinition = MultipoleType.MultipoleFrameDefinition.ZTHENBISECTOR;
                    }
                }
                for (int i = 0; i < numTypes; i++) {
                    atomTypes[i] = abs(atomTypes[i]);
                }
                double c = parseDouble(tokens[1 + numTypes]);
                input = br.readLine().split("#")[0];
                tokens = input.trim().split(" +");
                if (tokens.length != 3) {
                    logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
                    return null;
                }
                double[] dipole = new double[3];
                dipole[0] = parseDouble(tokens[0]);
                dipole[1] = parseDouble(tokens[1]);
                dipole[2] = parseDouble(tokens[2]);
                input = br.readLine().split("#")[0];
                tokens = input.trim().split(" +");
                if (tokens.length != 1) {
                    logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
                    return null;
                }
                double[][] quadrupole = new double[3][3];
                quadrupole[0][0] = parseDouble(tokens[0]);
                input = br.readLine().split("#")[0];
                tokens = input.trim().split(" +");
                if (tokens.length != 2) {
                    logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
                    return null;
                }
                quadrupole[1][0] = parseDouble(tokens[0]);
                quadrupole[1][1] = parseDouble(tokens[1]);
                input = br.readLine().split("#")[0];
                tokens = input.trim().split(" +");
                if (tokens.length != 3) {
                    logger.log(Level.WARNING, "Invalid MULTIPOLE type:\n{0}", input);
                    return null;
                }
                quadrupole[2][0] = parseDouble(tokens[0]);
                quadrupole[2][1] = parseDouble(tokens[1]);
                quadrupole[2][2] = parseDouble(tokens[2]);
                // Fill in symmetric components.
                quadrupole[0][1] = quadrupole[1][0];
                quadrupole[0][2] = quadrupole[2][0];
                quadrupole[1][2] = quadrupole[2][1];
                MultipoleType multipoleType = new MultipoleType(c, dipole,
                        quadrupole, atomTypes, frameDefinition, true);
                forceField.addForceFieldType(multipoleType);
                return multipoleType;
            } catch (NumberFormatException | IOException e) {
                String message = "Exception parsing MULTIPOLE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * Parse a single line multipole.
     *
     * @param input  Input String.
     * @param tokens Input tokens.
     * @since 1.0
     */
    private MultipoleType parseMultipoleType(String input, String[] tokens) {
        if (tokens.length < 14) {
            logger.log(Level.WARNING, "Invalid MULTIPOLE type:{0}", Arrays.toString(tokens));
        } else {
            try {
                int numTypes = tokens.length - 11;
                int[] atomTypes = new int[numTypes];
                for (int i = 0; i < numTypes; i++) {
                    atomTypes[i] = parseInt(tokens[i + 1]);
                }
                MultipoleType.MultipoleFrameDefinition frameDefinition = MultipoleType.MultipoleFrameDefinition.ZTHENX;
                if (atomTypes.length == 3 && (atomTypes[1] < 0 || atomTypes[2] < 0)) {
                    frameDefinition = MultipoleType.MultipoleFrameDefinition.BISECTOR;
                } else if (atomTypes.length == 4 && atomTypes[2] < 0 && atomTypes[3] < 0) {
                    if (atomTypes[1] < 0) {
                        frameDefinition = MultipoleType.MultipoleFrameDefinition.TRISECTOR;
                    } else {
                        frameDefinition = MultipoleType.MultipoleFrameDefinition.ZTHENBISECTOR;
                    }
                }
                for (int i = 0; i < numTypes; i++) {
                    atomTypes[i] = abs(atomTypes[i]);
                }
                double[] dipole = new double[3];
                double[][] quadrupole = new double[3][3];
                double c = parseDouble(tokens[1 + numTypes]);
                dipole[0] = parseDouble(tokens[2 + numTypes]);
                dipole[1] = parseDouble(tokens[3 + numTypes]);
                dipole[2] = parseDouble(tokens[4 + numTypes]);
                quadrupole[0][0] = parseDouble(tokens[5 + numTypes]);
                quadrupole[1][0] = parseDouble(tokens[6 + numTypes]);
                quadrupole[1][1] = parseDouble(tokens[7 + numTypes]);
                quadrupole[2][0] = parseDouble(tokens[8 + numTypes]);
                quadrupole[2][1] = parseDouble(tokens[9 + numTypes]);
                quadrupole[2][2] = parseDouble(tokens[10 + numTypes]);
                // Fill in symmetric components.
                quadrupole[0][1] = quadrupole[1][0];
                quadrupole[0][2] = quadrupole[2][0];
                quadrupole[1][2] = quadrupole[2][1];
                MultipoleType multipoleType = new MultipoleType(c, dipole,
                        quadrupole, atomTypes, frameDefinition, true);
                forceField.addForceFieldType(multipoleType);
                return multipoleType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing MULTIPOLE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private OutOfPlaneBendType parseOutOfPlaneBendType(String input, String[] tokens) {
        if (tokens.length < 6) {
            logger.log(Level.WARNING, "Invalid OPBEND type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);
                double forceConstant = parseDouble(tokens[5]);
                OutOfPlaneBendType outOfPlaneBendType = new OutOfPlaneBendType(atomClasses, forceConstant);
                forceField.addForceFieldType(outOfPlaneBendType);
                return outOfPlaneBendType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing OPBEND type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private PiTorsionType parsePiTorsionType(String input, String[] tokens) {
        if (tokens.length < 4) {
            logger.log(Level.WARNING, "Invalid PITORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[2];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                double forceConstant = parseDouble(tokens[3]);
                PiTorsionType piTorsionType = new PiTorsionType(atomClasses,
                        forceConstant);
                forceField.addForceFieldType(piTorsionType);
                return piTorsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing PITORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private PolarizeType parsePolarizeType(String input, String[] tokens) {
        if (tokens.length < 4) {
            logger.log(Level.WARNING, "Invalid POLARIZE type:\n{0}", input);
        } else {
            try {
                int atomType = parseInt(tokens[1]);
                double polarizability = parseDouble(tokens[2]);
                double thole = parseDouble(tokens[3]);
                int entries = tokens.length - 4;
                int[] polarizationGroup = null;
                if (entries > 0) {
                    polarizationGroup = new int[entries];
                    for (int i = 4; i < tokens.length; i++) {
                        polarizationGroup[i - 4] = parseInt(tokens[i]);
                    }
                }
                PolarizeType polarizeType = new PolarizeType(atomType,
                        polarizability, thole, polarizationGroup);
                forceField.addForceFieldType(polarizeType);
                return polarizeType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing POLARIZE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private StretchBendType parseStretchBendType(String input, String[] tokens) {
        if (tokens.length < 6) {
            logger.log(Level.WARNING, "Invalid STRBND type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[3];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                double[] forceConstants = new double[2];
                forceConstants[0] = parseDouble(tokens[4]);
                forceConstants[1] = parseDouble(tokens[5]);
                StretchBendType stretchBendType = new StretchBendType(atomClasses,
                        forceConstants);
                forceField.addForceFieldType(stretchBendType);
                return stretchBendType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing STRBND type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private ImproperTorsionType parseImpTorsType(String input, String[] tokens) {
        if (tokens.length < 8) {
            logger.log(Level.WARNING, "Invalid IMPTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);
                double k = parseDouble(tokens[5]);
                double phase = parseDouble(tokens[6]);
                int period = parseInt(tokens[7]);
                if (improperTorsionScale > 0.0) {
                    k = k * improperTorsionScale;
                }
                ImproperTorsionType improperTorsionType = new ImproperTorsionType(atomClasses,
                        k, phase, period);
                forceField.addForceFieldType(improperTorsionType);
                return improperTorsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing IMPTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private TorsionType parseTorsion(String input, String[] tokens) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid TORSION type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);
                int terms = (tokens.length - 5) / 3;
                double[] amplitude = new double[terms];
                double[] phase = new double[terms];
                int[] periodicity = new int[terms];
                int index = 5;
                for (int i = 0; i < terms; i++) {
                    amplitude[i] = parseDouble(tokens[index++]);
                    phase[i] = parseDouble(tokens[index++]);
                    periodicity[i] = parseInt(tokens[index++]);
                    if (torsionScale > 0.0) {
                        amplitude[i] = amplitude[i] * torsionScale;
                    }
                }
                TorsionType torsionType = new TorsionType(atomClasses, amplitude,
                        phase, periodicity);
                forceField.addForceFieldType(torsionType);
                return torsionType;
            } catch (NumberFormatException e) {
                String message = "NumberFormatException parsing TORSION type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            } catch (Exception e) {
                String message = "Exception parsing TORSION type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private TorsionType parseImproperTorsion(String input, String[] tokens) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid IMPROPER type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);
                double[] amplitude = new double[1];
                double[] phase = new double[1];
                int[] periodicity = new int[1];
                int index = 5;
                amplitude[0] = parseDouble(tokens[index++]);
                phase[0] = parseDouble(tokens[index++]);
                periodicity[0] = 1;
                TorsionType torsionType = new TorsionType(atomClasses, amplitude,
                        phase, periodicity, TorsionType.TorsionMode.IMPROPER);
                forceField.addForceFieldType(torsionType);
                return torsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing IMPROPER type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private StretchTorsionType parseStretchTorsion(String input, String[] tokens) {
        if (tokens.length < 13) {
            logger.log(Level.WARNING, "Invalid STRTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);

                double[] constants = new double[9];
                constants[0] = parseDouble(tokens[5]);
                constants[1] = parseDouble(tokens[6]);
                constants[2] = parseDouble(tokens[7]);
                constants[3] = parseDouble(tokens[8]);
                constants[4] = parseDouble(tokens[9]);
                constants[5] = parseDouble(tokens[10]);
                constants[6] = parseDouble(tokens[11]);
                constants[7] = parseDouble(tokens[12]);
                constants[8] = parseDouble(tokens[13]);

                StretchTorsionType stretchTorsionType = new StretchTorsionType(atomClasses, constants);
                forceField.addForceFieldType(stretchTorsionType);
                return stretchTorsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing STRTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private AngleTorsionType parseAngleTorsion(String input, String[] tokens) {
        if (tokens.length < 10) {
            logger.log(Level.WARNING, "Invalid ANGTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);

                double[] constants = new double[6];
                constants[0] = parseDouble(tokens[5]);
                constants[1] = parseDouble(tokens[6]);
                constants[2] = parseDouble(tokens[7]);
                constants[3] = parseDouble(tokens[8]);
                constants[4] = parseDouble(tokens[9]);
                constants[5] = parseDouble(tokens[10]);

                AngleTorsionType angleTorsionType = new AngleTorsionType(atomClasses, constants);
                forceField.addForceFieldType(angleTorsionType);
                return angleTorsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing ANGTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private TorsionTorsionType parseTorsionTorsionType(String input, String[] tokens, BufferedReader br) {
        if (tokens.length < 8) {
            logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[5];
                for (int i = 0; i < 5; i++) {
                    atomClasses[i] = parseInt(tokens[i + 1]);
                }
                int[] gridPoints = new int[2];
                gridPoints[0] = parseInt(tokens[6]);
                gridPoints[1] = parseInt(tokens[7]);
                int points = gridPoints[0] * gridPoints[1];
                double[] torsion1 = new double[points];
                double[] torsion2 = new double[points];
                double[] energy = new double[points];
                for (int i = 0; i < points; i++) {
                    input = br.readLine();
                    tokens = input.trim().split(" +");
                    if (tokens.length != 3) {
                        logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
                        return null;
                    }
                    torsion1[i] = parseDouble(tokens[0]);
                    torsion2[i] = parseDouble(tokens[1]);
                    energy[i] = parseDouble(tokens[2]);
                }
                TorsionTorsionType torsionTorsionType = new TorsionTorsionType(
                        atomClasses, gridPoints, torsion1, torsion2, energy);
                forceField.addForceFieldType(torsionTorsionType);
                return torsionTorsionType;
            } catch (NumberFormatException | IOException e) {
                String message = "Exception parsing TORTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private TorsionTorsionType parseTorsionTorsionType(String input, String[] tokens) {
        if (tokens.length < 8) {
            logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[5];
                for (int i = 0; i < 5; i++) {
                    atomClasses[i] = parseInt(tokens[i + 1]);
                }
                int[] gridPoints = new int[2];
                gridPoints[0] = parseInt(tokens[6]);
                gridPoints[1] = parseInt(tokens[7]);
                int points = gridPoints[0] * gridPoints[1];
                int numTokens = points * 3 + 8;
                if (tokens.length < numTokens) {
                    logger.log(Level.WARNING, "Invalid TORTORS type:\n{0}", input);
                    return null;
                }
                double[] torsion1 = new double[points];
                double[] torsion2 = new double[points];
                double[] energy = new double[points];
                int index = 8;
                for (int i = 0; i < points; i++) {
                    torsion1[i] = parseDouble(tokens[index++]);
                    torsion2[i] = parseDouble(tokens[index++]);
                    energy[i] = parseDouble(tokens[index++]);
                }
                TorsionTorsionType torsionTorsionType = new TorsionTorsionType(
                        atomClasses, gridPoints, torsion1, torsion2, energy);
                forceField.addForceFieldType(torsionTorsionType);
                return torsionTorsionType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing TORTORS type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private UreyBradleyType parseUreyBradleyType(String input, String[] tokens) {
        if (tokens.length < 5) {
            logger.log(Level.WARNING, "Invalid UREYBRAD type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[3];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                double forceConstant = parseDouble(tokens[4]);
                double distance = parseDouble(tokens[5]);
                UreyBradleyType ureyBradleyType = new UreyBradleyType(atomClasses,
                        forceConstant, distance);
                forceField.addForceFieldType(ureyBradleyType);
                return ureyBradleyType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing UREYBRAD type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private VDWType parseVDWType(String input, String[] tokens) {
        if (tokens.length < 4) {
            logger.log(Level.WARNING, "Invalid VDW type:\n{0}", input);
        } else {
            try {
                int atomType = parseInt(tokens[1]);
                double radius = parseDouble(tokens[2]);
                double wellDepth = parseDouble(tokens[3]);
                double reductionFactor = -1.0;
                if (tokens.length == 5) {
                    reductionFactor = parseDouble(tokens[4]);
                }
                if (convertRadiusToDiameter) {
                    radius = radius * 2.0;
                }
                if (convertSigmaToRMin) {
                    double twoSix = 1.122462048309372981;
                    radius = radius * twoSix;
                }
                VDWType vdwType = new VDWType(atomType, radius, wellDepth, reductionFactor);
                forceField.addForceFieldType(vdwType);
                return vdwType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing VDW type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private VDWType parseVDW14Type(String input, String[] tokens) {
        if (tokens.length < 4) {
            logger.log(Level.WARNING, "Invalid VDW type:\n{0}", input);
        } else {
            try {
                int atomType = parseInt(tokens[1]);
                double radius = parseDouble(tokens[2]);
                double wellDepth = parseDouble(tokens[3]);
                double reductionFactor = -1.0;
                if (tokens.length == 5) {
                    reductionFactor = parseDouble(tokens[4]);
                }
                if (convertRadiusToDiameter) {
                    radius = radius * 2.0;
                }
                if (convertSigmaToRMin) {
                    double twoSix = 1.122462048309372981;
                    radius = radius * twoSix;
                }
                VDWType vdwType = new VDWType(atomType, radius, wellDepth, reductionFactor, VDWType.VDWMode.VDW14);
                forceField.addForceFieldType(vdwType);
                return vdwType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing VDW14 type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    private SoluteType parseSolute(String input, String[] tokens) {
        if (tokens.length < 4) {
            logger.log(Level.WARNING, "Invalid SOLUTE type:\n{0}", input);
        } else {
            try {
                int atomType = parseInt(tokens[1].trim());
                String description = tokens[2].trim();
                double diameter = parseDouble(tokens[3].trim());
                SoluteType soluteType = new SoluteType(atomType, description, diameter);
                forceField.addForceFieldType(soluteType);
                return soluteType;
            } catch (NumberFormatException e) {
                String message = "Exception parsing SOLUTE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
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
