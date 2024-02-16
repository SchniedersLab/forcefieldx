//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
//package ffx.potential.groovy.test

import ffx.potential.MolecularAssembly
import ffx.potential.Utilities
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.bonded.Molecule
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.AngleType
import ffx.potential.parameters.AtomType
import ffx.potential.parameters.BondType
import ffx.potential.parameters.ChargeType
import ffx.potential.parameters.ForceField
import ffx.potential.parameters.TorsionType
import ffx.potential.parameters.VDWType
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.XYZFilter
import ffx.utilities.Keyword
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import java.util.logging.Level;
import org.w3c.dom.Document
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.w3c.dom.NodeList
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import javax.xml.parsers.DocumentBuilder
import javax.xml.parsers.DocumentBuilderFactory

import static java.lang.String.format
import static ffx.potential.bonded.Atom.ElementSymbol;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.PI;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;


/**
 * The XMLtoPatch script converts an Open Force Field XML parameter file to Patch format.
 *
 * @author Aaron J. Nessler
 * <br>
 * Usage:
 * <br>
 * ffxc test.XMLtoPatch &lt;filename&gt;
 */
@Command(description = " Convert XML file to Patch.", name = "ffxc test.XMLtoPatch")
class XMLtoPatch extends PotentialScript {

    /**
     * The final argument(s) should be two filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'XML file to be converted to Patch.')
    List<String> filenames = null;

    final double ANGperNM = 10.0;
    final double KJperKCal = 4.184;

    /**
     * Execute the script.
     */
    @Override
    XMLtoPatch run() {

        // Init the context and bind variables.
        if (!init()) {
            return this;
        }

        for (String filename : filenames) {
            if (FilenameUtils.getExtension(filename).containsIgnoreCase("XML")) {
                try {
                    CompositeConfiguration properties = Keyword.loadProperties(null);
                    properties.addProperty("renumberPatch", "FALSE");
                    properties.addProperty("FORCEFIELD", "AMBER-2023");
                    properties.addProperty("VDWTYPE", "LENNARD-JONES");
                    properties.addProperty("RADIUSRULE", "ARITHMETIC");
                    properties.addProperty("RADIUSTYPE", "SIGMA");
                    properties.addProperty("RADIUSSIZE", "DIAMETER");
                    properties.addProperty("EPSILONRULE", "GEOMETRIC");
                    properties.addProperty("VDW-14-SCALE", "0.500000");
                    properties.addProperty("CHG-14-SCALE", "1.2");
                    properties.addProperty("ELECTRIC", "332.0716");
                    properties.addProperty("DIELECTRIC", "1.0");
                    ForceField forceField = new ForceField(properties);

                    File inputFile = new File(filenames[0]);
                    String fileName = FilenameUtils.removeExtension(inputFile.getName());
                    DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
                    DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
                    Document doc = dBuilder.parse(inputFile);
                    // Node "System"
                    logger.info(format(" Root Element: %s", doc.getDocumentElement().getNodeName()));

                    StringBuilder sb = new StringBuilder();
                    sb.append("renumberPatch         False\n");
                    sb.append("\n\n     ##############################\n" +
                            "     ##                          ##\n" +
                            "     ##  Force Field Definition  ##\n" +
                            "     ##                          ##\n" +
                            "     ##############################\n\n\n" +
                            "FORCEFIELD                 AMBER-2023\n" +
                            "VDWTYPE                    LENNARD-JONES\n" +
                            "RADIUSRULE               ARITHMETIC\n" +
                            "RADIUSTYPE               SIGMA\n" +
                            "RADIUSSIZE                DIAMETER\n" +
                            "EPSILONRULE               GEOMETRIC\n" +
                            "VDW-14-SCALE               0.500000\n" +
                            "CHG-14-SCALE               1.2\n" +
                            "ELECTRIC               332.0716\n" +
                            "DIELECTRIC               1.0\n\n\n");
                    sb.append("\n\n     ###############################\n" +
                            "     ##                           ##\n" +
                            "     ##   Atom Type Definitions   ##\n" +
                            "     ##                           ##\n" +
                            "     ###############################\n\n");
                    StringBuilder sbVdw = new StringBuilder();
                    sbVdw.append("\n\n     ##############################\n" +
                            "     ##                          ##\n" +
                            "     ## Van der Waals Parameters ##\n" +
                            "     ##                          ##\n" +
                            "     ##############################\n\n");
                    StringBuilder sbBond = new StringBuilder();
                    sbBond.append("\n\n     ################################\n" +
                            "     ##                            ##\n" +
                            "     ## Bond Stretching Parameters ##\n" +
                            "     ##                            ##\n" +
                            "     ################################\n\n");
                    StringBuilder sbAngle = new StringBuilder();
                    sbAngle.append("\n\n     ################################\n" +
                            "     ##                            ##\n" +
                            "     ##  Angle Bending Parameters  ##\n" +
                            "     ##                            ##\n" +
                            "     ################################\n\n");
                    StringBuilder sbTor = new StringBuilder();
                    sbTor.append("\n\n     ################################\n" +
                            "     ##                            ##\n" +
                            "     ##    Torsional Parameters    ##\n" +
                            "     ##                            ##\n" +
                            "     ################################\n\n");
                    StringBuilder sbAPC = new StringBuilder();
                    sbAPC.append("\n\n     ######################################\n" +
                            "     ##                                  ##\n" +
                            "     ## Atomic Partial Charge Parameters ##\n" +
                            "     ##                                  ##\n" +
                            "     ######################################\n\n");
                    NodeList nodeList = doc.getChildNodes();
                    // Assumed one system label for now...
                    Node node = nodeList.item(0);
                    NodeList childNodes = node.getChildNodes();
                    logger.info(" Child Size:" + childNodes.length);
                    int particlesLength = -1;
                    String[] atomName;
                    double[] atomMass;
                    int[] atomNumber;
                    int[] valence;
                    for (Node child : childNodes) {
                        if (child.getNodeType() == Node.ELEMENT_NODE) {
                            Element e = (Element) child;
                            if(logger.isLoggable(Level.FINER)){
                                logger.finer(format(" Node Name: %s", child.getNodeName()));
                            }
                            switch (e.getNodeName()) {
                                case "PeriodicBoxVectors":
                                    NodeList boxList = e.getChildNodes();
                                    int boxLength = boxList.length;
                                    for (int i = 0; i < boxLength; i++) {
                                        Node box = boxList.item(i);
                                        if (box.getNodeName() == "A" || box.getNodeName() == "B" || box.getNodeName() == "C" && box.getNodeType() == Node.ELEMENT_NODE) {
                                            Element element = (Element) box;
                                            String x = element.getAttribute("x"); // in nm
                                            String y = element.getAttribute("y"); // in nm
                                            String z = element.getAttribute("z"); // in nm
                                            if(logger.isLoggable(Level.FINER)) {
                                                logger.finer(box.getNodeName());
                                                logger.finer(" x: " + x + " nm");
                                                logger.finer(" y: " + y + " nm");
                                                logger.finer(" z: " + z + " nm");
                                            }
                                        }
                                    }
                                    break;
                                case "Particles":
                                    NodeList particleList = e.getElementsByTagName("Particle");
                                    if (particlesLength >= 0 && particlesLength != particleList.length) {
                                        logger.warning(" Updated number of particles...");
                                    } else if (particlesLength < 0) {
                                        particlesLength = particleList.length;
                                    }
                                    logger.info(" Particle Num Mass: " + particlesLength);
                                    atomName = new String[particlesLength];
                                    atomMass = new double[particlesLength];
                                    atomNumber = new int[particlesLength];
                                    valence = new int[particlesLength];
                                    double[] atomicMass = AtomType.atomicMass
                                    int numAtomTable = atomicMass.length;
                                    int[] atomCount = new int[numAtomTable];
                                    for (int i = 0; i < particlesLength; i++) {
                                        Node particle = particleList.item(i);
                                        if (particle.getNodeName() == "Particle" && particle.getNodeType() == Node.ELEMENT_NODE) {
                                            Element element = (Element) particle;
                                            String mass = element.getAttribute("mass"); // Atomic mass units
                                            atomMass[i] = parseDouble(mass);
                                            if(logger.isLoggable(Level.FINER)) {
                                                logger.finer(" Mass: " + atomMass[i] + " AMU");
                                            }
                                            for (int j = 0; j < numAtomTable; j++) {
                                                if (abs(atomicMass[j] - atomMass[i]) < 0.25) {
                                                    atomNumber[i] = j + 1;
                                                    if(logger.isLoggable(Level.FINER)) {
                                                        logger.finer(" Atom Number: " + atomNumber[i]);
                                                    }
                                                    atomName[i] = ElementSymbol.values()[j].toString() + ++atomCount[j];
                                                    if(logger.isLoggable(Level.FINER)) {
                                                        logger.finer(" Atom Name: " + atomName[i]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    break;
                                case "Constraints":
                                    // Not encountered yet, but I know it is an option.
                                    logger.warning(" XML field \"Constraints\" has not been implemented.");
                                    break;
                                case "Forces":
                                    NodeList forceList = e.getChildNodes();
                                    int forceLength = forceList.length;
                                    for (int i = 0; i < forceLength; i++) {
                                        Node force = forceList.item(i);
                                        if (force.getNodeType() == Node.ELEMENT_NODE) {
                                            Element element = (Element) force;
                                            String elementName = element.getAttribute("name");
                                            if(logger.isLoggable(Level.FINER)) {
                                                logger.finer(format(" Node Name: %s", force.getNodeName()))
                                                logger.finer(format(" Attribute Name: %s", elementName));
                                            }

                                            switch (elementName) {
                                                case "NonbondedForce":
                                                    String alpha = element.getAttribute("alpha"); // PME separation parameter
                                                    String cutoff = element.getAttribute("cutoff"); // in nm
                                                    String dispersionCorrection = element.getAttribute("dispersionCorrection"); // Constant energy added beyond cutoff
                                                    String ewaldTolerance = element.getAttribute("ewaldTolerance"); // fractional tolerance
                                                    String exceptionsUsePeriodic = element.getAttribute("exceptionsUsePeriodic");
                                                    String forceGroup = element.getAttribute("forceGroup"); // 0 <= fg <= 31
                                                    String includeDirectSpace = element.getAttribute("includeDirectSpace");
                                                    String ljAlpha = element.getAttribute("ljAlpha"); // Leonard-Jones alpha
                                                    String ljnx = element.getAttribute("ljnx"); // LJ PME grid points
                                                    String ljny = element.getAttribute("ljny"); // LJ PME grid points
                                                    String ljnz = element.getAttribute("ljnz"); // LJ PME grid points
                                                    String method = element.getAttribute("method");
                                                    String name = element.getAttribute("name"); // arbitrary ID
                                                    String nx = element.getAttribute("nx"); // PME grid points
                                                    String ny = element.getAttribute("ny"); // PME grid points
                                                    String nz = element.getAttribute("nz"); // PME grid points
                                                    String recipForceGroup = element.getAttribute("recipForceGroup");// -1 --> same as direct
                                                    String rfDielectric = element.getAttribute("rfDielectric"); // reaction field approximation solvent dielectric
                                                    String switchingDistance = element.getAttribute("switchingDistance"); // For LJ in nm
                                                    String type = element.getAttribute("type");
                                                    String useSwitchingFunction = element.getAttribute("useSwitchingFunction");
                                                    String version = element.getAttribute("version");
                                                    if(logger.isLoggable(Level.FINE)) {
                                                        logger.fine(format(" PARAMETERS: alpha %s, nx %s, ny %s, nz %s, ljA %s, ljx %s, ljy %s, ljz %s, method %s, cutoff %s" +
                                                                ",\n dispersion %s, ewaldTol: %s, exceptionPME: %s, forceGroup %s, recipFG %s," +
                                                                "\n includeDirect %s, name %s, rfDielectric %s, switchDist %s, type %s, useSwitch %s, version %s",
                                                                alpha, nx, ny, nz, ljAlpha, ljnx, ljny, ljnz, method, cutoff,
                                                                dispersionCorrection, ewaldTolerance, exceptionsUsePeriodic, forceGroup,
                                                                recipForceGroup, includeDirectSpace, name, rfDielectric, switchingDistance, type, useSwitchingFunction, version));
                                                    }
                                                    NodeList pList = element.getElementsByTagName("Particle");
                                                    if (particlesLength >= 0 && particlesLength != pList.length) {
                                                        logger.warning(" Updated number of particles ...");
                                                        particlesLength = pList.length;
                                                    } else if (particlesLength < 0) {
                                                        particlesLength = pList.length;
                                                    }
                                                    logger.info(" Particle Length: " + particlesLength);
                                                    for (int j = 0; j < particlesLength; j++) {
                                                        Element particle = (Element) pList.item(j);
                                                        String eps = particle.getAttribute("eps"); // kJ/mol
                                                        String q = particle.getAttribute("q"); // in proton units
                                                        String sig = particle.getAttribute("sig"); // in nm
                                                        sbVdw.append(format("vdw     %3d     %12.8f     %7.4f\n", j + 1, parseDouble(sig) * ANGperNM, parseDouble(eps) / KJperKCal));
                                                        double sigma = parseDouble(sig) * ANGperNM;
                                                        double epsilon = parseDouble(eps) / KJperKCal;
                                                        forceField.addForceFieldType(new VDWType(j + 1, sigma, epsilon, -1.0));
                                                        sbAPC.append(format("charge          %3d          %7.4f\n", j + 1, parseDouble(q)));
                                                        forceField.addForceFieldType(new ChargeType(j + 1, parseDouble(q)));
                                                        if(logger.isLoggable(Level.FINER)) {
                                                            logger.finer(format(" Particle: eps %16s, q %16s, sig %16s", eps, q, sig));
                                                        }
                                                    }
                                                    NodeList exList = element.getElementsByTagName("Exception");
                                                    int numExceptions = exList.length;
                                                    logger.finer(" Exception Length: " + numExceptions);
                                                    for (int j = 0; j < numExceptions; j++) {
                                                        Element exception = (Element) exList.item(j);
                                                        String eps = exception.getAttribute("eps"); // kJ/mol
                                                        String p1 = exception.getAttribute("p1");
                                                        String p2 = exception.getAttribute("p2");
                                                        String q = exception.getAttribute("q"); // in proton units
                                                        String sig = exception.getAttribute("sig"); // in nm
                                                        if(logger.isLoggable(Level.FINER)) {
                                                            logger.finer(format(" Exception: eps %16s, p1 %2s, p2 %2s, q %16s, sig %16s", eps, p1, p2, q, sig));
                                                        }
                                                    }
                                                    break;
                                                case "PeriodicTorsionForce":
                                                    Map<String, Map<Integer, String>> torsionMap = new HashMap();
                                                    NodeList tList = element.getElementsByTagName("Torsion");
                                                    int numTorsions = tList.length;
                                                    logger.info(" Torsions Length: " + numTorsions);
                                                    for (int j = 0; j < numTorsions; j++) {
                                                        Element torsion = (Element) tList.item(j);
                                                        String k = torsion.getAttribute("k");
                                                        String p1 = torsion.getAttribute("p1");
                                                        String p2 = torsion.getAttribute("p2");
                                                        String p3 = torsion.getAttribute("p3");
                                                        String p4 = torsion.getAttribute("p4");
                                                        String periodicity = torsion.getAttribute("periodicity");
                                                        int iPeriodicity = parseInt(periodicity);
                                                        String phase = torsion.getAttribute("phase"); // Radians
                                                        String key = format("torsion          %2d     %2d     %2d     %2d          ",
                                                                parseInt(p1) + 1, parseInt(p2) + 1, parseInt(p3) + 1, parseInt(p4) + 1);
                                                        Map<Integer, String> previousMap = new HashMap<>();
                                                        if (torsionMap.containsKey(key)) {
                                                            previousMap = (Map<Integer, String>) torsionMap.get(key);
                                                        }
                                                        previousMap.put(iPeriodicity, format("%5.3f     %5.1f     ", parseDouble(k) / KJperKCal, parseDouble(phase) * (180 / PI)));
                                                        torsionMap.put(key, previousMap);
                                                        int[] atomClass = new int[]{parseInt(p1) + 1, parseInt(p2) + 1, parseInt(p3) + 1, parseInt(p4) + 1};
                                                        int[] period = new int[]{iPeriodicity};
                                                        double[] amplitude = new double[]{parseDouble(k) / KJperKCal};
                                                        double[] phaseShift = new double[]{parseDouble(phase) * (180 / PI)};
                                                        forceField.addForceFieldType(new TorsionType(atomClass, amplitude, phaseShift, period, TorsionType.TorsionMode.NORMAL))
                                                        if(logger.isLoggable(Level.FINE)) {
                                                            logger.info(format(" Torsion: k %16s, p1 %2s, p2 %2s, p3 %2s, p4 %2s, period %2s, phase %16s", k, p1, p2, p3, p4, periodicity, phase));
                                                        }
                                                    }
//                                            sbTor.append(format("torsion          %2d     %2d     %2d     %2d          %5.3f     %5.1f     %2d\n",
//                                                    parseInt(p1) + 1, parseInt(p2) + 1, parseInt(p3) + 1, parseInt(p4) + 1,
//                                                    parseDouble(k)/kJperKCal, parseDouble(phase) * (180/PI), iPeriodicity));
                                                    for (Map.Entry<String, Map<Integer, String>> entry : torsionMap.entrySet()) {
                                                        sbTor.append(format(entry.getKey()));
                                                        Map<Integer, String> map = entry.getValue();
                                                        for (Map.Entry<Integer, String> innerEntry : map.entrySet()) {
                                                            sbTor.append(format(innerEntry.getValue() + innerEntry.getKey() + "     "));
                                                        }
                                                        sbTor.append("\n");
                                                    }
                                                    break;
                                                case "HarmonicAngleForce":
                                                    NodeList aList = element.getElementsByTagName("Angle");
                                                    int numAngles = aList.length;
                                                    logger.info(" Angle Length: " + numAngles);
                                                    for (int j = 0; j < numAngles; j++) {
                                                        Element angle = (Element) aList.item(j);
                                                        String a = angle.getAttribute("a"); // Radians
                                                        String k = angle.getAttribute("k"); // kJ/mol/rad^2
                                                        String p1 = angle.getAttribute("p1");
                                                        String p2 = angle.getAttribute("p2");
                                                        String p3 = angle.getAttribute("p3");
                                                        // Angle force is ~2x what it should be via comparison... divide by 2... Conversion radius <--> diameter?
                                                        sbAngle.append(format("angle          %2d     %2d     %2d     %5.2f     %6.2f\n",
                                                                parseInt(p1) + 1, parseInt(p2) + 1, parseInt(p3) + 1, parseDouble(k) / KJperKCal / 2, parseDouble(a) * (180.0 / PI)));
                                                        int[] atomClass = new int[]{parseInt(p1) + 1, parseInt(p2) + 1, parseInt(p3) + 1}
                                                        double forceConstant = parseDouble(k) / KJperKCal / 2;
                                                        double[] angles = new double[]{parseDouble(a) * (180.0 / PI)};
                                                        forceField.addForceFieldType(new AngleType(atomClass, forceConstant, angles));
                                                        if(logger.isLoggable(Level.FINE)) {
                                                            logger.info(format(" Angle: a %16s, p1 %2s, p2 %2s, p3 %2s", a, k, p1, p2, p3));
                                                        }
                                                    }
                                                    break;
                                                case "HarmonicBondForce":
                                                    NodeList bList = element.getElementsByTagName("Bond");
                                                    int numBonds = bList.length;
                                                    logger.info(" Bond Length: " + numBonds);
                                                    if (valence == null) {
                                                        // This better not happen...
                                                        valence = new int[500];
                                                        logger.warning(" Connectivity size was undefined.");
                                                    }
                                                    logger.info("Connect Length: " + valence.length);
                                                    for (int j = 0; j < numBonds; j++) {
                                                        Element bond = (Element) bList.item(j);
                                                        String d = bond.getAttribute("d"); //length in nm
                                                        String k = bond.getAttribute("k"); // kJ/mol/nm^2
                                                        String p1 = bond.getAttribute("p1");
                                                        String p2 = bond.getAttribute("p2");
                                                        int p1int = parseInt(p1);
                                                        int p2int = parseInt(p2);
                                                        valence[p1int]++;
                                                        valence[p2int]++;
                                                        // TODO investigate comment below:
                                                        // Bond force is ~2x what it should be via comparison... divide by 2... Conversion radius <--> diameter?
                                                        sbBond.append(format("bond     %2d     %2d          %6.2f     %6.4f\n",
                                                                p1int + 1, p2int + 1, parseDouble(k) / KJperKCal / (ANGperNM * ANGperNM) / 2, parseDouble(d) * ANGperNM));
                                                        if(logger.isLoggable(Level.FINE)) {
                                                            logger.fine(format(" Bond: d %16s, k %16s, p1 %2s, p2 %2s", d, k, p1, p2));
                                                        }
                                                        int[] atomClass = new int[]{p1int + 1, p2int + 1};
                                                        double forceConstant = parseDouble(k) / KJperKCal / (ANGperNM * ANGperNM) / 2;
                                                        double distance = parseDouble(d) * ANGperNM;
                                                        forceField.addForceFieldType(new BondType(atomClass, forceConstant, distance));
                                                    }
                                                    break;
                                                default:
                                                    logger.warning(" Unknown force element: " + elementName);
                                                    break;
                                            }
                                        }
                                    }
                                    break;
                                default:
                                    logger.info("404: Node not found.");
                            }
                        }
                    }
                    for (int i = 0; i < particlesLength; i++) {
                        if (atomName == null || atomNumber == null || atomMass == null || valence == null) {
                            logger.warning("Atom record values were not detected from XML.");
                            break;
                        }
                        sb.append(format("atom     %3d     %3d     %5s     \"%5s\"          %3d     %6.3f     %1d\n",
                                i + 1, i + 1, atomName[i], atomName[i], atomNumber[i], atomMass[i], valence[i]))
                        forceField.addForceFieldType(new AtomType(i + 1, i + 1, atomName[i], "\"" + atomName[i] + "\"", atomNumber[i], atomMass[i], valence[i]))
                    }
                    if(logger.isLoggable(Level.FINER)) {
                        logger.finer("FORCE FIELD OBJECT: \n\n\n" + forceField.toStringBuffer().toString() + "END FORCE FIELD OBJECT.\n\n");
                    }
//                    //Create patch file.
//                    File ffOutFile = PDBFilter.version(new File(fileName + "_tmp.patch"));
//                    if (ffOutFile.createNewFile()) {
//                        BufferedWriter bw = new BufferedWriter(new FileWriter(ffOutFile, false));
//                        bw.write("renumberPatch                 " + properties.getString("renumberPatch") + "\n");
//                        bw.write("FORCEFIELD                 " + properties.getString("FORCEFIELD") + "\n");
//                        bw.write("VDWTYPE                 " + properties.getString("VDWTYPE") + "\n");
//                        bw.write("RADIUSRULE                 " + properties.getString("RADIUSRULE") + "\n");
//                        bw.write("RADIUSTYPE                 " + properties.getString("RADIUSTYPE") + "\n");
//                        bw.write("RADIUSSIZE                 " + properties.getString("RADIUSSIZE") + "\n");
//                        bw.write("EPSILONRULE                 " + properties.getString("EPSILONRULE") + "\n");
//                        bw.write("VDW-14-SCALE                 " + properties.getString("VDW-14-SCALE") + "\n");
//                        bw.write("CHG-14-SCALE                 " + properties.getString("CHG-14-SCALE") + "\n");
//                        bw.write("ELECTRIC                 " + properties.getString("ELECTRIC") + "\n");
//                        bw.write("DIELECTRIC                 " + properties.getString("DIELECTRIC") + "\n");
//                        bw.write(forceField.toStringBuffer().toString());
//                        bw.close();
//                    } else {
//                        logger.warning(" Could not create Patch file \"" + ffOutFile.getName() + "\"");
//                    }

                    sb.append(sbVdw.toString());
                    sb.append(sbBond.toString());
                    sb.append(sbAngle.toString());
                    sb.append(sbTor.toString());
                    sb.append(sbAPC.toString());
                    //Finalize print statements.
                    if (logger.isLoggable(Level.FINER)) {
                        logger.finer(sb.toString());
                    }
                    File outFile = PDBFilter.version(new File(fileName + ".patch"));
                    if (outFile.createNewFile()) {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(outFile, false));
                        bw.write(sb.toString());
                        bw.close();
                    } else {
                        logger.warning(" Could not create Patch file \"" + outFile.getName() + "\"");
                    }

                    logger.info(" END OF PARSE");
                } catch (Exception ex) {
                    logger.info(" Exception thrown during XML parse.")
                    ex.printStackTrace()
                }
            }
            // Handle SDF --> XYZ conversion
            if (FilenameUtils.getExtension(filename).equalsIgnoreCase("SDF")) {
                // Create XYZ file for inputted SDF.
                try (BufferedReader br = new BufferedReader(new FileReader(new File(filename)))) {
                    String data = br.readLine()
                    if (data == null) {
                        return this
                    }
                    while (data != null) {
                        logger.info(data);
                        String[] tokens = data.trim().split(" +");
                        int numToken = tokens.length;
                        int numAtoms;
                        int numBonds;
                        if (numToken == 10) {
                            numAtoms = tokens[0].toInteger();
                            numBonds = tokens[1].toInteger();
                            logger.info(format(" numAtom: %2d numBond: %2d", numAtoms, numBonds));
                            MolecularAssembly molecularAssembly = new MolecularAssembly(FilenameUtils.getBaseName(filename));
                            Molecule molecule = new Molecule(FilenameUtils.getBaseName(filename));
                            List<Atom> atomList = new ArrayList<>(numAtoms);
                            for (int i = 0; i < numAtoms; i++) {
                                data = br.readLine();
                                logger.info(data);
                                tokens = data.trim().split(" +");
                                String name = format("%s%d", tokens[3], i);
                                double[] coords = new double[]{tokens[0].toDouble(), tokens[1].toDouble(), tokens[2].toDouble()}
                                int valence = 0;
                                double[] atomicMass = AtomType.atomicMass
                                int atomicNumber = ElementSymbol.valueOf(tokens[3]).ordinal();
                                double atomicWeight = atomicMass[atomicNumber];
                                AtomType atomType = new AtomType(i+1, 0, name, "SDF", atomicNumber, atomicWeight, valence);
                                Atom atom = new Atom(i+1, name, atomType, coords);
                                atomList.add(atom)
                            }
                            List<Bond> bondList = new ArrayList<>(numBonds);
                            for (int i = 0; i < numBonds; i++) {
                                data = br.readLine();
                                tokens = data.trim().split(" +");
                                Atom newA1 = atomList.get(tokens[0].toInteger() - 1);
                                Atom newA2 = atomList.get(tokens[1].toInteger() - 1);
                                Bond bond = new Bond(newA1, newA2);
                                bondList.add(bond)
                            }
                            molecularAssembly.add(molecule);
                            Utilities.biochemistry(molecularAssembly, atomList);
                            XYZFilter xyzFilter = new XYZFilter(new File(FilenameUtils.getPath(filename) + FilenameUtils.getBaseName(filename) + ".xyz"), molecularAssembly, null, null);
                            xyzFilter.writeFile(xyzFilter.getFile(), true);
                        }
                        data = br.readLine()
                    }
                } catch (Exception ex) {
                    logger.info(" Exception thrown during target file parse.")
                    ex.printStackTrace()
                }
                // Write out a property file.
                File propertyFile = new File(FilenameUtils.getPath(filename) + FilenameUtils.getBaseName(filename) + ".properties");
                if (!propertyFile.exists()) {
                    try(BufferedWriter bw = new BufferedWriter(new FileWriter(propertyFile, false))) {
                        bw.write(format("parameters %s\n", FilenameUtils.getBaseName(filename)+".patch"));
                        logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath() + "\n");
                        bw.close();
                    } catch (Exception ex) {
                        logger.info("Failed to write files.\n" + ex);
                    }
                } else {
                    logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath() + "\n");
                }
            }
        }
        return this;
    }
}