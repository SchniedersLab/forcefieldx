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

import static java.lang.Integer.toString
import static java.lang.String.format
import static ffx.potential.bonded.Atom.ElementSymbol;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.PI;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;


/**
 * The ReadXML script converts an Open Force Field XML parameter file to Patch format.
 *
 * @author Jacob M. Miller
 * <br>
 * Usage:
 * <br>
 * ffxc test.ReadXML &lt;filename&gt;
 */
@Command(description = " Read OpenMM XML file.", name = "ffxc test.ReadXML")
class ReadXML extends PotentialScript {

    /**
     * The final argument(s) should be two filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'XML FF file to be read.')
    List<String> filenames = null

    final double ANGperNM = 10.0;
    final double KJperKCal = 4.184;

    LinkedHashMap<String, Integer> atomClassMap
    int atomClassNum = 0

    /**
     * Execute the script.
     */
    @Override
    ReadXML run() {

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        CompositeConfiguration properties = Keyword.loadProperties(null);
        properties.addProperty("renumberPatch", "FALSE")
        properties.addProperty("FORCEFIELD", "AMBER_1999_SB")
        properties.addProperty("VDWTYPE", "LENNARD-JONES")
        properties.addProperty("RADIUSRULE", "ARITHMETIC")
        properties.addProperty("RADIUSTYPE", "R-MIN")
//        properties.addProperty("RADIUSTYPE", "SIGMA")
        properties.addProperty("RADIUSSIZE", "RADIUS")
//        properties.addProperty("RADIUSSIZE", "DIAMETER")
        properties.addProperty("EPSILONRULE", "GEOMETRIC")
        properties.addProperty("VDW-14-SCALE", "2.0")
//        properties.addProperty("VDW-14-SCALE", "0.500000")
        properties.addProperty("CHG-14-SCALE", "1.2")
        properties.addProperty("ELECTRIC", "332.0522173") //was 332.0716
        properties.addProperty("DIELECTRIC", "1.0")
        ForceField forceField = new ForceField(properties);

        atomClassMap = new LinkedHashMap<>()
        File inputFile = new File(filenames[0])
        String fileName = FilenameUtils.removeExtension(inputFile.getName())
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance()
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder()
        Document doc = dBuilder.parse(inputFile)

        logger.info(format("Filename: %s%n",fileName))
        // Node "ForceField"
        logger.info(format(" Root Element: %s", doc.getDocumentElement().getNodeName()))

        NodeList nodeList = doc.getChildNodes()
        Node node = nodeList.item(0) // Assumed one system label for now (ForceField)
        NodeList childNodes = node.getChildNodes()
        logger.info(" Child Size:" + childNodes.length)

        for (Node child : childNodes) {
            if (child.hasChildNodes()) {
                switch (child.getNodeName()) {
                    case "AtomTypes":
                        NodeList types = child.getChildNodes()
                        logger.info(format("AtomTypes nodes: %d",types.length))
//                        int atomClassNum = 0
                        for (Node atom : types) {
                            if (atom.getNodeName() == "Type") {
                                logger.info(format("%s %s %s %s",atom.getAttribute("name"),atom.getAttribute("class"),atom.getAttribute("element"),atom.getAttribute("mass")))

                                int atomType = parseInt(atom.getAttribute("name"))
                                String className = atom.getAttribute("class")
                                String element = atom.getAttribute("element")
                                double mass = parseDouble(atom.getAttribute("mass"))

                                if (!atomClassMap.containsKey(className)) {
                                    atomClassMap.put(className, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",className,atomClassNum))
                                    atomClassNum++
                                }

                                int atomClass = atomClassMap.get(className)
                                int atomicNumber = ElementSymbol.valueOf(element).ordinal()+1
//                                logger.info(format("ELEMENT: %s %d",element,ElementSymbol.valueOf(element).ordinal()+1))

                                //TODO: get valence -> go thru residues and count up the number of bonds for each atom
                                // type.. will likely now have to add the atomtype to the forcefield outside of this loop
//                                forceField.addForceFieldType(new AtomType(atomType, atomClass, element, className, atomicNumber, mass, VALENCE))

                            } else if (atom.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "Residues":
                        NodeList residues = child.getChildNodes()
                        logger.info(format("Residues nodes: %d",residues.length))

                        for (Node res : residues) {
                            if (res.hasChildNodes()) {
                                logger.info(format("Residue: %s",res.getNodeName()))
                                NodeList resProps = res.getChildNodes()

                                for (Node resProp : resProps) {
                                    if (resProp.getNodeName() == "Atom") {
                                        logger.info(format("    Atom: %s %s",resProp.getAttribute("name"),resProp.getAttribute("type")))
                                    } else if (resProp.getNodeName() == "Bond") {
                                        logger.info(format("    Bond: %s %s",resProp.getAttribute("from"),resProp.getAttribute("to")))
                                    } else if (resProp.getNodeName() == "ExternalBond") {
                                        logger.info(format("    ExternalBond: %s",resProp.getAttribute("from")))
                                    } else if (resProp.hasAttributes()) {
                                        logger.info("CHECK")
                                    }
                                }
                            }
                        }
                        break

                    case "HarmonicBondForce":
                        NodeList bonds = child.getChildNodes()
                        logger.info(format("HarmonicBondForce nodes: %d",bonds.length))

                        for (Node bond : bonds) {
                            if (bond.getNodeName() == "Bond") {
                                logger.info(format("%s %s %s %s",bond.getAttribute("class1"),bond.getAttribute("class2"),bond.getAttribute("length"),bond.getAttribute("k")))
                                String class1 = bond.getAttribute("class1")
                                String class2 = bond.getAttribute("class2")
                                String bondLength = bond.getAttribute("length")
                                String k = bond.getAttribute("k")
//                                int[] classes = {atomClassMap.get(class1); atomClassMap.get(class2)}
//                                logger.info(format("BOND %s - %s (%s - %s",class1,class2,classes[0],classes[1]))
                                if (!atomClassMap.containsKey(class1)) {
                                    atomClassMap.put(class1, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",class1,atomClassNum))
                                    atomClassNum++
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",class2,atomClassNum))
                                    atomClassNum++
                                }

                                int[] classes = new int[2]
                                classes[0] = atomClassMap.get(class1)
                                classes[1] = atomClassMap.get(class2)
//                                logger.info(format("BOND %s - %s (%s - %s",class1,class2,class1Int,class2Int))
                                forceField.addForceFieldType(new BondType(classes, parseDouble(k), parseDouble(bondLength)))
                            } else if (bond.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "HarmonicAngleForce":
                        NodeList angles = child.getChildNodes()
                        logger.info(format("HarmonicAngleForce nodes: %d",angles.length))

                        for (Node angle : angles) {
                            if (angle.getNodeName() == "Angle") {
                                logger.info(format("%s %s %s %s %s",angle.getAttribute("class1"),angle.getAttribute("class2"),
                                        angle.getAttribute("class3"),angle.getAttribute("angle"),angle.getAttribute("k")))
                                String class1 = angle.getAttribute("class1")
                                String class2 = angle.getAttribute("class2")
                                String class3 = angle.getAttribute("class3")
                                String angleStr = angle.getAttribute("angle")
                                String k = angle.getAttribute("k")

                                if (!atomClassMap.containsKey(class1)) {
                                    atomClassMap.put(class1, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",class1,atomClassNum))
                                    atomClassNum++
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",class2,atomClassNum))
                                    atomClassNum++
                                }
                                if (!atomClassMap.containsKey(class3)) {
                                    atomClassMap.put(class3, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",class3,atomClassNum))
                                    atomClassNum++
                                }

                                int[] classes = new int[3]
                                classes[0] = atomClassMap.get(class1)
                                classes[1] = atomClassMap.get(class2)
                                classes[2] = atomClassMap.get(class3)
                                double[] angleArr = new double[1]
                                angleArr[0] = parseDouble(angleStr) * 180.0 / PI  // convert to degrees

                                forceField.addForceFieldType(new AngleType(classes, parseDouble(k), angleArr))
                            } else if (angle.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "PeriodicTorsionForce":
                        NodeList torsions = child.getChildNodes()
                        logger.info(format("PeriodicTorsionForce nodes: %d",torsions.length))

                        for (Node torsion : torsions) {
                            if (torsion.getNodeName() == "Proper" || torsion.getNodeName() == "Improper") {
                                logger.info(format("%s %s %s %s %s %s %s",torsion.getAttribute("class1"),
                                        torsion.getAttribute("class2"),torsion.getAttribute("class3"),
                                        torsion.getAttribute("class4"),torsion.getAttribute("periodicity1"),
                                        torsion.getAttribute("phase1"),torsion.getAttribute("k1")))
                                if (torsion.hasAttribute("periodicity2")) {
                                    logger.info(format("%s %s %s",torsion.getAttribute("periodicity2"),torsion.getAttribute("phase2"),torsion.getAttribute("k2")))
                                    if(torsion.hasAttribute("periodicity3")) {
                                        logger.info(format("%s %s %s",torsion.getAttribute("periodicity3"),torsion.getAttribute("phase3"),torsion.getAttribute("k3")))
                                        if(torsion.hasAttribute("periodicity4")) {
                                            logger.info(format("%s %s %s",torsion.getAttribute("periodicity4"),torsion.getAttribute("phase4"),torsion.getAttribute("k4")))
                                        }
                                    }
                                }

                                int numTerms = (torsion.getAttributes().length - 4) / 3  // get number of amplitudes/phases/periodicities
                                int[] classes = new int[4]
                                int[] periods = new int[numTerms]
                                double[] phases = new double[numTerms]
                                double[] amplitudes = new double[numTerms]
                                for (int i = 1; i <= 4; i++) {
                                    if (!atomClassMap.containsKey(torsion.getAttribute("class"+toString(i)))) {
                                        atomClassMap.put(torsion.getAttribute("class"+toString(i)), atomClassNum)
                                        logger.info(format("KEY: %s     VALUE: %d",torsion.getAttribute("class"+toString(i)),atomClassNum))
                                        atomClassNum++
                                    }
                                    classes[i-1] = atomClassMap.get(torsion.getAttribute("class"+toString(i)))

                                    if (torsion.hasAttribute("periodicity"+toString(i))) { // could replace with i < numTerms
                                        periods[i-1] = parseInt(torsion.getAttribute("periodicity"+toString(i)))
                                        phases[i-1] = parseDouble(torsion.getAttribute("phase"+toString(i))) * 180.0 / PI  // convert to degrees
                                        amplitudes[i-1] = parseDouble(torsion.getAttribute("k"+toString(i)))
                                    }
                                }

                                if (torsion.getNodeName() == "Proper") {
                                    forceField.addForceFieldType(new TorsionType(classes, amplitudes, phases, periods, TorsionType.TorsionMode.NORMAL))
                                } else {
                                    forceField.addForceFieldType(new TorsionType(classes, amplitudes, phases, periods, TorsionType.TorsionMode.IMPROPER))
                                }
                            }
                            else if (torsion.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "NonbondedForce":
                        NodeList nbForces = child.getChildNodes()
                        logger.info(format("NonbondedForce nodes: %d",nbForces.length))

                        for (Node nbF : nbForces) {
                            if (nbF.getNodeName() == "Atom") {
                                logger.info(format("%s %s %s %s",nbF.getAttribute("type"),nbF.getAttribute("charge"),nbF.getAttribute("sigma"),nbF.getAttribute("epsilon")))
                                int atomType = parseInt(nbF.getAttribute("type"))
                                double q = parseDouble(nbF.getAttribute("charge")) // in proton units
                                double sigma = parseDouble(nbF.getAttribute("sigma")) * ANGperNM  // nm to Ang
                                double eps = parseDouble(nbF.getAttribute("epsilon")) / KJperKCal  // kJ/mol? to KCal/mol

                                forceField.addForceFieldType(new ChargeType(atomType, q))
                                forceField.addForceFieldType(new VDWType(atomType, sigma, eps, -1.0))

                            } else if (nbF.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break
                }
            }
        }

        return this
    }
}