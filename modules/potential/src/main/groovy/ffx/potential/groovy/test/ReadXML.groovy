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
import ffx.potential.parameters.BioType
import ffx.potential.parameters.BondType
import ffx.potential.parameters.ChargeType
import ffx.potential.parameters.ForceField
import ffx.potential.parameters.TorsionType
import ffx.potential.parameters.VDWType
import ffx.potential.parsers.ForceFieldFilter
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.XYZFilter
import ffx.utilities.FFXKeyword
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
    LinkedHashMap<String, String> biotypeMap
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

        CompositeConfiguration properties = Keyword.loadProperties(null)
        properties.addProperty("renumberPatch", "FALSE")
        properties.addProperty("FORCEFIELD", "AMBER_1999_SB_XML")
        properties.addProperty("VDWINDEX", "TYPE")
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
        ForceField forceField = new ForceField(properties)

        atomClassMap = new LinkedHashMap<>()
        biotypeMap = new LinkedHashMap<>()

        File inputFile = new File(filenames[0])
        String fileName = FilenameUtils.removeExtension(inputFile.getName())
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance()
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder()
        Document doc = dBuilder.parse(inputFile)

        logger.info(format(" Filename: %s%n",fileName))
        logger.info(format(" Root Element: %s", doc.getDocumentElement().getNodeName())) // "ForceField"

        File biotypeClasses = new File(filenames[1])
        Scanner myReader = new Scanner(biotypeClasses)
        while (myReader.hasNextLine()) {
            String data = myReader.nextLine()
            String[] columnSplit = data.split("\t",-1)
            String str = columnSplit[0].replace(',', ' ').replace('"', ' ').trim()
            biotypeMap.put(str,columnSplit[1])
        }
        myReader.close()

        NodeList nodeList = doc.getChildNodes()
        Node node = nodeList.item(0) // Assumed one system label for now (ForceField)
        NodeList childNodes = node.getChildNodes()
//        logger.info(" Child Size:" + childNodes.length)

        int numAtomTypes
        int[] atomTypes
        int[] atomClasses
        String[] atomNames
        String[] atomEnvs // get from which residue it's in
        int[] atomicNumbers
        double[] atomicWeights
        int[] valence
        String[] biotypeAtomNames

        //TODO: change indices to start at 1 for TINKER format
        for (Node child : childNodes) {
            if (child.hasChildNodes()) {
                switch (child.getNodeName()) {
                    case "AtomTypes":
                        Element e = (Element) child;
                        numAtomTypes = e.getElementsByTagName("Type").length
                        atomTypes = new int[numAtomTypes]
                        atomClasses = new int[numAtomTypes]
                        atomNames = new String[numAtomTypes]
                        atomEnvs = new String[numAtomTypes]
                        atomicNumbers = new int[numAtomTypes]
                        atomicWeights = new int[numAtomTypes]
                        valence = new int[numAtomTypes]
                        biotypeAtomNames = new String[numAtomTypes]

                        NodeList types = child.getChildNodes()
//                        logger.info(format("AtomTypes nodes: %d",types.length))
                        int idx = 0
                        for (Node atom : types) {
                            if (atom.getNodeName() == "Type") {
//                                logger.info(format("%s %s %s %s",atom.getAttribute("name"),atom.getAttribute("class"),atom.getAttribute("element"),atom.getAttribute("mass")))
                                atomTypes[idx] = parseInt(atom.getAttribute("name"))
                                String className = atom.getAttribute("class")
                                String element = atom.getAttribute("element")
                                atomNames[idx] = className
                                atomicWeights[idx] = parseDouble(atom.getAttribute("mass"))

                                if (!atomClassMap.containsKey(className)) {
                                    atomClassMap.put(className, atomClassNum)
                                    logger.info(format("KEY: %s     VALUE: %d",className,atomClassNum))
                                    atomClassNum++
                                }

                                atomClasses[idx] = atomClassMap.get(className)
                                atomicNumbers[idx] = ElementSymbol.valueOf(element).ordinal()+1

                                idx++
//                                forceField.addForceFieldType(new AtomType(atomType, atomClass, element, className, atomicNumber, mass, VALENCE))

                            } else if (atom.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "Residues":
                        NodeList residues = child.getChildNodes()
//                        logger.info(format("Residues nodes: %d",residues.length))

                        for (Node res : residues) {
                            if (res.hasChildNodes()) {
                                String resName = res.getAttribute("name")

//                                logger.info(format("Residue: %s",res.getNodeName()))
                                NodeList resProps = res.getChildNodes()

                                Element e = (Element) res;
//                                logger.info(format("RES NODE LENGTH: %d",resProps.length))
//                                logger.info(format("TEST: %s",e.getElementsByTagName("Atom").length))
                                String[] nameMap = new String[e.getElementsByTagName("Atom").length]
                                int[] typeMap = new int[e.getElementsByTagName("Atom").length]
                                int[] bondCount = new int[e.getElementsByTagName("Atom").length]
                                int i = 0
                                for (Node resProp : resProps) {
                                    if (resProp.getNodeName() == "Atom") {
//                                        logger.info(format("    Atom: %s %s",resProp.getAttribute("name"),resProp.getAttribute("type")))
                                        String atomName = resProp.getAttribute("name")
                                        int atomType = parseInt(resProp.getAttribute("type"))
                                        nameMap[i] = atomName
                                        typeMap[i] = atomType
                                        i++

                                    } else if (resProp.getNodeName() == "Bond") {
//                                        logger.info(format("    Bond: %s %s",resProp.getAttribute("from"),resProp.getAttribute("to")))
                                        int atom1 = parseInt(resProp.getAttribute("from"))
                                        int atom2 = parseInt(resProp.getAttribute("to"))
                                        bondCount[atom1]++
                                        bondCount[atom2]++

                                    } else if (resProp.getNodeName() == "ExternalBond") {
//                                        logger.info(format("    ExternalBond: %s",resProp.getAttribute("from")))
                                        int extBondAtom = parseInt(resProp.getAttribute("from"))
                                        bondCount[extBondAtom]++ // could just use if statement above and get rid of this

                                    } else if (resProp.hasAttributes()) {
                                        logger.info("CHECK")
                                    }
                                }
                                for (int j = 0; j < typeMap.length; j++) {
//                                    logger.info(format("RES NODE %s: %d %d",res.getNodeName(),typeMap[j],bondCount[j]))
                                    valence[typeMap[j]] = bondCount[j] //TODO could check to see if it equals the previous value in the valence array
                                    if (atomEnvs[typeMap[j]] == "" || atomEnvs[typeMap[j]] == null) {
                                        atomEnvs[typeMap[j]] = resName
                                        biotypeAtomNames[typeMap[j]] = nameMap[j]
                                    } else if (atomEnvs[typeMap[j]] != resName || biotypeAtomNames[typeMap[j]] != nameMap[j]) {
                                        atomEnvs[typeMap[j]] = atomEnvs[typeMap[j]] + ";" + resName
                                        biotypeAtomNames[typeMap[j]] = biotypeAtomNames[typeMap[j]] + ";" + nameMap[j]
//                                        logger.info(format("ALREADY FILLED ORG: AtomType: %d Res: %s Name: %s",typeMap[j],atomEnvs[typeMap[j]],biotypeAtomNames[typeMap[j]]))
//                                        logger.info(format("ALREADY FILLED NEW: AtomType: %d Res: %s Name: %s",typeMap[j],resName,nameMap[j]))
                                        //TODO:
                                    }
                                }
                            }
                        }
                        break

                    case "HarmonicBondForce":
                        NodeList bonds = child.getChildNodes()
//                        logger.info(format("HarmonicBondForce nodes: %d",bonds.length))

                        for (Node bond : bonds) {
                            if (bond.getNodeName() == "Bond") {
//                                logger.info(format("%s %s %s %s",bond.getAttribute("class1"),bond.getAttribute("class2"),bond.getAttribute("length"),bond.getAttribute("k")))
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

                                double forceConstant = parseDouble(k) / KJperKCal / (ANGperNM * ANGperNM) / 2 // TODO why divided by 2?
                                double distance = parseDouble(bondLength) * ANGperNM

//                                logger.info(format("BOND %s - %s (%s - %s",class1,class2,class1Int,class2Int))
                                forceField.addForceFieldType(new BondType(classes, forceConstant, distance))
                            } else if (bond.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "HarmonicAngleForce":
                        NodeList angles = child.getChildNodes()
//                        logger.info(format("HarmonicAngleForce nodes: %d",angles.length))

                        for (Node angle : angles) {
                            if (angle.getNodeName() == "Angle") {
//                                logger.info(format("%s %s %s %s %s",angle.getAttribute("class1"),angle.getAttribute("class2"),
//                                        angle.getAttribute("class3"),angle.getAttribute("angle"),angle.getAttribute("k")))
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

                                double forceConstant = parseDouble(k) / KJperKCal / 2 //TODO why divided by two?

                                forceField.addForceFieldType(new AngleType(classes, forceConstant, angleArr))
                            } else if (angle.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "PeriodicTorsionForce":
                        NodeList torsions = child.getChildNodes()
//                        logger.info(format("PeriodicTorsionForce nodes: %d",torsions.length))

                        for (Node torsion : torsions) {
                            if (torsion.getNodeName() == "Proper" || torsion.getNodeName() == "Improper") {
                                /*
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
                                 */

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
                                        amplitudes[i-1] = parseDouble(torsion.getAttribute("k"+toString(i))) / KJperKCal // convert to kcal/mol
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
//                        logger.info(format("NonbondedForce nodes: %d",nbForces.length))

                        for (Node nbF : nbForces) {
                            if (nbF.getNodeName() == "Atom") {
//                                logger.info(format("%s %s %s %s",nbF.getAttribute("type"),nbF.getAttribute("charge"),nbF.getAttribute("sigma"),nbF.getAttribute("epsilon")))
                                int atomType = parseInt(nbF.getAttribute("type"))
                                double q = parseDouble(nbF.getAttribute("charge")) // in proton units
                                double sigma = parseDouble(nbF.getAttribute("sigma")) * ANGperNM  // nm to Ang
                                double eps = parseDouble(nbF.getAttribute("epsilon")) / KJperKCal  // kJ/mol? to KCal/mol

                                forceField.addForceFieldType(new ChargeType(atomType, q))
                                forceField.addForceFieldType(new VDWType(atomType, sigma, eps, -1.0))  // vdw by atom type

//                                int atomClass = atomClasses[atomType]
//                                forceField.addForceFieldType(new VDWType(atomClass, sigma, eps, -1.0))
                            } else if (nbF.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break
                }
            }
        }
        for (int j = 0; j < numAtomTypes; j++) {
            forceField.addForceFieldType(new AtomType(atomTypes[j], atomClasses[j], atomNames[j], atomEnvs[j], atomicNumbers[j], atomicWeights[j], valence[j]))
        }

        CompositeConfiguration props = Keyword.loadProperties(null)
        props.setProperty("parameters", "/iahome/j/jm/jmiller99/forcefieldx/modules/potential/src/main/java/ffx/potential/parameters/ff/AMBER_1999_SB")
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(props)
        ForceField predefFF = forceFieldFilter.parse()
        Map<String, BioType> biotypes = predefFF.getBioTypeMap()

        for (BioType biotype : biotypes.values()) {
            String resID = biotypeMap.get(biotype.moleculeName)
            logger.info(format("MoleculeName: %s -> resName: %s",biotype.moleculeName,resID))

            String btAtomName = biotype.atomName
            int tempFlag = 0
            for (int idx = 0; idx < numAtomTypes; idx++) {
                String[] diffRes = atomEnvs[idx].split(";")
                for (int j = 0; j < diffRes.length; j++) {
                    if (diffRes[j] == resID) {
                        String xmlAtomName = biotypeAtomNames[idx].split(";")[j] // split string at ';' and grab element corresponding to the correct residue
                        xmlAtomName = xmlAtomName.replace("'",'*') // biotype * = xml

                        if (btAtomName == xmlAtomName) {
                            // Add new BioType to forcefield with correct atom type corresponding to the xml
                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            //add continue or break here?
                        } else if (btAtomName + "1" == xmlAtomName || btAtomName + "2" == xmlAtomName || btAtomName + "3" == xmlAtomName) {
                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            // could compare atomtypes between the ones with numbers -> not best setup
                            // TODO: could add if statement to compare current atomType to new -> if they match fine, if not LOG
                        } else if (btAtomName == "HN") {
                            if (xmlAtomName == "H" || xmlAtomName == "H1" || xmlAtomName == "H2" || xmlAtomName == "H3") {
//                                if (xmlAtomName == "H" || xmlAtomName == "HN1" || xmlAtomName == "HN2") {
                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (biotype.moleculeName.contains(' Ion')) {
                            if (btAtomName + "+" == xmlAtomName || btAtomName + "-" == xmlAtomName) {
                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            } else if (btAtomName.length() == 2) {
                                String ionName = btAtomName[0] + btAtomName.toLowerCase()[1]
                                if (ionName + "+" == xmlAtomName || ionName + "-" == xmlAtomName) {
                                    logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                    forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                    tempFlag = 1
                                }
                            }
                        } else if (btAtomName + "11" == xmlAtomName || btAtomName + "12" == xmlAtomName || btAtomName + "13" == xmlAtomName || btAtomName + "21" == xmlAtomName || btAtomName + "22" == xmlAtomName) {
                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                        } else if (btAtomName == "H") {
//                            if (btAtomName + "H" == xmlAtomName.substring(0, 2)) {
                            if (btAtomName + "H31" == xmlAtomName || btAtomName + "H32" == xmlAtomName || btAtomName + "H33" == xmlAtomName) {
                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (btAtomName == "OP") {
                            if (xmlAtomName == "O1P" || xmlAtomName == "O2P") {
                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        }
                    }
                }
            }
            if (tempFlag == 0) {
                logger.info(format("NOTFOUND: biotype %d %s %s %d",biotype.index,biotype.atomName,biotype.moleculeName,biotype.atomType))
                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, -1)) // not in XML FF
            }
        }

        StringBuffer ffSB = forceField.toStringBuffer()

        BufferedWriter out = new BufferedWriter(new FileWriter("sampleFF.txt"))
        out.write(ffSB.toString())
        out.flush()
        out.close()

        return this
    }
}