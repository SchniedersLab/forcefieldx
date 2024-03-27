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

import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.*
import ffx.potential.parsers.ForceFieldFilter
import ffx.utilities.Keyword
import org.apache.commons.configuration2.CompositeConfiguration
import org.w3c.dom.Document
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.w3c.dom.NodeList
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import javax.xml.parsers.DocumentBuilder
import javax.xml.parsers.DocumentBuilderFactory

import static ffx.potential.bonded.Atom.ElementSymbol
import static java.lang.Double.parseDouble
import static java.lang.Integer.parseInt
import static java.lang.Integer.toString
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.PI
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

    LinkedHashMap<String, Integer> atomTypeMap
    LinkedHashMap<String, Integer> atomClassMap
    LinkedHashMap<String, String> biotypeMap

    /**
     * Execute the script.
     */
    @Override
    ReadXML run() {

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        atomTypeMap = new LinkedHashMap<>()
        atomClassMap = new LinkedHashMap<>()
        biotypeMap = new LinkedHashMap<>()

        // specify files coming in
        File biotypeClasses = new File(filenames[0])  // biotype description to xml Residue Map TODO - could we use residuesFinal.xml
        File inputFile = new File(filenames[1]) // main xml FF file
        File watFile = new File(filenames[2]) // water xml FF file (e.g. TIP3P)

        // Read and create hash map of biotype molecule descriptions and XML residue names
        Scanner myReader = new Scanner(biotypeClasses)
        while (myReader.hasNextLine()) {
            String data = myReader.nextLine()
            String[] columnSplit = data.split("\t",-1)
            String str = columnSplit[0].replace(',', ' ').replace('"', ' ').trim()
            biotypeMap.put(str,columnSplit[1])
        }
        myReader.close()

        // Instantiate Document building objects
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = null;
        Document doc = null;
        Document watDoc = null;

        // Parse XMLs into Document objects
        dBuilder = dbFactory.newDocumentBuilder()
        doc = dBuilder.parse(inputFile)
        watDoc = dBuilder.parse(watFile)

//        readAmber(doc, watDoc)
        readCharmm(doc, watDoc)

        return this
    }

    private void readAmber(Document doc, Document watDoc) {
        CompositeConfiguration properties = Keyword.loadProperties(null)
//        properties.addProperty("renumberPatch", "FALSE")
        properties.addProperty("FORCEFIELD", "AMBER_1999_SB_XML")
        properties.addProperty("VDWINDEX", "TYPE")
        properties.addProperty("VDWTYPE", "LENNARD-JONES")
        properties.addProperty("RADIUSRULE", "ARITHMETIC")
        properties.addProperty("RADIUSTYPE", "R-MIN")
//        properties.addProperty("RADIUSTYPE", "SIGMA")
        properties.addProperty("RADIUSSIZE", "RADIUS")
//        properties.addProperty("RADIUSSIZE", "DIAMETER")
        properties.addProperty("EPSILONRULE", "GEOMETRIC")
//        properties.addProperty("VDW-14-SCALE", "2.0")
        properties.addProperty("VDW-14-SCALE", "0.500000")
//        properties.addProperty("CHG-14-SCALE", "1.2")
        properties.addProperty("CHG-14-SCALE", "0.833333")
        properties.addProperty("ELECTRIC", "332.0522173") //was 332.0716
        properties.addProperty("DIELECTRIC", "1.0")
        ForceField forceField = new ForceField(properties)

        String[] xmlNodes = ["AtomTypes", "Residues", "HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce"]

        // Instantiate Atom Type variables
        int numAtomTypes
        int[] atomTypes
        int[] atomClasses
        String[] atomNames
        String[] atomEnvs // get from which residue it's in
        int[] atomicNumbers
        double[] atomicWeights
        int[] valence
        String[] biotypeAtomNames

        // Loop through all types of XML ForceField Nodes
        for (String xmlNode : xmlNodes) {
            NodeList parentNode = doc.getElementsByTagName(xmlNode) // grabs a parent node specified by the string

            if (parentNode.item(0).hasChildNodes()) {
                switch (parentNode.item(0).getNodeName()) {
                    case "AtomTypes":
                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("Type").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Type").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        Element e = (Element) parentNode.item(0);
                        numAtomTypes = e.getElementsByTagName("Type").length
                        atomTypes = new int[numAtomTypes]
                        atomClasses = new int[numAtomTypes]
                        atomNames = new String[numAtomTypes]
                        atomEnvs = new String[numAtomTypes]
                        atomicNumbers = new int[numAtomTypes]
                        atomicWeights = new int[numAtomTypes]
                        valence = new int[numAtomTypes]
                        biotypeAtomNames = new String[numAtomTypes]

                        NodeList types = parentNode.item(0).getChildNodes()
                        int idx = 0
                        for (Node atom : types) {
                            if (atom.getNodeName() == "Type") {
                                String atomType = atom.getAttribute("name")
                                if (!atomTypeMap.containsKey(atomType)) {
                                    atomTypeMap.put(atomType,atomTypeMap.size()+1)
                                }
                                atomTypes[idx] = atomTypeMap.get(atomType)
                                String className = atom.getAttribute("class")
                                String element = atom.getAttribute("element")
                                atomNames[idx] = className
                                atomicWeights[idx] = parseDouble(atom.getAttribute("mass"))

                                if (!atomClassMap.containsKey(className)) {
                                    atomClassMap.put(className, atomClassMap.size()+1)
                                }

                                atomClasses[idx] = atomClassMap.get(className)
                                atomicNumbers[idx] = ElementSymbol.valueOf(element).ordinal()+1

                                idx++

                            } else if (atom.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "Residues":
                        // Import additional (water) atom types
                        // Since "Deep" is true, it keeps all child nodes of the incoming Residue node
                        for (int i = 0; i < watDoc.getElementsByTagName("Residue").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Residue").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList residues = parentNode.item(0).getChildNodes()
                        for (Node res : residues) {
                            if (res.hasChildNodes()) {
                                String resName = res.getAttribute("name")

                                NodeList resProps = res.getChildNodes()

                                Element e = (Element) res;
                                String[] nameMap = new String[e.getElementsByTagName("Atom").length]
                                int[] typeMap = new int[e.getElementsByTagName("Atom").length]
                                int[] bondCount = new int[e.getElementsByTagName("Atom").length]
                                int i = 0
                                for (Node resProp : resProps) {
                                    if (resProp.getNodeName() == "Atom") {
                                        String atomName = resProp.getAttribute("name")
                                        int atomType = atomTypeMap.get(resProp.getAttribute("type")) - 1 // just used as index, so need this
                                        nameMap[i] = atomName
                                        typeMap[i] = atomType
                                        i++

                                    } else if (resProp.getNodeName() == "Bond") {
                                        int atom1 = -1
                                        int atom2 = -1
                                        if (resProp.hasAttribute("from")) {
                                            atom1 = parseInt(resProp.getAttribute("from"))
                                            atom2 = parseInt(resProp.getAttribute("to"))
                                        } else if (resProp.hasAttribute("atomName1")) {
                                            List<String> nameList = Arrays.asList(nameMap) // could just start with it as a List..
                                            atom1 = nameList.indexOf(resProp.getAttribute("atomName1"))
                                            atom2 = nameList.indexOf(resProp.getAttribute("atomName2"))
                                        }
                                        bondCount[atom1]++
                                        bondCount[atom2]++

                                    } else if (resProp.getNodeName() == "ExternalBond") {
                                        int extBondAtom = parseInt(resProp.getAttribute("from"))
                                        bondCount[extBondAtom]++ // could just use if statement above and get rid of this

                                    } else if (resProp.hasAttributes()) {
                                        logger.info("CHECK")
                                    }
                                }
                                for (int j = 0; j < typeMap.length; j++) {
                                    valence[typeMap[j]] = bondCount[j] //TODO could check to see if it equals the previous value in the valence array
                                    if (atomEnvs[typeMap[j]] == "" || atomEnvs[typeMap[j]] == null) {
                                        atomEnvs[typeMap[j]] = resName
                                        biotypeAtomNames[typeMap[j]] = nameMap[j]
                                    } else if (atomEnvs[typeMap[j]] != resName || biotypeAtomNames[typeMap[j]] != nameMap[j]) {
                                        atomEnvs[typeMap[j]] = atomEnvs[typeMap[j]] + ";" + resName
                                        biotypeAtomNames[typeMap[j]] = biotypeAtomNames[typeMap[j]] + ";" + nameMap[j]
                                    }
                                }
                            }
                        }
                        break

                    case "HarmonicBondForce":
                        // Import additional (water) atom types
                        // Created element to only grab "Bond" nodes within the "HarmonicBondForce" node (there are also Bond nodes in Residues)
                        Element ele = (Element) watDoc.getElementsByTagName("HarmonicBondForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Bond").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Bond").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        NodeList bonds = parentNode.item(0).getChildNodes()

                        for (Node bond : bonds) {
                            if (bond.getNodeName() == "Bond") {
                                String class1 = bond.getAttribute("class1")
                                String class2 = bond.getAttribute("class2")
                                String bondLength = bond.getAttribute("length")
                                String k = bond.getAttribute("k")
                                if (!atomClassMap.containsKey(class1)) {
                                    atomClassMap.put(class1, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassMap.size()+1)
                                }

                                int[] classes = new int[2]
                                classes[0] = atomClassMap.get(class1)
                                classes[1] = atomClassMap.get(class2)

                                double forceConstant = parseDouble(k) / KJperKCal / (ANGperNM * ANGperNM) / 2 // TODO why divided by 2?
                                double distance = parseDouble(bondLength) * ANGperNM

                                forceField.addForceFieldType(new BondType(classes, forceConstant, distance))
                            } else if (bond.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "HarmonicAngleForce":
                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("Angle").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Angle").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        NodeList angles = parentNode.item(0).getChildNodes()

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
                                    atomClassMap.put(class1, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class3)) {
                                    atomClassMap.put(class3, atomClassMap.size()+1)
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
                        // Import additional (water) atom types
                        // Note: don't need this for TIP3P.xml, but here if we're combining other xml's
                        for (int i = 0; i < watDoc.getElementsByTagName("Proper").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Proper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        for (int i = 0; i < watDoc.getElementsByTagName("Imroper").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Imroper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList torsions = parentNode.item(0).getChildNodes()

                        for (Node torsion : torsions) {
                            if (torsion.getNodeName() == "Proper" || torsion.getNodeName() == "Improper") {

                                int numTerms = (torsion.getAttributes().length - 4) / 3  // get number of amplitudes/phases/periodicities
                                int[] classes = new int[4]
                                int[] periods = new int[numTerms]
                                double[] phases = new double[numTerms]
                                double[] amplitudes = new double[numTerms]
                                for (int i = 1; i <= 4; i++) {
                                    String torsClass= torsion.getAttribute("class"+toString(i))
                                    if (!atomClassMap.containsKey(torsClass)) {
                                        if (torsClass == "") {
                                            atomClassMap.put(torsClass, 0)
                                        } else {
                                            atomClassMap.put(torsClass, atomClassMap.size()+1)
                                        }
                                    }
                                    classes[i-1] = atomClassMap.get(torsClass)

                                    if (torsion.hasAttribute("periodicity"+toString(i))) { // could replace with i < numTerms
                                        periods[i-1] = parseInt(torsion.getAttribute("periodicity"+toString(i)))
                                        phases[i-1] = parseDouble(torsion.getAttribute("phase"+toString(i))) * 180.0 / PI  // convert to degrees
                                        amplitudes[i-1] = parseDouble(torsion.getAttribute("k"+toString(i))) / KJperKCal // convert to kcal/mol
                                    }
                                }

                                if (torsion.getNodeName() == "Proper") {
                                    forceField.addForceFieldType(new TorsionType(classes, amplitudes, phases, periods, TorsionType.TorsionMode.NORMAL))
                                } else {
                                    int trigAtom = classes[0] // Trigonal atom is defined first in XML, needs to be third in FFX/TINKER
                                    classes[0] = classes[1]
                                    classes[1] = classes[2]
                                    classes[2] = trigAtom
                                    forceField.addForceFieldType(new ImproperTorsionType(classes, amplitudes[0], phases[0], periods[0]))
                                }
                            }
                            else if (torsion.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "NonbondedForce":
                        // Import additional (water) atom types
                        // Created element to only grab "Atom" nodes within the "NonbondedForce" node (there are also Atom nodes in Residues)
                        Element ele = (Element) watDoc.getElementsByTagName("NonbondedForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Atom").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Atom").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList nbForces = parentNode.item(0).getChildNodes()

                        for (Node nbF : nbForces) {
                            if (nbF.getNodeName() == "Atom") {
                                int atomType = atomTypeMap.get(nbF.getAttribute("type"))
                                double q = parseDouble(nbF.getAttribute("charge")) // in proton units
                                double sigma = parseDouble(nbF.getAttribute("sigma")) * ANGperNM / 2 * Math.pow(2.0,1.0/6.0) // nm to Ang and divide by 2 to get radius & r-min = 2^(1/6) * sigma
                                double eps = parseDouble(nbF.getAttribute("epsilon")) / KJperKCal  // kJ/mol to KCal/mol

                                forceField.addForceFieldType(new ChargeType(atomType, q))
                                forceField.addForceFieldType(new VDWType(atomType, sigma, eps, -1.0))  // vdw by atom type
                            } else if (nbF.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break
                }
            }
        }
        for (int j = 0; j < numAtomTypes; j++) {
            forceField.addForceFieldType(new AtomType(atomTypes[j], atomClasses[j], atomNames[j], "\"" + atomEnvs[j] + "\"", atomicNumbers[j], atomicWeights[j], valence[j]))
        }

        // Log atom type Map
//        Set<String> atStrSet = atomTypeMap.keySet()
//        for (String str : atStrSet) {
//            logger.info(format("%s - %d",str,atomTypeMap.get(str)))
//        }
        // Log atom class Map
        Set<String> acStrSet = atomClassMap.keySet()
        for (String str : acStrSet) {
            logger.info(format("%s - %d",str,atomClassMap.get(str)))
        }

        // Read in current AMBER_1999_SB in FFX to use its biotypes to set the XML FF's biotypes
        CompositeConfiguration props = Keyword.loadProperties(null)
        props.setProperty("parameters", "/iahome/j/jm/jmiller99/forcefieldx/modules/potential/src/main/java/ffx/potential/parameters/ff/AMBER_1999_SB")
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(props)
        ForceField predefFF = forceFieldFilter.parse()
        Map<String, BioType> biotypes = predefFF.getBioTypeMap()

        for (BioType biotype : biotypes.values()) {
            String resID = biotypeMap.get(biotype.moleculeName)

            String btAtomName = biotype.atomName
            int tempFlag = 0
            for (int idx = 0; idx < numAtomTypes; idx++) {
                String[] diffRes = atomEnvs[idx].split(";") //TODO should not do this for all atoms .. or shouldnt do it multiple times
                for (int j = 0; j < diffRes.length; j++) {
                    if (diffRes[j] == resID) {
                        String xmlAtomName = biotypeAtomNames[idx].split(";")[j] // split string at ';' and grab element corresponding to the correct residue
                        xmlAtomName = xmlAtomName.replace("'",'*') // biotype * = xml

                        if (btAtomName == xmlAtomName) {
                            // Add new BioType to forcefield with correct atom type corresponding to the xml
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            //add continue or break here?
                        } else if (btAtomName + "1" == xmlAtomName || btAtomName + "2" == xmlAtomName || btAtomName + "3" == xmlAtomName) {
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            // could compare atomtypes between the ones with numbers -> not best setup
                            // TODO: could add if statement to compare current atomType to new -> if they match fine, if not LOG
                        } else if (btAtomName == "HN") {
                            if (xmlAtomName == "H" || xmlAtomName == "H1" || xmlAtomName == "H2" || xmlAtomName == "H3") {
//                                if (xmlAtomName == "H" || xmlAtomName == "HN1" || xmlAtomName == "HN2") {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (biotype.moleculeName.contains(' Ion')) {
                            if (btAtomName + "+" == xmlAtomName || btAtomName + "-" == xmlAtomName) {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            } else if (btAtomName.length() == 2) {
                                String ionName = btAtomName[0] + btAtomName.toLowerCase()[1]
                                if (ionName + "+" == xmlAtomName || ionName + "-" == xmlAtomName) {
//                                    logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                    forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                    tempFlag = 1
                                }
                            }
                        } else if (btAtomName + "11" == xmlAtomName || btAtomName + "12" == xmlAtomName || btAtomName + "13" == xmlAtomName || btAtomName + "21" == xmlAtomName || btAtomName + "22" == xmlAtomName) {
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                        } else if (btAtomName == "H") {
//                            if (btAtomName + "H" == xmlAtomName.substring(0, 2)) {
                            if (btAtomName + "H31" == xmlAtomName || btAtomName + "H32" == xmlAtomName || btAtomName + "H33" == xmlAtomName) {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (btAtomName == "OP") {
                            if (xmlAtomName == "O1P" || xmlAtomName == "O2P") {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
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

        // build header
        StringBuilder head = new StringBuilder()
//        head.append("forcefield\t\tAMBER-FF99SB-XML\n\n" +
//                    "vdwtype\t\t\tLENNARD-JONES\n" +
//                    "vdwindex\t\tTYPE\n" +
//                    "radiusrule\t\tARITHMETIC\n" +
//                    "radiustype\t\tR-MIN\n" + // or SIGMA
//                    "radiussize\t\tRADIUS\n" + // or DIAMETER
//                    "epsilonrule\t\tGEOMETRIC\n" +
//                    "vdw-14-scale\t\t0.500000\n" + // or 2.0
//                    "chg-14-scale\t\t0.833333\n" + // or 1.2
//                    "electric\t\t332.0522173\n" +
//                    "dielectric\t\t1.0\n")
        head.append("vdwtype LENNARD-JONES\n" +
                "vdwindex TYPE\n" +
                "radiusrule ARITHMETIC\n" +
                "radiustype R-MIN\n" + // or SIGMA
                "radiussize RADIUS\n" + // or DIAMETER
                "epsilonrule GEOMETRIC\n" +
                "vdw-14-scale 0.500000\n" + // or 2.0
                "chg-14-scale 0.833333\n" + // or 1.2
                "electric 332.0522173\n" +
                "dielectric 1.0\n")
// TODO: can get vdw-14-scale and chg-14-scale from the xml Node NonbondedForce attributes "lj14scale" and "coulomb14scale", respectively
        StringBuffer ffSB = forceField.toStringBuffer() // convert FF to string buffer
//        String ffStr = ffSB.replaceAll('improper','imptors') // replace improper FF terms with 'imptors'

        // write forcefield file
        BufferedWriter out = new BufferedWriter(new FileWriter("sampleFF.txt"))
        out.write(head.toString())
        out.write(ffSB.toString())
        out.flush()
        out.close()
    }

    private void readCharmm(Document doc, Document watDoc) {
        CompositeConfiguration properties = Keyword.loadProperties(null)
        properties.addProperty("FORCEFIELD", "CHARMM_36_CMAP_XML")
        properties.addProperty("VDWTYPE", "LENNARD-JONES")
        properties.addProperty("RADIUSRULE", "ARITHMETIC")
        properties.addProperty("RADIUSTYPE", "R-MIN")
        properties.addProperty("RADIUSSIZE", "RADIUS")
        properties.addProperty("EPSILONRULE", "GEOMETRIC")
        properties.addProperty("VDW-14-SCALE", "1.0")
        properties.addProperty("CHG-14-SCALE", "1.0")
        properties.addProperty("ELECTRIC", "332.0716")
        properties.addProperty("DIELECTRIC", "1.0")
        ForceField forceField = new ForceField(properties)

//        String[] xmlNodes = ["AtomTypes", "Residues", "Patches", "HarmonicBondForce", "HarmonicAngleForce",
//                             "AmoebaUreyBradleyForce", "PeriodicTorsionForce", "CustomTorsionForce",
//                             "CMAPTorsionForce", "NonbondedForce", "LennardJonesForce"]
        String[] xmlNodes = ["AtomTypes", "CMAPTorsionForce", "Residues"]

        // Instantiate Atom Type variables
        int numAtomTypes
        int[] atomTypes
        int[] atomClasses
        String[] atomNames
        String[] atomEnvs // get from which residue it's in
        int[] atomicNumbers
        double[] atomicWeights
        int[] valence
        String[] biotypeAtomNames

        // Loop through all types of XML ForceField Nodes
        for (String xmlNode : xmlNodes) {
            NodeList parentNode = doc.getElementsByTagName(xmlNode) // grabs a parent node specified by the string

            if (parentNode.item(0).hasChildNodes()) {
                switch (parentNode.item(0).getNodeName()) {
                    case "AtomTypes":
                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("Type").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Type").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        Element e = (Element) parentNode.item(0);
                        numAtomTypes = e.getElementsByTagName("Type").length
                        atomTypes = new int[numAtomTypes]
                        atomClasses = new int[numAtomTypes]
                        atomNames = new String[numAtomTypes]
                        atomEnvs = new String[numAtomTypes]
                        atomicNumbers = new int[numAtomTypes]
                        atomicWeights = new int[numAtomTypes]
                        valence = new int[numAtomTypes]
                        biotypeAtomNames = new String[numAtomTypes]
                        atomRes = new String[numAtomTypes]

                        NodeList types = parentNode.item(0).getChildNodes()
                        int idx = 0
                        for (Node atom : types) {
                            if (atom.getNodeName() == "Type") {
                                String name = atom.getAttribute("name")
                                if (!atomTypeMap.containsKey(name)) {
                                    atomTypeMap.put(name,atomTypeMap.size()+1)
                                }
                                atomTypes[idx] = atomTypeMap.get(name)
                                String className = atom.getAttribute("class")

                                String element
                                if (atom.hasAttribute("element")) {
                                    element = atom.getAttribute("element")
                                } else {
                                    element = null
                                }

                                atomNames[idx] = className
                                atomEnvs[idx] = name
                                atomicWeights[idx] = parseDouble(atom.getAttribute("mass"))

                                if (!atomClassMap.containsKey(className)) {
                                    atomClassMap.put(className, atomClassMap.size()+1)
                                }

                                atomClasses[idx] = atomClassMap.get(className)
                                if (element != null) {
                                    atomicNumbers[idx] = ElementSymbol.valueOf(element).ordinal()+1
                                } else {
                                    atomicNumbers[idx] = 0 // dummy atom
                                }

                                idx++

                            } else if (atom.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "Residues": //TODO
                        // Import additional (water) atom types
                        // Since "Deep" is true, it keeps all child nodes of the incoming Residue node
                        for (int i = 0; i < watDoc.getElementsByTagName("Residue").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Residue").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList residues = parentNode.item(0).getChildNodes()
                        for (Node res : residues) {
                            if (res.hasChildNodes()) {
                                String resName = res.getAttribute("name")

                                NodeList resProps = res.getChildNodes()

                                Element e = (Element) res;
                                String[] nameMap = new String[e.getElementsByTagName("Atom").length]
                                int[] typeMap = new int[e.getElementsByTagName("Atom").length]
                                int[] bondCount = new int[e.getElementsByTagName("Atom").length]
                                int i = 0
                                for (Node resProp : resProps) {
                                    if (resProp.getNodeName() == "Atom") {
                                        String atomName = resProp.getAttribute("name")
                                        int atomType = atomTypeMap.get(resProp.getAttribute("type")) - 1 // just used as index, so need this
                                        nameMap[i] = atomName
                                        typeMap[i] = atomType
                                        i++

                                    } else if (resProp.getNodeName() == "Bond") {
                                        int atom1 = -1
                                        int atom2 = -1
                                        if (resProp.hasAttribute("from")) {
                                            atom1 = parseInt(resProp.getAttribute("from"))
                                            atom2 = parseInt(resProp.getAttribute("to"))
                                        } else if (resProp.hasAttribute("atomName1")) {
                                            List<String> nameList = Arrays.asList(nameMap) // could just start with it as a List..
                                            atom1 = nameList.indexOf(resProp.getAttribute("atomName1"))
                                            atom2 = nameList.indexOf(resProp.getAttribute("atomName2"))
                                        }
                                        bondCount[atom1]++
                                        bondCount[atom2]++

                                    } else if (resProp.getNodeName() == "ExternalBond") {
                                        int extBondAtom = parseInt(resProp.getAttribute("from"))
                                        bondCount[extBondAtom]++ // could just use if statement above and get rid of this

                                    } else if (resProp.hasAttributes()) {
                                        logger.info("CHECK")
                                    }
                                }
                                for (int j = 0; j < typeMap.length; j++) {
                                    valence[typeMap[j]] = bondCount[j] //TODO could check to see if it equals the previous value in the valence array
                                    if (biotypeAtomNames[typeMap[j]] == "" || biotypeAtomNames[typeMap[j]] == null) {
                                        atomRes[typeMap[j]] = resName
                                        biotypeAtomNames[typeMap[j]] = nameMap[j]
                                    } else if (atomRes[typeMap[j]] != resName || biotypeAtomNames[typeMap[j]] != nameMap[j]) {
                                        atomRes[typeMap[j]] = atomRes[typeMap[j]] + ";" + resName
                                        biotypeAtomNames[typeMap[j]] = biotypeAtomNames[typeMap[j]] + ";" + nameMap[j]
                                    }
                                }
                            }
                        }
                        break

                    case "Patches": //TODO
                        logger.info("At Patches")

                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("Patch").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Patch").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        break

                    case "HarmonicBondForce":
                        // Import additional (water) atom types
                        // Created element to only grab "Bond" nodes within the "HarmonicBondForce" node (there are also Bond nodes in Residues)
                        Element ele = (Element) watDoc.getElementsByTagName("HarmonicBondForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Bond").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Bond").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        NodeList bonds = parentNode.item(0).getChildNodes()

                        for (Node bond : bonds) {
                            if (bond.getNodeName() == "Bond") {
                                String class1 = bond.getAttribute("class1")
                                String class2 = bond.getAttribute("class2")
                                String bondLength = bond.getAttribute("length")
                                String k = bond.getAttribute("k")
                                if (!atomClassMap.containsKey(class1)) {
                                    atomClassMap.put(class1, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassMap.size()+1)
                                }

                                int[] classes = new int[2]
                                classes[0] = atomClassMap.get(class1)
                                classes[1] = atomClassMap.get(class2)

                                double forceConstant = parseDouble(k) / KJperKCal / (ANGperNM * ANGperNM) / 2 // TODO why divided by 2?
                                double distance = parseDouble(bondLength) * ANGperNM

                                forceField.addForceFieldType(new BondType(classes, forceConstant, distance))
                            } else if (bond.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "HarmonicAngleForce":
                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("Angle").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("Angle").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        NodeList angles = parentNode.item(0).getChildNodes()

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
                                    atomClassMap.put(class1, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class2)) {
                                    atomClassMap.put(class2, atomClassMap.size()+1)
                                }
                                if (!atomClassMap.containsKey(class3)) {
                                    atomClassMap.put(class3, atomClassMap.size()+1)
                                }

                                int[] classes = new int[3]
                                classes[0] = atomClassMap.get(class1)
                                classes[1] = atomClassMap.get(class2)
                                classes[2] = atomClassMap.get(class3)
                                double[] angleArr = new double[1]
                                angleArr[0] = parseDouble(angleStr) * 180.0 / PI  // convert to degrees

                                double forceConstant = parseDouble(k) / KJperKCal / 2

                                forceField.addForceFieldType(new AngleType(classes, forceConstant, angleArr))
                            } else if (angle.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "AmoebaUreyBradleyForce": //TODO
                        logger.info("At AmoebaUreyBradleyForce")

                        // Import additional (water) atom types
                        for (int i = 0; i < watDoc.getElementsByTagName("UreyBradley").length; i++) {
                            Node watNode = doc.importNode(watDoc.getElementsByTagName("UreyBradley").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        break

                    case "PeriodicTorsionForce":
                        // Import additional (water) atom types
                        // Note: don't need this for water/ions xml, but here if we're combining other xml's
                        Element ele = (Element) watDoc.getElementsByTagName("PeriodicTorsionForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Proper").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Proper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        for (int i = 0; i < ele.getElementsByTagName("Imroper").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Imroper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList torsions = parentNode.item(0).getChildNodes()

                        for (Node torsion : torsions) {
                            if (torsion.getNodeName() == "Proper" || torsion.getNodeName() == "Improper") {

                                int numTerms = (torsion.getAttributes().length - 4) / 3  // get number of amplitudes/phases/periodicities
                                int[] classes = new int[4]
                                int[] periods = new int[numTerms]
                                double[] phases = new double[numTerms]
                                double[] amplitudes = new double[numTerms]
                                for (int i = 1; i <= 4; i++) {
                                    String torsClass= torsion.getAttribute("class"+toString(i))
                                    if (!atomClassMap.containsKey(torsClass)) {
                                        if (torsClass == "") {
                                            atomClassMap.put(torsClass, 0)
                                        } else {
                                            atomClassMap.put(torsClass, atomClassMap.size()+1)
                                        }
                                    }
                                    classes[i-1] = atomClassMap.get(torsClass)

                                    if (torsion.hasAttribute("periodicity"+toString(i))) { // could replace with i < numTerms
                                        periods[i-1] = parseInt(torsion.getAttribute("periodicity"+toString(i)))
                                        phases[i-1] = parseDouble(torsion.getAttribute("phase"+toString(i))) * 180.0 / PI  // convert to degrees
                                        amplitudes[i-1] = parseDouble(torsion.getAttribute("k"+toString(i))) / KJperKCal // convert to kcal/mol
                                    }
                                }

                                if (torsion.getNodeName() == "Proper") {
                                    forceField.addForceFieldType(new TorsionType(classes, amplitudes, phases, periods, TorsionType.TorsionMode.NORMAL))
                                } else {
                                    int trigAtom = classes[0] // Trigonal atom is defined first in XML, needs to be third in FFX/TINKER
                                    classes[0] = classes[1]
                                    classes[1] = classes[2]
                                    classes[2] = trigAtom
                                    forceField.addForceFieldType(new ImproperTorsionType(classes, amplitudes[0], phases[0], periods[0]))
                                }
                            }
                            else if (torsion.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "CustomTorsionForce":
                        logger.info("At CustomTorsionForce")

                        Element ele = (Element) watDoc.getElementsByTagName("CustomTorsionForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Proper").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Proper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }
                        for (int i = 0; i < ele.getElementsByTagName("Improper").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Improper").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList customTors = parentNode.item(0).getChildNodes()

                        for (Node torsion : customTors) {
                            if (torsion.getNodeName() == "Improper") {
                                int[] classes = new int[4]
                                for (int i = 1; i <= 4; i++) {
                                    String torsClass= torsion.getAttribute("class"+toString(i))
                                    if (!atomClassMap.containsKey(torsClass)) {
                                        if (torsClass == "") {
                                            atomClassMap.put(torsClass, 0)
                                        } else {
                                            atomClassMap.put(torsClass, atomClassMap.size()+1)
                                        }
                                    }
                                    classes[i-1] = atomClassMap.get(torsClass)
                                }
                                int period = 2 // always 2 for improper torsions, not defined in XML
                                double phase = parseDouble(torsion.getAttribute("theta0")) * 180.0 / PI // convert to degrees
                                double amplitude = parseDouble(torsion.getAttribute("k")) / KJperKCal // convert to kcal/mol

                                int trigAtom = classes[0] // Trigonal atom is defined first in XML, needs to be third in FFX/TINKER
                                classes[0] = classes[1]
                                classes[1] = classes[2]
                                classes[2] = trigAtom
                                forceField.addForceFieldType(new ImproperTorsionType(classes, amplitude, phase, period))
                            }
                            else if (torsion.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }

                        break

                    case "CMAPTorsionForce": //TODO
                        logger.info("At CMAPTorsionForce")
                        // Don't need to add waters/ions because they are do not have phi/psi angles
                        // <Map> : long lists doubles - look for new line chars?
                        // <Torsion> : map="0-11", class1-5

                        NodeList cmaps = parentNode.item(0).getChildNodes()
                        for (Node cmap : cmaps) {
//                            public TorsionTorsionType(int[] atomClasses, int[] gridPoints, double[] torsion1, double[] torsion2, double[] energy)
                            // int[5] atomClasses come from <Torsion>
                            // int[2] gridPoints - 25 25 in CHARMM_22_CMAP; should be length(torsion1) length(torsion2)
                            // double[] torsion1
                            // double[] torsion2
                            // double[] energy
                            if (cmap.getNodeName() == "Map") {
                                String nodeData = cmap.getChildNodes().item(0).getNodeValue() // gets data in <Map> node
                                //TODO
                            } else if (cmap.getNodeName() == "Torsion") {
                                int[] classes = new int[5]
                                String map = cmap.getAttribute("map")
                                classes[0] = cmap.getAttribute("class1")
                                classes[1] = cmap.getAttribute("class2")
                                classes[2] = cmap.getAttribute("class3")
                                classes[3] = cmap.getAttribute("class4")
                                classes[4] = cmap.getAttribute("class5")
                                //TODO
                            } else if (cmap.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "NonbondedForce":
                        // Import additional (water) atom types
                        // Created element to only grab "Atom" nodes within the "NonbondedForce" node (there are also Atom nodes in Residues)
                        Element ele = (Element) watDoc.getElementsByTagName("NonbondedForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Atom").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Atom").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList nbForces = parentNode.item(0).getChildNodes()

                        for (Node nbF : nbForces) {
                            if (nbF.getNodeName() == "Atom") {
                                int atomClass = atomClassMap.get(nbF.getAttribute("class"))
                                double sigma = parseDouble(nbF.getAttribute("sigma")) * ANGperNM / 2 * Math.pow(2.0,1.0/6.0) // nm to Ang and divide by 2 to get radius & r-min = 2^(1/6) * sigma
                                double eps = parseDouble(nbF.getAttribute("epsilon")) / KJperKCal  // kJ/mol to KCal/mol

                                forceField.addForceFieldType(new VDWType(atomClass, sigma, eps, -1.0))  // vdw by atom type
                            } else if (nbF.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break

                    case "LennardJonesForce":
                        logger.info("At LennardJonesForce")

                        // Import additional (water) atom types
                        // Created element to only grab "Atom" nodes within the "LennardJonesForce" node (there are also Atom nodes in other nodes)
                        Element ele = (Element) watDoc.getElementsByTagName("LennardJonesForce").item(0)
                        for (int i = 0; i < ele.getElementsByTagName("Atom").length; i++) {
                            Node watNode = doc.importNode(ele.getElementsByTagName("Atom").item(i), true)
                            parentNode.item(0).appendChild(watNode)
                        }

                        NodeList ljForces = parentNode.item(0).getChildNodes()

                        for (Node ljF : ljForces) {
                            if (ljF.getNodeName() == "Atom") {
                                int atomClass = atomClassMap.get(ljF.getAttribute("class"))
                                double sigma = parseDouble(ljF.getAttribute("sigma")) * ANGperNM / 2 * Math.pow(2.0,1.0/6.0) // nm to Ang and divide by 2 to get radius & r-min = 2^(1/6) * sigma
                                double eps = parseDouble(ljF.getAttribute("epsilon")) / KJperKCal  // kJ/mol to KCal/mol

                                if (ljF.hasAttribute("epsilon14")) {
                                    double eps14 = parseDouble(ljF.getAttribute("epsilon14")) / KJperKCal // kJ/mol to KCal/mol
                                    double sigma14 = parseDouble(ljF.getAttribue("sigma14")) * ANGperNM / 2 * Math.pow(2.0,1.0/6.0) // nm to Ang and divide by 2 to get radius & r-min = 2^(1/6) * sigma
                                    forceField.addForceFieldType(new VDWType(atomClass, sigma14, eps14, -1.0, VDWMode.VDW14))
                                }

                                // TODO - i think these should just replace the ones by NonbondedForce, they are incorrect -> if L.J.F. dont read N.B.F?
                                forceField.addForceFieldType(new VDWType(atomClass, sigma, eps, -1.0))  // vdw by atom type
                            } else if (ljF.getNodeName() == "NBFixPair") {
                                int[] classes = new int[2]
                                classes[0] = atomClassMap.get(ljF.getAttribute("class1"))
                                classes[1] = atomClassMap.get(ljF.getAttribute("class2"))
                                double sigma = parseDouble(ljF.getAttribute("sigma")) * ANGperNM / 2 * Math.pow(2.0,1.0/6.0) // nm to Ang and divide by 2 to get radius & r-min = 2^(1/6) * sigma
                                double eps = parseDouble(ljF.getAttribute("epsilon")) / KJperKCal // kJ/mol to KCal/mol

                                forceField.addForceFieldType(new VDWPairType(classes, sigma, eps))
                            } else if (ljF.hasAttributes()) {
                                logger.info("CHECK")
                            }
                        }
                        break
                }
            }
        }
        for (int j = 0; j < numAtomTypes; j++) {
            forceField.addForceFieldType(new AtomType(atomTypes[j], atomClasses[j], atomNames[j], "\"" + atomEnvs[j] + "\"", atomicNumbers[j], atomicWeights[j], valence[j]))
        }

        // Log atom type Map
//        Set<String> atStrSet = atomTypeMap.keySet()
//        for (String str : atStrSet) {
//            logger.info(format("%s - %d",str,atomTypeMap.get(str)))
//        }
        // Log atom class Map
//        Set<String> acStrSet = atomClassMap.keySet()
//        for (String str : acStrSet) {
//            logger.info(format("%s - %d",str,atomClassMap.get(str)))
//        }

        // Read in current AMBER_1999_SB in FFX to use its biotypes to set the XML FF's biotypes
        CompositeConfiguration props = Keyword.loadProperties(null)
        props.setProperty("parameters", "/iahome/j/jm/jmiller99/forcefieldx/modules/potential/src/main/java/ffx/potential/parameters/ff/AMBER_1999_SB")
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(props)
        ForceField predefFF = forceFieldFilter.parse()
        Map<String, BioType> biotypes = predefFF.getBioTypeMap()

        for (BioType biotype : biotypes.values()) {
            String resID = biotypeMap.get(biotype.moleculeName)

            String btAtomName = biotype.atomName
            int tempFlag = 0
            for (int idx = 0; idx < numAtomTypes; idx++) {
                String[] diffRes = atomRes[idx].split(";") //TODO should not do this for all atoms .. or shouldnt do it multiple times
                for (int j = 0; j < diffRes.length; j++) {
                    if (diffRes[j] == resID) {
                        String xmlAtomName = biotypeAtomNames[idx].split(";")[j] // split string at ';' and grab element corresponding to the correct residue
                        xmlAtomName = xmlAtomName.replace("'",'*') // biotype * = xml

                        if (btAtomName == xmlAtomName) {
                            // Add new BioType to forcefield with correct atom type corresponding to the xml
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            //add continue or break here?
                        } else if (btAtomName + "1" == xmlAtomName || btAtomName + "2" == xmlAtomName || btAtomName + "3" == xmlAtomName) {
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s",biotype.moleculeName,diffRes[j],btAtomName,xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                            // could compare atomtypes between the ones with numbers -> not best setup
                            // TODO: could add if statement to compare current atomType to new -> if they match fine, if not LOG
                        } else if (btAtomName == "HN") {
                            if (xmlAtomName == "H" || xmlAtomName == "H1" || xmlAtomName == "H2" || xmlAtomName == "H3") {
//                                if (xmlAtomName == "H" || xmlAtomName == "HN1" || xmlAtomName == "HN2") {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (biotype.moleculeName.contains(' Ion')) {
                            if (btAtomName + "+" == xmlAtomName || btAtomName + "-" == xmlAtomName) {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            } else if (btAtomName.length() == 2) {
                                String ionName = btAtomName[0] + btAtomName.toLowerCase()[1]
                                if (ionName + "+" == xmlAtomName || ionName + "-" == xmlAtomName) {
//                                    logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                    forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                    tempFlag = 1
                                }
                            }
                        } else if (btAtomName + "11" == xmlAtomName || btAtomName + "12" == xmlAtomName || btAtomName + "13" == xmlAtomName || btAtomName + "21" == xmlAtomName || btAtomName + "22" == xmlAtomName) {
//                            logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                            forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                            tempFlag = 1
                        } else if (btAtomName == "H") {
//                            if (btAtomName + "H" == xmlAtomName.substring(0, 2)) {
                            if (btAtomName + "H31" == xmlAtomName || btAtomName + "H32" == xmlAtomName || btAtomName + "H33" == xmlAtomName) {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
                                forceField.addForceFieldType(new BioType(biotype.index, biotype.atomName, biotype.moleculeName, atomTypes[idx]))
                                tempFlag = 1
                            }
                        } else if (btAtomName == "OP") {
                            if (xmlAtomName == "O1P" || xmlAtomName == "O2P") {
//                                logger.info(format("BT RES BTname XMLname: %s   %s      %s  %s", biotype.moleculeName, diffRes[j], btAtomName, xmlAtomName))
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

        // build header
        StringBuilder head = new StringBuilder()
//        head.append("forcefield\t\tAMBER-FF99SB-XML\n\n" +
//                    "vdwtype\t\t\tLENNARD-JONES\n" +
//                    "vdwindex\t\tTYPE\n" +
//                    "radiusrule\t\tARITHMETIC\n" +
//                    "radiustype\t\tR-MIN\n" + // or SIGMA
//                    "radiussize\t\tRADIUS\n" + // or DIAMETER
//                    "epsilonrule\t\tGEOMETRIC\n" +
//                    "vdw-14-scale\t\t0.500000\n" + // or 2.0
//                    "chg-14-scale\t\t0.833333\n" + // or 1.2
//                    "electric\t\t332.0522173\n" +
//                    "dielectric\t\t1.0\n")
        head.append("vdwtype LENNARD-JONES\n" +
                "radiusrule ARITHMETIC\n" +
                "radiustype R-MIN\n" + // or SIGMA
                "radiussize RADIUS\n" + // or DIAMETER
                "epsilonrule GEOMETRIC\n" +
                "vdw-14-scale 0.500000\n" + // or 2.0
                "chg-14-scale 0.833333\n" + // or 1.2
                "electric 332.0522173\n" +
                "dielectric 1.0\n")
// TODO: can get vdw-14-scale and chg-14-scale from the xml Node NonbondedForce attributes "lj14scale" and "coulomb14scale", respectively
        StringBuffer ffSB = forceField.toStringBuffer() // convert FF to string buffer
//        String ffStr = ffSB.replaceAll('improper','imptors') // replace improper FF terms with 'imptors'

        // write forcefield file
        BufferedWriter out = new BufferedWriter(new FileWriter("sampleFF.txt"))
        out.write(head.toString())
        out.write(ffSB.toString())
        out.flush()
        out.close()
    }
}