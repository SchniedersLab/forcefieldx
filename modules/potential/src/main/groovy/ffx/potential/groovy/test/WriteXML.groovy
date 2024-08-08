package ffx.potential.groovy.test

import ffx.potential.cli.PotentialScript

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

import ffx.potential.parameters.*
import ffx.potential.parsers.ForceFieldFilter
import ffx.utilities.Keyword
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.filefilter.FalseFileFilter
import org.biojava.nbio.core.util.FileDownloadUtils
import org.w3c.dom.Document
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.w3c.dom.NodeList
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.management.InvalidAttributeValueException
import javax.xml.crypto.dsig.dom.DOMSignContext
import javax.xml.parsers.DocumentBuilder
import javax.xml.parsers.DocumentBuilderFactory
import javax.xml.transform.OutputKeys
import javax.xml.transform.Transformer
import javax.xml.transform.TransformerFactory
import javax.xml.transform.dom.DOMSource
import javax.xml.transform.stream.StreamResult
import java.text.SimpleDateFormat

import static ffx.potential.bonded.Atom.ElementSymbol
import static java.lang.Double.parseDouble
import static java.lang.Integer.parseInt
import static java.lang.Integer.toString
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.PI

/**
 * The ReadXML script converts a prm file to an XML file usable by OpenMM.
 *
 * @author Jacob M. Miller
 * <br>
 * Usage:
 * <br>
 * ffxc test.ReadXML &lt;filename&gt;
 */
@Command(description = "Write an OpenMM XML file.", name = "ffxc test.WriteXML")
class WriteXML extends PotentialScript {

    /**
     * -c or --charmm enables reading the XML in as a Charmm forcefield.
     */
    @Option(names = ['-c', '--charmm'], paramLabel = "", defaultValue = "false",
    description = "Parse forcefield file as Charmm instead of Amber")
    private boolean charmmRead = false;

    /**
     * The final argument(s) should be three filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'XML FF file to be read.')
    List<String> arguments = null

    final double ANGperNM = 10.0
    final double KJperKCal = 4.184

//    LinkedHashMap<String, Integer> atomTypeMap
//    LinkedHashMap<String, Integer> atomClassMap
//    LinkedHashMap<String, String> biotypeMap

    /**
     * Execute the script.
     */
    @Override
    WriteXML run() {

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

//        atomTypeMap = new LinkedHashMap<>()
//        atomClassMap = new LinkedHashMap<>()
//        biotypeMap = new LinkedHashMap<>()

        // specify files coming in
//        File biotypeClasses = new File(arguments[0])  // biotype description to xml Residue Map TODO - could we use residuesFinal.xml
//        File inputFile = new File(arguments[1]) // main xml FF file
//        File watFile = new File(arguments[2]) // water xml FF file (e.g. TIP3P)
        File prmFile = new File(arguments[0])
        File resFile = new File("/Users/jakemiller/ffx_testing/writeXML/residuesFinal.xml")

        // CREATE XML ROOT NODES AND FILL IN BASIC INFO
        // Instantiate Document building objects
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance()
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder()
        Document doc = dBuilder.newDocument()
        Document resDoc = dBuilder.parse(resFile)

        // Create Root Node "ForceField"
        Element rootElement = doc.createElement("ForceField")
        doc.appendChild(rootElement)

        // Create Info Section
        Element infoNode = doc.createElement("Info")
        rootElement.appendChild(infoNode)
        // If wanted attribute in this node: infoNode.setAttribute("id", "1001") e.g.
        // an Attribute is inside the node definition - <Node Att="n">...</Node>

        Element srcNode = doc.createElement("Source")
        srcNode.setTextContent(prmFile.name)
        infoNode.appendChild(srcNode)

        Element dateNode = doc.createElement("DateGenerated")
        SimpleDateFormat ft = new SimpleDateFormat("yyyy-MM-dd")
        String date = ft.format(new Date())
        dateNode.setTextContent(date)
        infoNode.appendChild(dateNode)

        Element refNode = doc.createElement("Reference")
        infoNode.appendChild(refNode)

        // GO THROUGH FORCE FIELD TYPES AND ADD TO XML
        // Creat PRM file's Force Field
        CompositeConfiguration props = Keyword.loadProperties(null)
        props.setProperty("parameters", prmFile.path)
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(props)
        ForceField ff = forceFieldFilter.parse()

        // Add Atom Types
        Element atomTypesNode = doc.createElement("AtomTypes")
        rootElement.appendChild(atomTypesNode)
        Map<String, AtomType> atMap = ff.getTypes(ForceField.ForceFieldType.ATOM)
        for (AtomType at : atMap.values()) {
            Element atomType = doc.createElement("Type")
            at.toXML(atomType)
            atomTypesNode.appendChild(atomType)
        }

        // Add Residues??       HashMap<String, AtomType> x = ff.getAtomTypes("GLY") TODO try for residues
        buildResidueDict(resDoc)

        // Add Bond Types
        Element bondTypesNode = doc.createElement("AmoebaBondedForce")
        rootElement.appendChild(bondTypesNode)
        bondTypesNode.setAttribute("bond-cubic", format("%f",ff.getDouble("bond-cubic")*10.0))
        bondTypesNode.setAttribute("bond-quartic", format("%f",ff.getDouble("bond-quartic")*100.0))
        Map<String, BondType> btMap = ff.getTypes(ForceField.ForceFieldType.BOND)
        for (BondType bt : btMap.values()) {
            Element bondType = doc.createElement("Bond")
            bt.toXML(bondType)
            bondTypesNode.appendChild(bondType)
        }

        // Add Angle Types
        Element angleTypesNode = doc.createElement("AmoebaAngleForce")
        rootElement.appendChild(angleTypesNode)
        angleTypesNode.setAttribute("angle-cubic", ff.getDouble("angle-cubic").toString())
        angleTypesNode.setAttribute("angle-quartic", ff.getDouble("angle-quartic").toString())
        angleTypesNode.setAttribute("angle-pentic", ff.getDouble("angle-pentic").toString())
        angleTypesNode.setAttribute("angle-sextic", ff.getDouble("angle-sextic").toString())
        Map<String, AngleType> angMap = ff.getTypes(ForceField.ForceFieldType.ANGLE)
        for (AngleType ang : angMap.values()) {
            Element angleType = doc.createElement("Angle")
            ang.toXML(angleType)
            angleTypesNode.appendChild(angleType)
        }

        // Add OutOfPlane Bend Types
        Element oopbTypesNode = doc.createElement("AmoebaOutOfPlaneBendForce")
        rootElement.appendChild(oopbTypesNode)
        oopbTypesNode.setAttribute("type", ff.getString("opbendtype"))
        oopbTypesNode.setAttribute("opbend-cubic", ff.getDouble("opbend-cubic").toString())
        oopbTypesNode.setAttribute("opbend-quartic", ff.getDouble("opbend-quartic").toString())
        oopbTypesNode.setAttribute("opbend-pentic", ff.getDouble("opbend-pentic").toString())
        oopbTypesNode.setAttribute("opbend-sextic", ff.getDouble("opbend-sextic").toString())
        Map<String, OutOfPlaneBendType> oopbMap = ff.getTypes(ForceField.ForceFieldType.OPBEND)
        for (OutOfPlaneBendType oopb : oopbMap.values()) {
            Element oopbType = doc.createElement("Angle")
            oopb.toXML(oopbType)
            oopbTypesNode.appendChild(oopbType)
        }

        // Add Torsion Types
        Element torsTypesNode = doc.createElement("PeriodicTorsionForce")
        rootElement.appendChild(torsTypesNode)
        Map<String, TorsionType> torsMap = ff.getTypes(ForceField.ForceFieldType.TORSION)
        double torsionUnit = ff.getDouble("torsionunit")
        for (TorsionType tors : torsMap.values()) {
            Element torsType = doc.createElement("Proper")
            tors.toXML(torsType,torsionUnit)
            torsTypesNode.appendChild(torsType)
        }

        // Add PiOrbitalTorsion Types
        Element piTorsTypesNode = doc.createElement("AmoebaPiTorsionForce") // todo should we not create node if not PiTorsTypes in FF?
        rootElement.appendChild(piTorsTypesNode)
        piTorsTypesNode.setAttribute("piTorsionUnit", ff.getDouble("pitorsunit", 1.0).toString()) // todo should have default for all instances
        Map<String, PiOrbitalTorsionType> ptMap = ff.getTypes(ForceField.ForceFieldType.PITORS)
        for (PiOrbitalTorsionType pt : ptMap.values()) {
            Element piTorsType = doc.createElement("PiTorsion")
            pt.toXML(piTorsType)
            piTorsTypesNode.appendChild(piTorsType)
        }

        // Add StretchTorsionType's
        Element strTorsTypesNode = doc.createElement("AmoebaStretchTorsionForce")
        rootElement.appendChild(strTorsTypesNode)
        Map<String, StretchTorsionType> stMap = ff.getTypes(ForceField.ForceFieldType.STRTORS)
        for (StretchTorsionType st : stMap.values()) {
            Element strTorsType = doc.createElement("Torsion")
            st.toXML(strTorsType)
            strTorsTypesNode.appendChild(strTorsType)
        }

        // Add AngleTorsionType's
        Element angTorsTypesNode = doc.createElement("AmoebaAngleTorsionForce")
        rootElement.appendChild(angTorsTypesNode)
        Map<String, AngleTorsionType> angTorsMap = ff.getTypes(ForceField.ForceFieldType.ANGTORS)
        for (AngleTorsionType angTors : angTorsMap.values()) {
            Element angTorsType = doc.createElement("Torsion")
            angTors.toXML(angTorsType)
            angTorsTypesNode.appendChild(angTorsType)
        }

        // Add StretchBendType's
        Element strBendTypesNode = doc.createElement("AmoebaStretchBendForce")
        rootElement.appendChild(strBendTypesNode)
        strBendTypesNode.setAttribute("stretchBendUnit", format("%f",ff.getDouble("strbndunit", PI/180)*180/PI)) // TODO OpenMM has value of 1.0 (so multipled by 180/pi, but unnecessary?)
        Map<String, StretchBendType> sbMap = ff.getTypes(ForceField.ForceFieldType.STRBND)
        for (StretchBendType sb : sbMap.values()) {
            Element strBendType = doc.createElement("StretchBend")
            sb.toXML(strBendType)
            strBendTypesNode.appendChild(strBendType)
        }

        // Add TorsionTorsionType's
        Element tortorsNode = doc.createElement("AmoebaTorsionTorsionForce")
        rootElement.appendChild(tortorsNode)
        Map<String, TorsionTorsionType> ttMap = ff.getTypes(ForceField.ForceFieldType.TORTORS)
        int i = 0
        for (TorsionTorsionType tt : ttMap.values()) {
            Element tortors = doc.createElement("TorsionTorsion")
            Element ttGrid = doc.createElement("TorsionTorsionGrid")
            tt.toXML(doc, tortors, ttGrid) // send doc to method to create additional child nodes under ttGrid
            tortors.setAttribute("grid", toString(i))
            ttGrid.setAttribute("grid", toString(i))
            tortorsNode.appendChild(tortors)
            tortorsNode.appendChild(ttGrid)
            i++
        }

        // Add VDWType
        Element vdwTypesNode = doc.createElement("AmoebaVdwForce")
        rootElement.appendChild(vdwTypesNode)
        vdwTypesNode.setAttribute("type", ff.getString("vdwtype"))
        vdwTypesNode.setAttribute("radiusrule", ff.getString("radiusrule"))
        vdwTypesNode.setAttribute("radiustype", ff.getString("radiustype"))
        vdwTypesNode.setAttribute("radiussize", ff.getString("radiussize"))
        vdwTypesNode.setAttribute("epsilonrule", ff.getString("epsilonrule"))
        vdwTypesNode.setAttribute("vdw-13-scale", format("%f",ff.getDouble("vdw-13-scale")))
        vdwTypesNode.setAttribute("vdw-14-scale", format("%f",ff.getDouble("vdw-14-scale")))
        vdwTypesNode.setAttribute("vdw-15-scale", format("%f",ff.getDouble("vdw-15-scale")))
        Map<String, VDWType> vdwMap = ff.getTypes(ForceField.ForceFieldType.VDW)
        for (VDWType vdw : vdwMap.values()) {
            Element vdwType = doc.createElement("Vdw")
            vdw.toXML(vdwType)
            vdwTypesNode.appendChild(vdwType)
        }
        Map<String, VDWPairType> vdwPMap = ff.getTypes(ForceField.ForceFieldType.VDWPR)
        for (VDWPairType vdwP : vdwPMap.values()) {
            Element vdwPairType = doc.createElement("Pair")
            vdwP.toXML(vdwPairType)
            vdwTypesNode.appendChild(vdwPairType)
        }

        // Add MultipoleType and PolarizeType
        Element mpTypesNode = doc.createElement("AmoebaMultipoleForce")
        rootElement.appendChild(mpTypesNode)
        mpTypesNode.setAttribute("direct11Scale", format("%f",ff.getDouble("direct-11-scale")))
        mpTypesNode.setAttribute("direct12Scale", format("%f",ff.getDouble("direct-12-scale")))
        mpTypesNode.setAttribute("direct13Scale", format("%f",ff.getDouble("direct-13-scale")))
        mpTypesNode.setAttribute("direct14Scale", format("%f",ff.getDouble("direct-14-scale")))
        mpTypesNode.setAttribute("mpole12Scale", format("%f",ff.getDouble("mpole-12-scale")))
        mpTypesNode.setAttribute("mpole13Scale", format("%f",ff.getDouble("mpole-13-scale")))
        mpTypesNode.setAttribute("mpole14Scale", format("%f",ff.getDouble("mpole-14-scale")))
        mpTypesNode.setAttribute("mpole15Scale", format("%f",ff.getDouble("mpole-15-scale")))
        mpTypesNode.setAttribute("mutual11Scale", format("%f",ff.getDouble("mutual-11-scale")))
        mpTypesNode.setAttribute("mutual12Scale", format("%f",ff.getDouble("mutual-12-scale")))
        mpTypesNode.setAttribute("mutual13Scale", format("%f",ff.getDouble("mutual-13-scale")))
        mpTypesNode.setAttribute("mutual14Scale", format("%f",ff.getDouble("mutual-14-scale")))
        mpTypesNode.setAttribute("polar12Scale", format("%f",ff.getDouble("polar-12-scale")))
        mpTypesNode.setAttribute("polar13Scale", format("%f",ff.getDouble("polar-13-scale")))
        mpTypesNode.setAttribute("polar14Intra", format("%f",ff.getDouble("polar-14-intra")))
        mpTypesNode.setAttribute("polar14Scale", format("%f",ff.getDouble("polar-14-scale")))
        mpTypesNode.setAttribute("polar15Scale", format("%f",ff.getDouble("polar-15-scale"))) // todo why dont they have all (e.g. intra)
        Map<String, MultipoleType> mpMap = ff.getTypes(ForceField.ForceFieldType.MULTIPOLE)
        for (MultipoleType mp : mpMap.values()) {
            Element mpType = doc.createElement("Multipole")
            mp.toXML(mpType)
            mpTypesNode.appendChild(mpType)
        }
        Map<String, PolarizeType> polMap = ff.getTypes(ForceField.ForceFieldType.POLARIZE)
        for (PolarizeType pol : polMap.values()) {
            Element polType = doc.createElement("Polarize")
            pol.toXML(polType)
            mpTypesNode.appendChild(polType)
        }

        // Add UreyBradleyType
        Element ubTypesNode = doc.createElement("AmoebaUreyBradleyForce")
        rootElement.appendChild(ubTypesNode)
        ubTypesNode.setAttribute("cubic", format("%f",ff.getDouble("urey-cubic", 0.0)))
        ubTypesNode.setAttribute("quartic", format("%f",ff.getDouble("urey-quartic", 0.0)))
        Map<String, UreyBradleyType> ubMap = ff.getTypes(ForceField.ForceFieldType.UREYBRAD)
        for (UreyBradleyType ub : ubMap.values()) {
            Element ubType = doc.createElement("UreyBradley")
            ub.toXML(ubType)
            ubTypesNode.appendChild(ubType)
        }


        // Write XML to baseName.xml
        writeXML(doc, prmFile.baseName)

        return this
    }

    // Build entry for protein residue
    private static void buildProteinResidue(Map residueDict, Map atoms, Map bondInfo, String abbr, String loc,
                                            String tinkerLookupName, boolean include, String residueName, String type) {
        // Add an inner hash map with all of this residue's info
        residueDict.put(abbr, new HashMap<>())
        residueDict.get(abbr).put("atoms", atoms)
        residueDict.get(abbr).put("type", type)
        residueDict.get(abbr).put("loc", loc)
        residueDict.get(abbr).put("tinkerLookupName", tinkerLookupName)
        residueDict.get(abbr).put("residueName", residueName)
        residueDict.get(abbr).put("include", include)

        // # for each bond, add entry to
        //    #   residueDict[abbr]['atoms'][atom]['bonds']
        //    #   residueDict[abbr]['atoms'][bondedAtom]['bonds']
        for (String atom : bondInfo.keySet()) {
            if (residueDict.get(abbr).get("atoms").containsKey(atom)) {
                if (!residueDict.get(abbr).get("atoms").get(atom).containsKey("bonds")) {
                    residueDict.get(abbr).get("atoms").get(atom).put("bonds", new HashMap<>())
                }

                for (String bondedAtom : bondInfo.get(atom)) { // todo as is or keySet?

                }
            }
            logger.info(atom.getAttribute("from"))
        }
    }

    private static void buildResidueDict(Document resDoc) {
//        Map<String, Map<String, Object>> residueDict = new HashMap<>()
        Map residueDict = new HashMap()
        NodeList residues = resDoc.getElementsByTagName("Residue")
        for (Node residue : residues) {
            String abbr = residue.getAttribute("abbreviation")
            String loc = residue.getAttribute("loc")
            String type = residue.getAttribute("type")
            String tinkerName = residue.getAttribute("tinkerLookupName")
            String residueName = residue.getAttribute("fullName")

            Element e = (Element) residue
            Map atoms = getXmlAtoms(e.getElementsByTagName("Atom"))
            Map bondInfo = getXmlBonds(e.getElementsByTagName("Bond"))

            if (type == "AmoebaWater") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, "x", tinkerName, true, "HOH", "water")
            } else if (type == "protein") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, "m", tinkerName, true, residueName, "protein")
            } else if (type == "dna" || type == "rna") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, loc, tinkerName, true, residueName, type)
            }


        }
    }

    // Default 'constructor' for atoms
    private static Map getDefaultAtom() {
        Map atom = new HashMap()
        atom.put("tinkerLookupName", "XXX")
        atom.put("type", "-1")
        atom.put("bonds", new HashMap())
        return atom
    }

    // Get atom Map from xml atom list
    private static Map getXmlAtoms(NodeList atoms) {
        Map atomInfo = new HashMap()
        for (Node atom : atoms) {
            String name = atom.getAttribute("name")
            atomInfo.put(name, getDefaultAtom())
            atomInfo.get(name).get("tinkerLookupName").put(atom.getAttribute("tinkerLookupName"))
        }
        return atomInfo
    }

    // Get bond Map from xml bond list
    private static Map getXmlBonds(NodeList bonds) {
        Map bondInfo = new HashMap()
        for (Node bond : bonds) {
            String atom1 = bond.getAttribute("from")
            String atom2 = bond.getAttribute("to")
            if (!bondInfo.containsKey(atom1)) {
                bondInfo.put(atom1, new HashMap())
            }
            if (!bondInfo.containsKey(atom2)) {
                bondInfo.put(atom2, new HashMap())
            }

            bondInfo.get(atom1).put(atom2, "1")
            bondInfo.get(atom2).put(atom1, "1")
        }
        return bondInfo
    }

    private static void writeXML(Document doc, String outputName) {
        TransformerFactory tfFactory = TransformerFactory.newInstance()
        Transformer transformer = tfFactory.newTransformer()

        transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes") // if want to hide first line XML version/encoding
        transformer.setOutputProperty(OutputKeys.INDENT, "yes") // print with indentation

        DOMSource src = new DOMSource(doc)
        StreamResult result = new StreamResult(outputName + ".xml")
        transformer.transform(src,result)
    }
}