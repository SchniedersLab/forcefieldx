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
import org.w3c.dom.Document
import org.w3c.dom.Element
import org.w3c.dom.Node
import org.w3c.dom.NodeList
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import javax.xml.parsers.DocumentBuilder
import javax.xml.parsers.DocumentBuilderFactory
import javax.xml.transform.OutputKeys
import javax.xml.transform.Transformer
import javax.xml.transform.TransformerFactory
import javax.xml.transform.dom.DOMSource
import javax.xml.transform.stream.StreamResult
import java.text.SimpleDateFormat

import static java.lang.Integer.toString
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.PI

/**
 * The PrmToXml script converts a prm file to an XML file usable by OpenMM.
 *
 * @author Jacob M. Miller
 * <br>
 * Usage:
 * <br>
 * ffxc test.ReadXML &lt;filename&gt;
 */
@Command(description = "Write an OpenMM XML file.", name = "ffxc test.WriteXML")
class PrmToXml extends PotentialScript {

    /**
     * The final argument(s) should be three filenames.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'XML FF file to be read.')
    List<String> arguments = null

    /**
     * Execute the script.
     */
    @Override
    PrmToXml run() {

        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // specify files
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
        refNode.setTextContent("Yue Shi, Zhen Xia, Jiajing Zhang, Robert Best, Chuanjie Wu, Jay W. Ponder, and Pengyu Ren. Polarizable Atomic Multipole-Based AMOEBA Force Field for Proteins. Journal of Chemical Theory and Computation, 9(9):4046–4063, 2013.")
        Element refNode2 = doc.createElement("Reference")
        refNode2.setTextContent("Changsheng Zhang, Chao Lu, Zhifeng Jing, Chuanjie Wu, Jean-Philip Piquemal, Jay W Ponder, and Pengyu Ren. AMOEBA Polarizable Atomic Multipole Force Field for Nucleic Acids. Journal of Chemical Theory and Computation, 14(4):2084-2108, 2018.")
        infoNode.appendChild(refNode)
        infoNode.appendChild(refNode2)

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

        // Add Residues
        Map bioTypeDict = new HashMap()
        Map<String, BioType> bioTypes = ff.getBioTypeMap()
        Map moleculeDict = new HashMap()
        for (BioType bioType : bioTypes.values()) {
            String lookUp = bioType.atomName + "_" + bioType.moleculeName
            if (bioTypeDict.containsKey(lookUp)) {
                // Workaround for Tinker using the same name but different types for H2', H2'', and for H5', H5''
                lookUp = bioType.atomName + "*_" + bioType.moleculeName
            }
            bioTypeDict.put(lookUp, new HashMap())
            bioTypeDict.get(lookUp).put("index", bioType.index) // 0
            bioTypeDict.get(lookUp).put("atomName", bioType.atomName) // 1
            bioTypeDict.get(lookUp).put("moleculeName", bioType.moleculeName) // 2
            bioTypeDict.get(lookUp).put("atomType", bioType.atomType) // 3

            if (!moleculeDict.containsKey(bioType.moleculeName)) {
                moleculeDict.put(bioType.moleculeName, new HashMap())
                moleculeDict.get(bioType.moleculeName).put("added", false) // todo can change this to true in original methods? can add to make getExtraBiotypes simpler
                moleculeDict.get(bioType.moleculeName).put("atoms", new HashMap())
            }
            moleculeDict.get(bioType.moleculeName).get("atoms").put(bioType.atomName, new HashMap())
            moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("index", bioType.index)
            moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("atomType", bioType.atomType)
            moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("bonds", bioType.bonds)
//            moleculeDict.get(bioType.moleculeName).put(bioType.atomType, bioType.atomName)
        }
        Map residueDict = buildResidueDict(resDoc, doc)
//        List<String> extraBTs = getExtraBiotypes(bioTypeDict, residueDict) // todo DELETE
        setBioTypes(bioTypeDict, residueDict) // set biotypes for all atoms

        Element residuesNode = createResidueNodes(doc, residueDict, bioTypeDict)
        buildExtraResidues(moleculeDict, doc, residuesNode)
        rootElement.appendChild(residuesNode)


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

    // Add bond to atomMap; atoms are added to atomMap if missing
    private static void addBond(Map atomDict, String atom1, String atom2) {
        if (!atomDict.containsKey(atom1)) {
            atomDict.put(atom1, getDefaultAtom())
        }

        if (!atomDict.containsKey(atom2)) {
            atomDict.put(atom2, getDefaultAtom())
        }

        atomDict.get(atom2).get("bonds").put(atom1, 1)
        atomDict.get(atom1).get("bonds").put(atom2, 1)
    }

    // Build residues not found in OpenMM finalResidues.xml
    private static void buildExtraResidues(Map moleculeDict, Document doc, Element residuesNode) {
        for (String mol : moleculeDict.keySet()) {
            Element resNode = doc.createElement("Residue")
            resNode.setAttribute("name", mol)
            residuesNode.appendChild(resNode)
            for (String atom : moleculeDict.get(mol).get("atoms").keySet()) {
                Element atomNode = doc.createElement("Atom")
                atomNode.setAttribute("name", atom)
                atomNode.setAttribute("type", toString(moleculeDict.get(mol).get("atoms").get(atom).get("atomType")))
                resNode.appendChild(atomNode)

                String[] bonds = moleculeDict.get(mol).get("atoms").get(atom).get("bonds")
                for (String bond : bonds) { // todo just store here and have another loop or..
                    Element bondNode = doc.createElement("Bond")
                    bondNode.setAttribute("from", atom)
                    bondNode.setAttribute("to", bond)
                    resNode.appendChild(bondNode)
                }
            }
        }

        // TODO - won't be able to determine an external bond? - won't be any for guest molecule
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

                for (String bondedAtom : bondInfo.get(atom).keySet()) {
                    if (residueDict.get(abbr).get("atoms").containsKey(bondedAtom)) {
                        if (!residueDict.get(abbr).get("atoms").get(bondedAtom).containsKey("bonds")) {
                            residueDict.get(abbr).get("atoms").get(bondedAtom).put("bonds", new HashMap<>())
                        }
                        residueDict.get(abbr).get("atoms").get(bondedAtom).get("bonds").put(atom, 1)
                        residueDict.get(abbr).get("atoms").get(atom).get("bonds").put(bondedAtom, 1)
                    } else {
                        logger.info(format("ERROR: bonded atom=%s not in residue=%s",atom,abbr))
                    }
                }
            } else {
                logger.info(format("ERROR: bonded atom=%s nt in residue=%s",atom,abbr))
            }
        }
    }

    // Build residue hash based on xml file
    private static Map buildResidueDict(Document resDoc, Document xmlDoc) {
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

            // if residue is an amino acid, then create C-ter and N-ter residues, in addition to non-terminal residue,
            // and include appropriate atoms (HXT is excluded from all residues)
            if (type == "AmoebaWater") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, "x", tinkerName, true, "HOH", "water")

            } else if (type == "protein") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, "m", tinkerName, true, residueName, "protein")

                String cResidueName = "C" + abbr
                residueDict.put(cResidueName, copyProteinResidue(residueDict.get(abbr)))
                residueDict.get(cResidueName).put("loc", "c")
                if (residueDict.get(abbr).get("tinkerLookupName").indexOf("(") > -1) {
                    int begin = residueDict.get(abbr).get("tinkerLookupName").indexOf("(")
                    int end = residueDict.get(abbr).get("tinkerLookupName").indexOf(")") + 1
                    String sub = residueDict.get(abbr).get("tinkerLookupName").substring(begin,end)
                    if (sub == "(HD)" || sub == "(HE)") {
                        residueDict.get(cResidueName).replace("tinkerLookupName", "C-Terminal " + "HIS " + sub)
                    } else {
                        residueDict.get(cResidueName).replace("tinkerLookupName", "C-Terminal " + abbr + " " + sub)
                    }
                    logger.info(format("tinkerLookupName %s %s",abbr,residueDict.get(cResidueName).get("tinkerLookupName")))
                } else {
                    residueDict.get(cResidueName).replace("tinkerLookupName", "C-Terminal " + abbr)
                }
                residueDict.get(cResidueName).get("atoms").put("OXT", copyAtom(residueDict.get(abbr).get("atoms").get("O")))
                residueDict.get(cResidueName).get("atoms").get("OXT").replace("tinkerLookupName", "OXT")
                residueDict.get(cResidueName).get("atoms").get("O").replace("tinkerLookupName", "OXT")
                residueDict.get(cResidueName).put("parent", residueDict.get(abbr))

                String nResidueName = "N" + abbr
                residueDict.put(nResidueName, copyProteinResidue(residueDict.get(abbr)))
                residueDict.get(nResidueName).put("loc", "n")
                residueDict.get(nResidueName).put("tinkerLookupName", "N-Terminal " + abbr)
                residueDict.get(nResidueName).put("parent", residueDict.get(abbr))

                if (abbr == "PRO") {
                    // <Atom name="H" tinkerLookupName="HN" bonds="1" />
                    residueDict.get(nResidueName).get("atoms").put("H2", getDefaultAtom())
                    residueDict.get(nResidueName).get("atoms").put("H3", getDefaultAtom())

                    residueDict.get(nResidueName).get("atoms").get("H2").put("tinkerLookupName", "HN")
                    residueDict.get(nResidueName).get("atoms").get("H3").put("tinkerLookupName", "HN")

                    addBond(residueDict.get(nResidueName).get("atoms"), "H2", "N")
                    addBond(residueDict.get(nResidueName).get("atoms"), "H2", "N")

                } else {
                    residueDict.get(nResidueName).get("atoms").put("H2", copyAtom(residueDict.get(abbr).get("atoms").get("H")))
                    residueDict.get(nResidueName).get("atoms").put("H3", copyAtom(residueDict.get(abbr).get("atoms").get("H")))
                }

            } else if (type == "dna" || type == "rna") {
                buildProteinResidue(residueDict, atoms, bondInfo, abbr, loc, tinkerName, true, residueName, type)
            }
        }

//        logger.info("Start Lookup XML FFFFinal\n\n")
        boolean printResMapping = false
        if (printResMapping) {
            Element residuesNode = xmlDoc.createElement("Residues")
            for (String resName : residueDict.keySet().sort()) {
                if (residueDict.get(resName).containsKey("include") && residueDict.get(resName).get("include")) {
                    Element resNode = xmlDoc.createElement("Residue")
                    residuesNode.appendChild(resNode)

                    String type = residueDict.get(resName).get("type")
                    String loc = residueDict.get(resName).get("loc")
                    String tinkerLookupName = residueDict.get(resName).get("tinkerLookupName")
                    String fullName = residueDict.get(resName).get("residueName")
                    resNode.setAttribute("abbreviation", resName)
                    resNode.setAttribute("loc", loc)
                    resNode.setAttribute("type", type)
                    resNode.setAttribute("tinkerLookupName", tinkerLookupName)
                    resNode.setAttribute("fullName", fullName)

                    Map atomsInfo = residueDict.get(resName).get("atoms")
                    Map includedBonds = new HashMap()
                    for (String atomName : atomsInfo.keySet().sort()) {
                        Element atomNode = xmlDoc.createElement("Atom")
                        resNode.appendChild(atomNode)

                        tinkerLookupName = atomsInfo.get(atomName).get("tinkerLookupName")
                        atomNode.setAttribute("name", atomName)
                        atomNode.setAttribute("tinkerLookupName", tinkerLookupName)

                        Map bondsInfo = atomsInfo.get(atomName).get("bonds")
                        for (String bondedAtom : bondsInfo.keySet()) {
                            if (!includedBonds.containsKey(bondedAtom) || !includedBonds.get(bondedAtom).containsKey(atomName)) {
                                Element bondNode = xmlDoc.createElement("Bond")
                                resNode.appendChild(bondNode)

                                bondNode.setAttribute("from", atomName)
                                bondNode.setAttribute("to", bondedAtom)

                                if (!includedBonds.containsKey(atomName)) {
                                    includedBonds.put(atomName, new HashMap())
                                }
                                if (!includedBonds.containsKey(bondedAtom)) {
                                    includedBonds.put(bondedAtom, new HashMap())
                                }
                                includedBonds.get(atomName).put(bondedAtom, 1)
                                includedBonds.get(bondedAtom).put(atomName, 1)
                            }
                        }
                    }
                }
            }
        }
        return residueDict
    }

    // Copy a atom (Map copy, including the 'bonds')
    private static Map copyAtom(Map atom) {
        Map atomCopy = new HashMap()
        for (String key : atom.keySet()) {
            if (key != "bonds") {
                atomCopy.put(key, atom.get(key))
            } else {
                atomCopy.put("bonds", copyBonds(atom.get(key)))
            }
        }
        return atomCopy
    }

    // Copy a bond (Map copy)
    private static Map copyBonds(Map bonds) {
        Map bondCopy = new HashMap()
        for (String key : bonds.keySet()) {
            bondCopy.put(key, bonds.get(key))
        }
        return bondCopy
    }

    // Copy a residue, including atom list
    private static Map copyProteinResidue(Map residue) {
        Map residueCopy = new HashMap()
        residueCopy.put("atoms", new HashMap())
        residueCopy.put("type", residue.get("type"))
        residueCopy.put("loc", residue.get("loc"))
        residueCopy.put("tinkerLookupName", residue.get("tinkerLookupName"))
        residueCopy.put("residueName", residue.get("residueName"))
        residueCopy.put("include", residue.get("include"))

        for (String atom : residue.get("atoms").keySet()) {
            residueCopy.get("atoms").put(atom, copyAtom(residue.get("atoms").get(atom)))
        }
        return residueCopy
    }

    // Create an "End cap" residue (additional residues)
    private static Element createEndCap(Document doc, String[] resInfo, String[] names, String[] btNames, String[] froms, String[] tos, String extBondInt, Map bioTypes) {
        Element resNode = doc.createElement("Residue")
        resNode.setAttribute("name", resInfo[0])
        for (int i=0; i < names.length; i++) {
            Element atomNode = doc.createElement("Atom")
            atomNode.setAttribute("name", names[i])
            atomNode.setAttribute("type", toString(bioTypes.get(btNames[i]+"_"+resInfo[1]).get("atomType")))
            resNode.appendChild(atomNode)

            if (i < froms.length) {
                Element bondNode = doc.createElement("Bond")
                bondNode.setAttribute("from", froms[i])
                bondNode.setAttribute("to", tos[i])
                resNode.appendChild(bondNode)
            } else {
                Element extBond = doc.createElement("ExternalBond")
                extBond.setAttribute("from", extBondInt)
                resNode.appendChild(extBond)
            }
        }
        return resNode
    }

    // Create the Residues Node and add Residue sub-nodes
    private static Element createResidueNodes(Document doc, Map residueDict, Map bioTypes) {
        Element residuesNode = doc.createElement("Residues")
        for (String resname : residueDict.keySet().sort()) {
            Map res = residueDict.get(resname)
            logger.info(format("INCLUDE: %s", res.get("include")))
            if (res.get("include")) {
                Element resNode = doc.createElement("Residue")
                residuesNode.appendChild(resNode)
                resNode.setAttribute("name", resname)
                if (resname == "HOH") {
                    resNode.setAttribute("rigidWater", "false")
                }
                Map atomIndex = new HashMap()
                int atomCount = 0
                for (String atom : res.get("atoms").keySet().sort()) {
                    int type = res.get("atoms").get(atom).get("type")
                    if (type < 0) {
                        logger.info(format("ERROR: type=%d for atom=%s of residue=%s",type,atom,resname))
                    }
                    Element atomNode = doc.createElement("Atom")
                    atomNode.setAttribute("name", atom)
                    atomNode.setAttribute("type", toString(type))
                    resNode.appendChild(atomNode)
                    atomIndex.put(atom, atomCount)
                    atomCount++
                }

                Map includedBonds = new HashMap()
                for (String atomName : res.get("atoms").keySet().sort()) {
                    Map bondsInfo = res.get("atoms").get(atomName).get("bonds")
                    for (String bondedAtom : bondsInfo.keySet()) {
                        if (!includedBonds.containsKey(bondedAtom) || !includedBonds.get(bondedAtom).containsKey(atomName)) {
                            Element bondNode = doc.createElement("Bond")
                            bondNode.setAttribute("from", toString(atomIndex.get(atomName)))
                            bondNode.setAttribute("to", toString(atomIndex.get(bondedAtom)))
                            resNode.appendChild(bondNode)
                            if (!includedBonds.containsKey(atomName)) {
                                includedBonds.put(atomName, new HashMap())
                            }
                            if (!includedBonds.containsKey(bondedAtom)) {
                                includedBonds.put(bondedAtom, new HashMap())
                            }
                            includedBonds.get(atomName).put(bondedAtom, 1)
                            includedBonds.get(bondedAtom).put(atomName, 1)
                        }
                    }
                }

                if (res.get("type") == "rna" || res.get("type") == "dna") {
                    if (res.get("loc") == "middle" || res.get("loc") == "3") {
                        Element extBond = doc.createElement("ExternalBond")
                        extBond.setAttribute("from", toString(atomIndex.get("P")))
                        resNode.appendChild(extBond)
                    }
                    if (res.get("loc") == "middle" || res.get("loc") == "5") {
                        Element extBond = doc.createElement("ExternalBond")
                        extBond.setAttribute("from", toString(atomIndex.get("O3'")))
                        resNode.appendChild(extBond)
                    }
                } else {
                    if (res.get("loc") == "m") {
                        Element extBond1 = doc.createElement("ExternalBond")
                        Element extBond2 = doc.createElement("ExternalBond")
                        extBond1.setAttribute("from", toString(atomIndex.get("N")))
                        extBond2.setAttribute("from", toString(atomIndex.get("C")))
                        resNode.appendChild(extBond1)
                        resNode.appendChild(extBond2)
                    }
                    if (res.get("loc") == "n") {
                        Element extBond = doc.createElement("ExternalBond")
                        extBond.setAttribute("from", toString(atomIndex.get("C")))
                        resNode.appendChild(extBond)
                    }
                    if (res.get("loc") == "c") {
                        Element extBond = doc.createElement("ExternalBond")
                        extBond.setAttribute("from", toString(atomIndex.get("N")))
                        resNode.appendChild(extBond)
                    }
                    if (resname.indexOf("CYX") > -1) {
                        Element extBond = doc.createElement("ExternalBond")
                        extBond.setAttribute("from", toString(atomIndex.get("SG")))
                        resNode.appendChild(extBond)
                    }
                }
            }
        }

        // End caps
        // Residue == ACE
        if (bioTypes.containsKey("H_Acetyl N-Terminus")) {
            String[] resInfo = ["ACE", "Acetyl N-Terminus"] // resname and biotype molecule name
            String[] names = ["HH31", "CH3", "HH32", "HH33", "C", "O"]
            String[] btNames = ["H", "CH3", "H", "H", "C", "O"]
            String[] froms = ["0", "1", "1", "1", "4"]
            String[] tos = ["1", "2", "3", "4", "5"]
            String extBondInt = "4"
            Element resNode = createEndCap(doc, resInfo, names, btNames, froms, tos, extBondInt, bioTypes)
            residuesNode.appendChild(resNode)
        }

        // Residue == NME
        if (bioTypes.containsKey("N_N-MeAmide C-Terminus")) {
            String[] resInfo = ["NME", "N-MeAmide C-Terminus"] // resname and biotype molecule name
            String[] names = ["N", "H", "CH3", "HH31", "HH32", "HH33"]
            String[] btNames = ["N", "HN", "CH3", "H", "H", "H"]
            String[] froms = ["0", "0", "2", "2", "2"]
            String[] tos = ["1", "2", "3", "4", "5"]
            String extBondInt = "0"
            Element resNode = createEndCap(doc, resInfo, names, btNames, froms, tos, extBondInt, bioTypes)
            residuesNode.appendChild(resNode)
        }

        // Ions TODO finish

        return residuesNode
    }

    // Default 'constructor' for atoms
    private static Map getDefaultAtom() {
        Map atom = new HashMap()
        atom.put("tinkerLookupName", "XXX")
        atom.put("type", -1)
        atom.put("bonds", new HashMap())
        return atom
    }

    // Find any biotypes that aren't in residueDict
    private static List<String> getExtraBiotypes(Map bioTypeDict, Map residueDict) {
        List<String> extraBTs = new ArrayList<>()
        for (String lookupName : bioTypeDict.keySet()) {
            String[] splits = lookupName.split("_")
            boolean done = false
            String molName = splits[1]
            if (splits.length > 2) {
                for (int i = 2; i < splits.length; i++) {
                    molName = molName + "_" + splits[i]
                }
            }
            for (String resname : residueDict.keySet()) {
                if (residueDict.get(resname).get("tinkerLookupName") == molName) {
//                        logger.info(format("TEST: %s", lookupName))
                    done = true
                    break
                }
            }
            if (!done) {
//                logger.info(format("EXTRA: %s", lookupName))
                extraBTs.add(lookupName)
            }
        }
        logger.info(format("LASTTEST: %d", extraBTs.size()))
        return extraBTs
    }

    // Get atom Map from xml atom list
    private static Map getXmlAtoms(NodeList atoms) {
        Map atomInfo = new HashMap()
        for (Node atom : atoms) {
            String name = atom.getAttribute("name")
            atomInfo.put(name, getDefaultAtom())
//            atomInfo.get(name).get("tinkerLookupName").put(atom.getAttribute("tinkerLookupName"))
            atomInfo.get(name).replace("tinkerLookupName", atom.getAttribute("tinkerLookupName"))
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

            bondInfo.get(atom1).put(atom2, 1)
            bondInfo.get(atom2).put(atom1, 1)
        }
        return bondInfo
    }

    // Set biotype for each atom in residueDict
    private static void setBioTypes(Map bioTypes, Map residueDict) {
        for (String resname : residueDict.keySet()) {
            Map res = residueDict.get(resname)
            for (String atom : res.get("atoms").keySet()) {
                String atomLookup = res.get("atoms").get(atom).get("tinkerLookupName")
                List<String> resLookup = new ArrayList<>()
                if (res.get("type") == "dna") {
                    if (res.get("loc") == "5" || res.get("loc") == "N") {
                        resLookup.add("5'-Hydroxyl DNA")
                    }
                    if (res.get("loc") == "3" || res.get("loc") == "N") {
                        resLookup.add("3'-Hydroxyl DNA")
                    }
                    resLookup.add("Phosphodiester DNA")
                }
                if (res.get("type") == "rna") {
                    if (res.get("loc") == "5" || res.get("loc") == "N") {
                        resLookup.add("5'-Hydroxyl RNA")
                    }
                    if (res.get("loc") == "3" || res.get("loc") == "N") {
                        resLookup.add("3'-Hydroxyl RNA")
                    }
                    resLookup.add("Phosphodiester RNA")
                }
                resLookup.add(res.get("tinkerLookupName"))
                String lookupName = ""
                for (String suffix : resLookup) {
                    lookupName = atomLookup + "_" + suffix
                    if (bioTypes.containsKey(lookupName)) {
                        break
                    }
                }
                if (bioTypes.containsKey(lookupName)) {
                    res.get("atoms").get(atom).put("type", bioTypes.get(lookupName).get("atomType"))
//                    usedAtomType.put(bioTypes.get(lookupName).get("atomType")) TODO THIS IS AN OPTION
                } else {
                    logger.info(format("For %s lookupName=%s not in biotype",atom,lookupName))
                    if (res.containsKey("parent")) {
                        lookupName = res.get("atoms").get(atom).get("tinkerLookupName") + "_" + res.get("parent").get("tinkerLookupName")
                        if (bioTypes.containsKey(lookupName)) {
                            res.get("atoms").get(atom).put("type", bioTypes.get(lookupName).get("atomType"))
                        } else {
                            logger.info(format("MISSING lookupName=%s from biotype",lookupName))
                            res.replace("include", false)
                        }
                    } else {
                        res.replace("include", false)
                    }
                }
            }
        }
    }
    //TODO - THOUGHT: could have a tracker for all atoms seeing if they get placed into a biotype, if not, must create new residue

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