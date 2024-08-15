// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
// ******************************************************************************
package ffx.potential.parsers;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;

import ffx.potential.parameters.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Logger;

/**
 * The OpenMMXmlFilter class creates an OpenMM-style XML file from a ForceField object. Work is underway to parse
 * OpenMM XML files into ForceField objects (AMBER & CHARMM).
 *
 * @author Jacob M. Miller
 */
public class OpenMMXmlFilter {

    /** A Logger for the XmlFilter class. */
    private static final Logger logger = Logger.getLogger(OpenMMXmlFilter.class.getName());

    /** The ForceField object that will be written to an XML. */
    private final ForceField forceField;

    /**
     * Constructor for outputting XML.
     *
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     */
    public OpenMMXmlFilter(ForceField forceField) {
        this.forceField = forceField;
    }

    /**
     * Create an OpenMM XML file for the given force field.
     *
     * @throws Exception
     */
    public void toXML() throws Exception {
        // CREATE XML ROOT NODES AND FILL IN BASIC INFO
        // Instantiate Document building objects
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.newDocument();

        // Create Root Node "ForceField"
        Element rootElement = doc.createElement("ForceField");
        doc.appendChild(rootElement);

        // Create Info Section
        Element infoNode = doc.createElement("Info");
        rootElement.appendChild(infoNode);

        // Set Source node to the force field name
        Element srcNode = doc.createElement("Source");
        srcNode.setTextContent(forceField.getString("forcefield", "UNKNOWN"));
        infoNode.appendChild(srcNode);

        // Set the date
        Element dateNode = doc.createElement("DateGenerated");
        SimpleDateFormat ft = new SimpleDateFormat("yyyy-MM-dd");
        String date = ft.format(new Date());
        dateNode.setTextContent(date);
        infoNode.appendChild(dateNode);

        // Adds references for the AMOEBA force field.
        Element refNode = doc.createElement("Reference");
        refNode.setTextContent("Yue Shi, Zhen Xia, Jiajing Zhang, Robert Best, Chuanjie Wu, Jay W. Ponder, and Pengyu Ren. Polarizable Atomic Multipole-Based AMOEBA Force Field for Proteins. Journal of Chemical Theory and Computation, 9(9):4046–4063, 2013.");
        Element refNode2 = doc.createElement("Reference");
        refNode2.setTextContent("Changsheng Zhang, Chao Lu, Zhifeng Jing, Chuanjie Wu, Jean-Philip Piquemal, Jay W Ponder, and Pengyu Ren. AMOEBA Polarizable Atomic Multipole Force Field for Nucleic Acids. Journal of Chemical Theory and Computation, 14(4):2084-2108, 2018.");
        infoNode.appendChild(refNode);
        infoNode.appendChild(refNode2);

        // Add Atom Types
        Element atomTypesNode = doc.createElement("AtomTypes");
        rootElement.appendChild(atomTypesNode);
        Map<String, AtomType> atMap = (Map<String, AtomType>) forceField.getTypes(ForceField.ForceFieldType.ATOM);
        for (AtomType at : atMap.values()) {
            Element atomType = doc.createElement("Type");
            at.toXML(atomType);
            atomTypesNode.appendChild(atomType);
        }


        // Add Residues
        Map<String, BioType> bioTypes = forceField.getBioTypeMap();
        if (bioTypes.values().size() >= 1) {
            Map<String, Map<String, Map<String, Map<String, Object>>>> moleculeDict = new HashMap<>();
            for (BioType bioType : bioTypes.values()) {
                if (bioType.bonds != null) {
                    if (!moleculeDict.containsKey(bioType.moleculeName)) {
                        moleculeDict.put(bioType.moleculeName, new HashMap<>());
//                    moleculeDict.get(bioType.moleculeName).put("added", false); // todo can change this to true in original methods? can add to make getExtraBiotypes simpler
                        moleculeDict.get(bioType.moleculeName).put("atoms", new HashMap<>());
                    }
                    moleculeDict.get(bioType.moleculeName).get("atoms").put(bioType.atomName, new HashMap<>());
                    moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("index", bioType.index);
                    moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("atomType", bioType.atomType);
                    moleculeDict.get(bioType.moleculeName).get("atoms").get(bioType.atomName).put("bonds", bioType.bonds);
                }
            }
            Element residuesNode = doc.createElement("Residues");
            buildExtraResidues(moleculeDict, doc, residuesNode);
            rootElement.appendChild(residuesNode);
        } else {
            logger.info("WARNING: No BioTypes defined for force field, therefore no OpenMM residues are created.");
        }

        // Add Bond Types
        Map<String, BondType> btMap = (Map<String, BondType>) forceField.getTypes(ForceField.ForceFieldType.BOND);
        if (btMap.values().size() >= 1) {
            Element bondTypesNode = doc.createElement("AmoebaBondedForce");
            rootElement.appendChild(bondTypesNode);
            bondTypesNode.setAttribute("bond-cubic", format("%f",forceField.getDouble("bond-cubic", 0.0)*10.0));
            bondTypesNode.setAttribute("bond-quartic", format("%f",forceField.getDouble("bond-quartic", 0.0)*100.0));
            for (BondType bt : btMap.values()) {
                Element bondType = doc.createElement("Bond");
                bt.toXML(bondType);
                bondTypesNode.appendChild(bondType);
            }
        }

        // Add Angle Types
        Map<String, AngleType> angMap = (Map<String, AngleType>) forceField.getTypes(ForceField.ForceFieldType.ANGLE);
        if (angMap.values().size() >= 1) {
            Element angleTypesNode = doc.createElement("AmoebaAngleForce");
            rootElement.appendChild(angleTypesNode);
            angleTypesNode.setAttribute("angle-cubic", String.valueOf(forceField.getDouble("angle-cubic", 0.0)));
            angleTypesNode.setAttribute("angle-quartic", String.valueOf(forceField.getDouble("angle-quartic", 0.0)));
            angleTypesNode.setAttribute("angle-pentic", String.valueOf(forceField.getDouble("angle-pentic", 0.0)));
            angleTypesNode.setAttribute("angle-sextic", String.valueOf(forceField.getDouble("angle-sextic", 0.0)));
            for (AngleType ang : angMap.values()) {
                Element angleType = doc.createElement("Angle");
                ang.toXML(angleType);
                angleTypesNode.appendChild(angleType);
            }
        }

        // Add OutOfPlane Bend Types
        Map<String, OutOfPlaneBendType> oopbMap = (Map<String, OutOfPlaneBendType>) forceField.getTypes(ForceField.ForceFieldType.OPBEND);
        if (oopbMap.values().size() >= 1) {
            Element oopbTypesNode = doc.createElement("AmoebaOutOfPlaneBendForce");
            rootElement.appendChild(oopbTypesNode);
            oopbTypesNode.setAttribute("type", forceField.getString("opbendtype", "ALLINGER"));
            oopbTypesNode.setAttribute("opbend-cubic", String.valueOf(forceField.getDouble("opbend-cubic", 0.0)));
            oopbTypesNode.setAttribute("opbend-quartic", String.valueOf(forceField.getDouble("opbend-quartic", 0.0)));
            oopbTypesNode.setAttribute("opbend-pentic", String.valueOf(forceField.getDouble("opbend-pentic", 0.0)));
            oopbTypesNode.setAttribute("opbend-sextic", String.valueOf(forceField.getDouble("opbend-sextic", 0.0)));
            for (OutOfPlaneBendType oopb : oopbMap.values()) {
                Element oopbType = doc.createElement("Angle");
                oopb.toXML(oopbType);
                oopbTypesNode.appendChild(oopbType);
            }
        }

        // Add Torsion Types
        Map<String, TorsionType> torsMap = (Map<String, TorsionType>) forceField.getTypes(ForceField.ForceFieldType.TORSION);
        if (torsMap.values().size() >= 1) {
            Element torsTypesNode = doc.createElement("PeriodicTorsionForce");
            rootElement.appendChild(torsTypesNode);
            double torsionUnit = forceField.getDouble("torsionunit", 1.0);
            for (TorsionType tors : torsMap.values()) {
                Element torsType = doc.createElement("Proper");
                tors.toXML(torsType,torsionUnit);
                torsTypesNode.appendChild(torsType);
            }
        }

        // Add PiOrbitalTorsion Types
        Map<String, PiOrbitalTorsionType> ptMap = (Map<String, PiOrbitalTorsionType>) forceField.getTypes(ForceField.ForceFieldType.PITORS);
        if (ptMap.values().size() >= 1) {
            Element piTorsTypesNode = doc.createElement("AmoebaPiTorsionForce");
            rootElement.appendChild(piTorsTypesNode);
            piTorsTypesNode.setAttribute("piTorsionUnit", String.valueOf(forceField.getDouble("pitorsunit", 1.0)));
            for (PiOrbitalTorsionType pt : ptMap.values()) {
                Element piTorsType = doc.createElement("PiTorsion");
                pt.toXML(piTorsType);
                piTorsTypesNode.appendChild(piTorsType);
            }
        }

        // Add StretchTorsionType's
        Map<String, StretchTorsionType> stMap = (Map<String, StretchTorsionType>) forceField.getTypes(ForceField.ForceFieldType.STRTORS);
        if (stMap.values().size() >= 1) {
            Element strTorsTypesNode = doc.createElement("AmoebaStretchTorsionForce");
            rootElement.appendChild(strTorsTypesNode);
            for (StretchTorsionType st : stMap.values()) {
                Element strTorsType = doc.createElement("Torsion");
                st.toXML(strTorsType);
                strTorsTypesNode.appendChild(strTorsType);
            }
        }

        // Add AngleTorsionType's
        Map<String, AngleTorsionType> angTorsMap = (Map<String, AngleTorsionType>) forceField.getTypes(ForceField.ForceFieldType.ANGTORS);
        if (angTorsMap.values().size() >= 1) {
            Element angTorsTypesNode = doc.createElement("AmoebaAngleTorsionForce");
            rootElement.appendChild(angTorsTypesNode);
            for (AngleTorsionType angTors : angTorsMap.values()) {
                Element angTorsType = doc.createElement("Torsion");
                angTors.toXML(angTorsType);
                angTorsTypesNode.appendChild(angTorsType);
            }
        }

        // Add StretchBendType's
        Map<String, StretchBendType> sbMap = (Map<String, StretchBendType>) forceField.getTypes(ForceField.ForceFieldType.STRBND);
        if (sbMap.values().size() >= 1) {
            Element strBendTypesNode = doc.createElement("AmoebaStretchBendForce");
            rootElement.appendChild(strBendTypesNode);
            strBendTypesNode.setAttribute("stretchBendUnit", String.valueOf(forceField.getDouble("strbndunit", PI/180)*180/PI));
            for (StretchBendType sb : sbMap.values()) {
                Element strBendType = doc.createElement("StretchBend");
                sb.toXML(strBendType);
                strBendTypesNode.appendChild(strBendType);
            }
        }

        // Add TorsionTorsionType's
        Map<String, TorsionTorsionType> ttMap = (Map<String, TorsionTorsionType>) forceField.getTypes(ForceField.ForceFieldType.TORTORS);
        if (ttMap.values().size() >= 1) {
            Element tortorsNode = doc.createElement("AmoebaTorsionTorsionForce");
            rootElement.appendChild(tortorsNode);
            int i = 0;
            for (TorsionTorsionType tt : ttMap.values()) {
                Element tortors = doc.createElement("TorsionTorsion");
                Element ttGrid = doc.createElement("TorsionTorsionGrid");
                tt.toXML(doc, tortors, ttGrid); // send doc to method to create additional child nodes under ttGrid
                tortors.setAttribute("grid", String.valueOf(i));
                ttGrid.setAttribute("grid", String.valueOf(i));
                tortorsNode.appendChild(tortors);
                tortorsNode.appendChild(ttGrid);
                i++;
            }
        }

        // Add VDWType
        Map<String, VDWType> vdwMap = (Map<String, VDWType>) forceField.getTypes(ForceField.ForceFieldType.VDW);
        Map<String, VDWPairType> vdwPMap = (Map<String, VDWPairType>) forceField.getTypes(ForceField.ForceFieldType.VDWPR);
        if (vdwMap.values().size() >= 1 || vdwPMap.values().size() >= 1) {
            Element vdwTypesNode = doc.createElement("AmoebaVdwForce");
            rootElement.appendChild(vdwTypesNode);
            vdwTypesNode.setAttribute("type", forceField.getString("vdwtype", "LENNARD-JONES"));
            vdwTypesNode.setAttribute("radiusrule", forceField.getString("radiusrule", "ARITHMETIC"));
            vdwTypesNode.setAttribute("radiustype", forceField.getString("radiustype", "R-MIN"));
            vdwTypesNode.setAttribute("radiussize", forceField.getString("radiussize", "RADIUS"));
            vdwTypesNode.setAttribute("epsilonrule", forceField.getString("epsilonrule", "GEOMETRIC"));
            // todo vdw-12-scale?
            vdwTypesNode.setAttribute("vdw-13-scale", String.valueOf(forceField.getDouble("vdw-13-scale", 0.0)));
            vdwTypesNode.setAttribute("vdw-14-scale", String.valueOf(forceField.getDouble("vdw-14-scale", 1.0)));
            vdwTypesNode.setAttribute("vdw-15-scale", String.valueOf(forceField.getDouble("vdw-15-scale", 1.0)));
            for (VDWType vdw : vdwMap.values()) {
                Element vdwType = doc.createElement("Vdw");
                vdw.toXML(vdwType);
                vdwTypesNode.appendChild(vdwType);
            }
            for (VDWPairType vdwP : vdwPMap.values()) {
                Element vdwPairType = doc.createElement("Pair");
                vdwP.toXML(vdwPairType);
                vdwTypesNode.appendChild(vdwPairType);
            }
        }

        // Add MultipoleType and PolarizeType
        Map<String, MultipoleType> mpMap = (Map<String, MultipoleType>) forceField.getTypes(ForceField.ForceFieldType.MULTIPOLE);
        Map<String, PolarizeType> polMap = (Map<String, PolarizeType>) forceField.getTypes(ForceField.ForceFieldType.POLARIZE);
        if (mpMap.values().size() >= 1 || polMap.values().size() >= 1) {
            Element mpTypesNode = doc.createElement("AmoebaMultipoleForce");
            rootElement.appendChild(mpTypesNode);
            mpTypesNode.setAttribute("direct11Scale", String.valueOf(forceField.getDouble("direct-11-scale", 0.0)));
            mpTypesNode.setAttribute("direct12Scale", String.valueOf(forceField.getDouble("direct-12-scale", 1.0)));
            mpTypesNode.setAttribute("direct13Scale", String.valueOf(forceField.getDouble("direct-13-scale", 1.0)));
            mpTypesNode.setAttribute("direct14Scale", String.valueOf(forceField.getDouble("direct-14-scale", 1.0)));
            mpTypesNode.setAttribute("mpole12Scale", String.valueOf(forceField.getDouble("mpole-12-scale", 0.0)));
            mpTypesNode.setAttribute("mpole13Scale", String.valueOf(forceField.getDouble("mpole-13-scale", 0.0)));
            mpTypesNode.setAttribute("mpole14Scale", String.valueOf(forceField.getDouble("mpole-14-scale", 1.0)));
            mpTypesNode.setAttribute("mpole15Scale", String.valueOf(forceField.getDouble("mpole-15-scale", 1.0)));
            mpTypesNode.setAttribute("mutual11Scale", String.valueOf(forceField.getDouble("mutual-11-scale", 1.0)));
            mpTypesNode.setAttribute("mutual12Scale", String.valueOf(forceField.getDouble("mutual-12-scale", 1.0)));
            mpTypesNode.setAttribute("mutual13Scale", String.valueOf(forceField.getDouble("mutual-13-scale", 1.0)));
            mpTypesNode.setAttribute("mutual14Scale", String.valueOf(forceField.getDouble("mutual-14-scale", 1.0)));
            mpTypesNode.setAttribute("polar12Scale", String.valueOf(forceField.getDouble("polar-12-scale", 0.0)));
            mpTypesNode.setAttribute("polar13Scale", String.valueOf(forceField.getDouble("polar-13-scale", 0.0)));
            mpTypesNode.setAttribute("polar14Intra", String.valueOf(forceField.getDouble("polar-14-intra", 0.5)));
            mpTypesNode.setAttribute("polar14Scale", String.valueOf(forceField.getDouble("polar-14-scale", 1.0)));
            mpTypesNode.setAttribute("polar15Scale", String.valueOf(forceField.getDouble("polar-15-scale", 1.0)));
            // todo why dont they have all (e.g. intra) [processTinkerForceField.py]
            for (MultipoleType mp : mpMap.values()) {
                Element mpType = doc.createElement("Multipole");
                mp.toXML(mpType);
                mpTypesNode.appendChild(mpType);
            }
            for (PolarizeType pol : polMap.values()) {
                Element polType = doc.createElement("Polarize");
                pol.toXML(polType);
                mpTypesNode.appendChild(polType);
            }
        }

        // Add UreyBradleyType
        Map<String, UreyBradleyType> ubMap = (Map<String, UreyBradleyType>) forceField.getTypes(ForceField.ForceFieldType.UREYBRAD);
        if (ubMap.values().size() >= 1) {
            Element ubTypesNode = doc.createElement("AmoebaUreyBradleyForce");
            rootElement.appendChild(ubTypesNode);
            ubTypesNode.setAttribute("cubic", String.valueOf(forceField.getDouble("urey-cubic", 0.0)));
            ubTypesNode.setAttribute("quartic", String.valueOf(forceField.getDouble("urey-quartic", 0.0)));
            for (UreyBradleyType ub : ubMap.values()) {
                Element ubType = doc.createElement("UreyBradley");
                ub.toXML(ubType);
                ubTypesNode.appendChild(ubType);
            }
        }

        // Write XML to 'force field name'.xml
        writeXML(doc, forceField.getString("forcefield", "UNKNOWN"));
    }

    /**
     * Build residues not found in OpenMM finalResidues.xml todo - should this just replace that?
     * @param moleculeDict Map containing molecules defined in BioTypes along with individual atom information.
     * @param doc Document with XML Nodes.
     * @param residuesNode Residues Node that will contain individual Residue nodes with Atoms and Bonds.
     */
    private static void buildExtraResidues(Map<String, Map<String, Map<String, Map<String, Object>>>> moleculeDict, Document doc, Element residuesNode) {
        // Go through all molecules defined in the moleculeDict (keys)
        for (String mol : moleculeDict.keySet()) {
            Element resNode = doc.createElement("Residue"); // create Residue node for each molecule
            resNode.setAttribute("name", mol);
            residuesNode.appendChild(resNode);

            List<String> order = new ArrayList<>();
            List<String> froms = new ArrayList<>();
            List<String> tos = new ArrayList<>();
            // Go through all atoms inside the molecule (keys)
            for (String atom : moleculeDict.get(mol).get("atoms").keySet()) {
                Element atomNode = doc.createElement("Atom");
                atomNode.setAttribute("name", atom);
                atomNode.setAttribute("type", String.valueOf(moleculeDict.get(mol).get("atoms").get(atom).get("atomType")));
                resNode.appendChild(atomNode);

                String[] bonds = (String[]) moleculeDict.get(mol).get("atoms").get(atom).get("bonds");
                // Go through all atoms (bonds) that the defined atom is bonded to
                for (String bond : bonds) {
                    if (!moleculeDict.get(mol).get("atoms").get(atom).containsKey("usedBonds")) {
                        moleculeDict.get(mol).get("atoms").get(atom).put("usedBonds", new ArrayList<>());
                    }
                    if (!moleculeDict.get(mol).get("atoms").get(bond).containsKey("usedBonds")) {
                        moleculeDict.get(mol).get("atoms").get(bond).put("usedBonds", new ArrayList<>());
                    }

                    // get the bonds already done for each atom
                    List<String> atomBonds = (ArrayList<String>) moleculeDict.get(mol).get("atoms").get(atom).get("usedBonds");
                    List<String> bondBonds = (ArrayList<String>) moleculeDict.get(mol).get("atoms").get(bond).get("usedBonds");

                    // record the bond in from/to lists if the bond has not already been used
                    if (!atomBonds.contains(bond) && !bondBonds.contains(atom)) {
                        froms.add(atom);
                        tos.add(bond);
                    }

                    // store that the bond has been identified
                    atomBonds.add(bond);
                    bondBonds.add(atom);
                    moleculeDict.get(mol).get("atoms").get(atom).put("usedBonds", atomBonds);
                    moleculeDict.get(mol).get("atoms").get(bond).put("usedBonds", bondBonds);
                }
                order.add(atom); // add atom in the order they have been added to the Residue node
            }

            // add Bond nodes using the order the Atom nodes were placed
            for (int i = 0; i < froms.size(); i++) {
                Element bondNode = doc.createElement("Bond");
                bondNode.setAttribute("from", String.valueOf(order.indexOf(froms.get(i))));
                bondNode.setAttribute("to", String.valueOf(order.indexOf(tos.get(i))));
                resNode.appendChild(bondNode);
            }
        }

        // TODO - won't be able to determine an external bond? - won't be any for guest molecule
        // this would be only problem if didn't use OMM defined Residue records
    }

    /**
     * Create an OpenMM-style XML file from the Document object that was created in toXML().
     * @param doc Document object containing XML nodes.
     * @param outputName String with name of the new file.
     * @throws TransformerException
     */
    private static void writeXML(Document doc, String outputName) throws TransformerException {
        TransformerFactory tfFactory = TransformerFactory.newInstance();
        Transformer transformer = tfFactory.newTransformer();

        transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes"); // if want to hide first line XML version/encoding
        transformer.setOutputProperty(OutputKeys.INDENT, "yes"); // print with indentation

        DOMSource src = new DOMSource(doc);
        StreamResult result = new StreamResult(outputName + ".xml");
        transformer.transform(src,result);
    }
}
