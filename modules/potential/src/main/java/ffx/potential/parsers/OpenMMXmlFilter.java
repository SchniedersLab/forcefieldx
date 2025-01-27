// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.parsers;

import ffx.potential.parameters.AngleTorsionType;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BioType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiOrbitalTorsionType;
import ffx.potential.parameters.StretchBendType;
import ffx.potential.parameters.StretchTorsionType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * The OpenMMXmlFilter class creates an OpenMM-style XML file from a ForceField object. Work is underway to parse
 * OpenMM XML files into ForceField objects (AMBER & CHARMM).
 *
 * @author Jacob M. Miller
 */
public class OpenMMXmlFilter {

  /**
   * A Logger for the XmlFilter class.
   */
  private static final Logger logger = Logger.getLogger(OpenMMXmlFilter.class.getName());

  /**
   * The ForceField object that will be written to an XML.
   */
  private final ForceField forceField;

  /**
   * The path to write the XML file to. Should not include '.xml' on end.
   */
  private String outputName;

  /**
   * Constructor for outputting XML.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public OpenMMXmlFilter(ForceField forceField) {
    this.forceField = forceField;
  }

  /**
   * Constructor for outputting XML with output path specified.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param saveName   a String with the output path.
   */
  public OpenMMXmlFilter(ForceField forceField, String saveName) {
    this.forceField = forceField;
    this.outputName = saveName;
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

    // Add Atom Types
    Element atomTypes = AtomType.getXMLAtomTypes(doc, forceField);
    if (atomTypes != null) {
      rootElement.appendChild(atomTypes);
    } else {
      throw new Exception(" The force field defines no atom types.");
    }

    // Add Residues
    Map<String, BioType> bioTypes = forceField.getBioTypeMap();
    if (!bioTypes.values().isEmpty()) {
      Map<String, ArrayList<BioType>> moleculeDict = new HashMap<>();
      for (BioType bioType : bioTypes.values()) {
        if (!moleculeDict.containsKey(bioType.moleculeName)) {
          moleculeDict.put(bioType.moleculeName, new ArrayList<>());
        }
        ArrayList<BioType> moleculeBioTypes = moleculeDict.get(bioType.moleculeName);
        moleculeBioTypes.add(bioType);
      }
      Element residuesNode = doc.createElement("Residues");
      buildExtraResidues(moleculeDict, doc, residuesNode);
      rootElement.appendChild(residuesNode);
    } else {
      logger.info("WARNING: No BioTypes defined for force field, therefore no OpenMM residues are created.");
    }

    // Add Bond Types
    Element bondForce = BondType.getXMLForce(doc, forceField);
    if (bondForce != null) {
      rootElement.appendChild(bondForce);
    }

    // Add Angle Types
    Element angleForce = AngleType.getXMLForce(doc, forceField);
    if (angleForce != null) {
      rootElement.appendChild(angleForce);
    }

    // Add OutOfPlaneBend Types
    Element outOfPlaneBendForce = OutOfPlaneBendType.getXMLForce(doc, forceField);
    if (outOfPlaneBendForce != null) {
      rootElement.appendChild(outOfPlaneBendForce);
    }

    // Add Torsion Types
    Element torsionForce = TorsionType.getXMLForce(doc, forceField);
    if (torsionForce != null) {
      rootElement.appendChild(torsionForce);
    }

    // Add PiOrbitalTorsion Types
    Element piOrbitalTorsionForce = PiOrbitalTorsionType.getXMLElement(doc, forceField);
    if (piOrbitalTorsionForce != null) {
      rootElement.appendChild(piOrbitalTorsionForce);
    }

    // Add StretchTorsion Types
    Element stretchTorsionForce = StretchTorsionType.getXMLForce(doc, forceField);
    if (stretchTorsionForce != null) {
      rootElement.appendChild(stretchTorsionForce);
    }

    // Add AngleTorsion Types
    Element angleTorsionForce = AngleTorsionType.getXMLForce(doc, forceField);
    if (angleTorsionForce != null) {
      rootElement.appendChild(angleTorsionForce);
    }

    // Add StretchBend Types
    Element stretchBendForce = StretchBendType.getXMLForce(doc, forceField);
    if (stretchBendForce != null) {
      rootElement.appendChild(stretchBendForce);
    }

    // Add TorsionTorsion Type
    Element torsionTorsionForce = TorsionTorsionType.getXMLForce(doc, forceField);
    if (torsionTorsionForce != null) {
      rootElement.appendChild(torsionTorsionForce);
    }

    // Add VDW and VDWPair Types
    Element vdwForce = VDWType.getXMLForce(doc, forceField);
    if (vdwForce != null) {
      rootElement.appendChild(vdwForce);
    }

    // Add Multipole and Polarize Types
    Element multipoleForce = MultipoleType.getXMLForce(doc, forceField);
    if (multipoleForce != null) {
      rootElement.appendChild(multipoleForce);
    }

    // Add UreyBradley Types
    Element ureyBradleyForce = UreyBradleyType.getXMLForce(doc, forceField);
    if (ureyBradleyForce != null) {
      rootElement.appendChild(ureyBradleyForce);
    }

    // Write XML to 'force field name'.xml
    writeXML(doc);
  }

  /**
   * Build residues not found in OpenMM finalResidues.xml
   *
   * @param moleculeDict Map containing molecules and their BioTypes.
   * @param doc          Document with XML Nodes.
   * @param residuesNode Residues Node that will contain individual Residue nodes with Atoms and Bonds.
   */
  private static void buildExtraResidues(Map<String, ArrayList<BioType>> moleculeDict, Document doc, Element residuesNode) {
    // Loop over molecules
    for (String mol : moleculeDict.keySet()) {
      // Create Residue node for each molecule
      Element resNode = doc.createElement("Residue");
      resNode.setAttribute("name", mol);
      residuesNode.appendChild(resNode);
      List<String> order = new ArrayList<>();
      List<String> froms = new ArrayList<>();
      List<String> tos = new ArrayList<>();

      // Go through all atoms inside the molecule (keys)
      int atomCount = 0;
      int bondCount = 0;
      for (BioType bioType : moleculeDict.get(mol)) {
        String atom = bioType.atomName;
        Element atomNode = doc.createElement("Atom");
        atomNode.setAttribute("name", atom);
        atomNode.setAttribute("type", String.valueOf(bioType.atomType));
        resNode.appendChild(atomNode);
        atomCount++;
        String[] bonds = bioType.bonds;
        // Go through all atoms (bonds) that the defined atom is bonded to
        for (String bond : bonds) {
          BioType bondBioType = null;
          for (BioType bioType2 : moleculeDict.get(mol)) {
            if (bioType2.atomName.equals(bond)) {
              bondBioType = bioType2;
              break;
            }
          }

          // Record each bond and avoid duplicates by only adding the bond if the atomType is less than the bondAtomType.
          if (bioType.index < bondBioType.index) {
            froms.add(atom);
            tos.add(bond);
            bondCount++;
          }
        }
        // Add each atom in the order they have been added to the Residue node.
        order.add(atom);
      }

      // add Bond nodes using the order the Atom nodes were placed
      for (int i = 0; i < froms.size(); i++) {
        Element bondNode = doc.createElement("Bond");
        bondNode.setAttribute("from", String.valueOf(order.indexOf(froms.get(i))));
        bondNode.setAttribute("to", String.valueOf(order.indexOf(tos.get(i))));
        resNode.appendChild(bondNode);
      }
      logger.info(" Added " + atomCount + " atoms and "
          + bondCount + " bonds to residue " + mol + ".");
    }
  }

  /**
   * Create an OpenMM-style XML file from the Document object that was created in toXML().
   *
   * @param doc Document object containing XML nodes.
   * @throws TransformerException
   */
  private void writeXML(Document doc) throws TransformerException {
    String saveName;
    if (outputName != null) {
      saveName = outputName;
    } else {
      saveName = forceField.getString("forcefield", "UNKNOWN");
    }

    TransformerFactory tfFactory = TransformerFactory.newInstance();
    Transformer transformer = tfFactory.newTransformer();

    transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes"); // if want to hide first line XML version/encoding
    transformer.setOutputProperty(OutputKeys.INDENT, "yes"); // print with indentation

    DOMSource src = new DOMSource(doc);
    StreamResult result = new StreamResult(saveName + ".xml");
    transformer.transform(src, result);
  }
}
