// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.parsers;

import java.io.File;
import java.io.IOException;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

/**
 * PDBMLFilter class.
 *
 * @author Michael J. Schnieders
 */
public class PDBMLFilter implements ErrorHandler {

  private static final String JAXP_SCHEMA_LANGUAGE =
      "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
  private static final String W3C_XML_SCHEMA = "http://www.w3.org/2001/XMLSchema";
  private static final String JAXP_SCHEMA_SOURCE =
      "http://java.sun.com/xml/jaxp/properties/schemaSource";
  private static final String PDBX_XSD = "pdbx.xsd";
  private File pdbxFile = null;

  /** PDBML Constructor. */
  public PDBMLFilter() {}

  /**
   * Create an instance of the PDBML Filter
   *
   * @param args an array of {@link java.lang.String} objects.
   */
  public static void main(String[] args) {
    if (args.length != 1) {
      System.out.println("Usage: java PDBXFilter name.xml");
    }
    PDBMLFilter pdbx = new PDBMLFilter();
    File pdbxFile = new File(args[0]);
    pdbx.parse(pdbxFile);
  }

  /** {@inheritDoc} */
  @Override
  public void error(SAXParseException e) throws SAXException {
    parseMessage(e);
  }

  /** {@inheritDoc} */
  @Override
  public void fatalError(SAXParseException e) throws SAXException {
    parseMessage(e);
  }

  /**
   * parseMessage.
   *
   * @param e a {@link java.lang.Exception} object.
   */
  public void parseMessage(Exception e) {
    System.out.println("Could not parse: " + pdbxFile + "\n" + e);
  }

  /** {@inheritDoc} */
  @Override
  public void warning(SAXParseException e) throws SAXException {
    parseMessage(e);
  }

  private void parse(File pdbxml) {
    pdbxFile = pdbxml;
    try {
      DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
      dbf.setNamespaceAware(true);
      dbf.setValidating(true);
      dbf.setAttribute(JAXP_SCHEMA_LANGUAGE, W3C_XML_SCHEMA);
      dbf.setAttribute(JAXP_SCHEMA_SOURCE, new File(PDBX_XSD));
      DocumentBuilder db = dbf.newDocumentBuilder();
      db.setErrorHandler(this);
      Document doc = db.parse(pdbxFile);
      NodeList list = doc.getElementsByTagName("PDBx:atom_site");
      System.out.println("Number of Atoms: " + list.getLength());
    } catch (ParserConfigurationException
        | IllegalArgumentException
        | SAXException
        | IOException e) {
      parseMessage(e);
    }
  }
}
