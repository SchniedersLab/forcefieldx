/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
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

public class PDBMLFilter implements ErrorHandler {

    private File pdbxFile = null;
    static final String JAXP_SCHEMA_LANGUAGE = "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
    static final String W3C_XML_SCHEMA = "http://www.w3.org/2001/XMLSchema";
    static final String JAXP_SCHEMA_SOURCE = "http://java.sun.com/xml/jaxp/properties/schemaSource";
    static final String PDBX_XSD = "pdbx.xsd";

    /**
     * PDBML Constructor.
     */
    public PDBMLFilter() {
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
        } catch (ParserConfigurationException e) {
            System.out.println("Could not parse: " + pdbxFile + "\n" + e);
        } catch (IllegalArgumentException e) {
            System.out.println("Could not parse: " + pdbxFile + "\n" + e);
        } catch (SAXException e) {
            System.out.println("Could not parse: " + pdbxFile + "\n" + e);
        } catch (IOException e) {
            System.out.println("Could not parse: " + pdbxFile + "\n" + e);
        }
    }

    /**
     * Create an instance of the PDBML Filter
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
     */
    public static void main(String[] args) throws Exception {
        if (args.length != 1) {
            System.out.println("Usage: java PDBXFilter name.xml");
        }
        PDBMLFilter pdbx = new PDBMLFilter();
        File pdbxFile = new File(args[0]);
        pdbx.parse(pdbxFile);
    }

    /**
     * {@inheritDoc}
     */
    public void error(SAXParseException e) throws SAXException {
        System.out.println("Could not parse: " + pdbxFile + "\n" + e);
    }

    /**
     * {@inheritDoc}
     */
    public void fatalError(SAXParseException e) throws SAXException {
        System.out.println("Could not parse: " + pdbxFile + "\n" + e);
    }

    /**
     * {@inheritDoc}
     */
    public void warning(SAXParseException e) throws SAXException {
        System.out.println("Could not parse: " + pdbxFile + "\n" + e);
    }
}
