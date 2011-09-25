
/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Engineering Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 * <p>Institution: Labs of Axel T. Brunger (Stanford), Vijay S. Pande (Stanford) and Jay W. Ponder (WUSTL)</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
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

	/** {@inheritDoc} */
	public void error(SAXParseException e) throws SAXException {
		System.out.println("Could not parse: " + pdbxFile + "\n" + e);
	}

	/** {@inheritDoc} */
	public void fatalError(SAXParseException e) throws SAXException {
		System.out.println("Could not parse: " + pdbxFile + "\n" + e);
	}

	/** {@inheritDoc} */
	public void warning(SAXParseException e) throws SAXException {
		System.out.println("Could not parse: " + pdbxFile + "\n" + e);
	}
}
