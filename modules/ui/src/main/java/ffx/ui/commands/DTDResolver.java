/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui.commands;

import java.net.URL;
import java.util.logging.Logger;

import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;

/**
 * The DTDResolver class just points the DOM DocumentBuilder to the XML Document
 * Type Definition files.
 */
public class DTDResolver implements EntityResolver {
	private static final Logger logger = Logger.getLogger(DTDResolver.class.getName());

	public InputSource resolveEntity(String publicId, String systemId) {
		if (systemId.lastIndexOf("keywords.dtd") >= 0) {
			URL keyURL = getClass().getClassLoader().getResource(
					"ffx/xml/keywords.dtd");
			try {
				return new InputSource(keyURL.openStream());
			} catch (Exception e) {
				logger.warning("" + e);
				return null;
			}
		} else if (systemId.lastIndexOf("commands.dtd") >= 0) {
			URL commandURL = getClass().getClassLoader().getResource(
					"ffx/xml/commands.dtd");
			try {
				return new InputSource(commandURL.openStream());
			} catch (Exception e) {
				return null;
			}
		}
		return null;
	}
}
