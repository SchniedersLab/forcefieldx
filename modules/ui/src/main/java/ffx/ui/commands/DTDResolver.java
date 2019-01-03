/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.ui.commands;

import java.net.URL;
import java.util.logging.Logger;

import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;

/**
 * The DTDResolver class just points the DOM DocumentBuilder to the XML Document
 * Type Definition (DTD) files.
 *
 * @author Michael J. Schnieders
 *
 */
public class DTDResolver implements EntityResolver {

    private static final Logger logger = Logger.getLogger(DTDResolver.class.getName());

    /**
     * {@inheritDoc}
     */
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
