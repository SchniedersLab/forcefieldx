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
package ffx.potential.parameters;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.parameters.ForceField.ForceFieldType;

/**
 * The BaseType class.
 */
public abstract class BaseType {

    private static final Logger log = Logger.getLogger(BaseType.class.getName());

    public final ForceFieldType forceFieldType;
    public final String key;

    public BaseType(ForceFieldType forceFieldType, int keys[]) {
        this.forceFieldType = forceFieldType;
        if (keys == null) {
            key = null;
        } else {
            StringBuffer keyBuffer = new StringBuffer(Integer.toString(keys[0]));
            for (int i = 1; i < keys.length; i++) {
                keyBuffer.append(" " + keys[i]);
            }
            key = keyBuffer.toString();
        }
    }

    public BaseType(ForceFieldType forceFieldType, String key) {
        this.forceFieldType = forceFieldType;
        this.key = key;
    }

    public void log() {
        if (log.isLoggable(Level.FINE)) {
            log.fine(this.toString());
        }
    }

    /**
     * Print the Type to System.out.
     */
    public void print() {
        log.info(toString());
    }
}
