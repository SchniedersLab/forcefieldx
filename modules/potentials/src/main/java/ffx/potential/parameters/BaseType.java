/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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
 * All force field types should extend the BaseType class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public abstract class BaseType {

    private static final Logger logger = Logger.getLogger(BaseType.class.getName());
    protected final ForceFieldType forceFieldType;
    protected String key;

    /**
     * Public constructor.
     *
     * @param forceFieldType a {@link ffx.potential.parameters.ForceField.ForceFieldType} object.
     * @param keys an array of int.
     * @since 1.0
     */
    public BaseType(ForceFieldType forceFieldType, int keys[]) {
        this.forceFieldType = forceFieldType;
        if (keys == null) {
            key = null;
        } else {
            StringBuilder keyBuffer = new StringBuilder(Integer.toString(keys[0]));
            for (int i = 1; i < keys.length; i++) {
                keyBuffer.append(" ");
                keyBuffer.append(keys[i]);
            }
            key = keyBuffer.toString();
        }
    }

    /**
     * Public constructor.
     *
     * @param forceFieldType a {@link ffx.potential.parameters.ForceField.ForceFieldType} object.
     * @param key a {@link java.lang.String} object.
     * @since 1.0
     */
    public BaseType(ForceFieldType forceFieldType, String key) {
        this.forceFieldType = forceFieldType;
        this.key = key;
    }

    /**
     * <p>Setter for the field <code>key</code>.</p>
     *
     * @param keys an array of int.
     */
    public void setKey(int[] keys) {
        StringBuilder keyBuffer = new StringBuilder(Integer.toString(keys[0]));
        for (int i = 1; i < keys.length; i++) {
            keyBuffer.append(" ");
            keyBuffer.append(keys[i]);
        }
        key = keyBuffer.toString();
    }

    /**
     * <p>Setter for the field <code>key</code>.</p>
     *
     * @param key a {@link java.lang.String} object.
     */
    public void setKey(String key) {
        this.key = key;
    }

    /**
     * Get the <code>key</code> for this Type.
     *
     * @return the key
     * @since 1.0
     */
    public String getKey() {
        return key;
    }

    /**
     * Log <code>this</code> type.
     *
     * @since 1.0
     */
    public void log() {
        if (logger.isLoggable(Level.INFO)) {
            logger.info(toString());
        }
    }

    /**
     * Print <code>this</code> to System.out.
     *
     * @since 1.0
     */
    public void print() {
        System.out.println(toString());
    }

    /**
     * {@inheritDoc}
     *
     * Basic toString method.
     * @since 1.0
     */
    @Override
    public String toString() {
        return forceFieldType + " " + key;
    }
}
