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
package ffx.potential.parameters;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.potential.parameters.ForceField.ForceFieldType;

/**
 * All force field types should extend the BaseType class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class BaseType {

    private static final Logger logger = Logger.getLogger(BaseType.class.getName());

    /**
     * The ForceFieldType of this term.
     */
    final ForceFieldType forceFieldType;
    /**
     * The look-up key for this term, which is usually a concatenation of atom classes or atom types.
     */
    protected String key;

    /**
     * Public constructor.
     *
     * @param forceFieldType a {@link ffx.potential.parameters.ForceField.ForceFieldType} object.
     * @param keys           an array of int.
     * @since 1.0
     */
    public BaseType(ForceFieldType forceFieldType, int[] keys) {
        this.forceFieldType = forceFieldType;
        setKey(keys);
    }

    /**
     * Public constructor.
     *
     * @param forceFieldType a {@link ffx.potential.parameters.ForceField.ForceFieldType} object.
     * @param key            a {@link java.lang.String} object.
     * @since 1.0
     */
    public BaseType(ForceFieldType forceFieldType, String key) {
        this.forceFieldType = forceFieldType;
        this.key = key;
    }

    /**
     * <p>
     * Setter for the field <code>key</code>.</p>
     *
     * @param keys an array of int.
     */
    public void setKey(int[] keys) {
        if (keys == null) {
            key = null;
            return;
        }

        StringBuilder keyBuffer = new StringBuilder();
        for (int k : keys) {
            keyBuffer.append(k);
            keyBuffer.append(" ");
        }
        key = keyBuffer.toString().trim();
    }

    /**
     * <p>
     * Setter for the field <code>key</code>.</p>
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
     * <p>
     * Basic toString method.
     *
     * @since 1.0
     */
    @Override
    public String toString() {
        return forceFieldType + " " + key;
    }
}
