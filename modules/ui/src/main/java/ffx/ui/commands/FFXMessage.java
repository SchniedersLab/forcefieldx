/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.ui.commands;

import java.io.Serializable;

/**
 * The FFXMessage class is used to pass simple messages between a TinkerServer
 * and its FFXClient(s).
 *
 * @author Michael J. Schnieders
 *
 */
public class FFXMessage implements Serializable {

    private static final long serialVersionUID = 1L;
    /**
     * Constant <code>SYSTEM=0</code>
     */
    public static int SYSTEM = 0;
    /**
     * Constant <code>UPDATE=1</code>
     */
    public static int UPDATE = 1;
    /**
     * Constant <code>CLOSING=2</code>
     */
    public static int CLOSING = 2;
    /**
     * Constant <code>OK=3</code>
     */
    public static int OK = 3;
    private int message = 0;
    private int step = -1;
    private double time = -1.0;
    private int type = 0;

    /**
     * <p>
     * Constructor for FFXMessage.</p>
     *
     * @param m a int.
     */
    public FFXMessage(int m) {
        message = m;
    }

    /**
     * <p>
     * Getter for the field <code>message</code>.</p>
     *
     * @return a int.
     */
    public int getMessage() {
        return message;
    }

    /**
     * <p>
     * Getter for the field <code>step</code>.</p>
     *
     * @return a int.
     */
    public int getStep() {
        return step;
    }

    /**
     * <p>
     * Getter for the field <code>time</code>.</p>
     *
     * @return a double.
     */
    public double getTime() {
        return time;
    }

    /**
     * <p>
     * Getter for the field <code>type</code>.</p>
     *
     * @return a int.
     */
    public int getType() {
        return type;
    }

    /**
     * <p>
     * print</p>
     */
    public void print() {
        System.out.println(toString());
    }

    /**
     * <p>
     * Setter for the field <code>message</code>.</p>
     *
     * @param m a int.
     */
    public void setMessage(int m) {
        message = m;
    }

    /**
     * <p>
     * Setter for the field <code>step</code>.</p>
     *
     * @param s a int.
     */
    public void setStep(int s) {
        step = s;
    }

    /**
     * <p>
     * Setter for the field <code>time</code>.</p>
     *
     * @param t a double.
     */
    public void setTime(double t) {
        time = t;
    }

    /**
     * <p>
     * Setter for the field <code>type</code>.</p>
     *
     * @param t a int.
     */
    public void setType(int t) {
        type = t;
    }

    /**
     * <p>
     * toString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toString() {
        if (message == 0) {
            return new String("SYSTEM").intern();
        } else if (message == 1) {
            return new String("UPDATE").intern();
        } else if (message == 2) {
            return new String("CLOSING").intern();
        } else {
            return new String("OK").intern();
        }
    }
}
