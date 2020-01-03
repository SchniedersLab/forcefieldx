//******************************************************************************
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
//******************************************************************************
package ffx.potential.utils;

/**
 * This Exception class indicates an error in calculating energy or gradients.
 * Expected behavior is that it will be caught by Potential.energy(), resulting
 * in any necessary cleanup. Then, if the causeSevere flag is set true, FFE will
 * issue a logger.severe (resulting in exit); else, FFE will simply rethrow the
 * exception. The default is to rethrow the exception.
 *
 * @author Jacob Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class EnergyException extends ArithmeticException {

    private final boolean causeSevere;

    /**
     * <p>Constructor for EnergyException.</p>
     */
    public EnergyException() {
        super();
        causeSevere = false;
    }

    /**
     * <p>Constructor for EnergyException.</p>
     *
     * @param str a {@link java.lang.String} object.
     */
    public EnergyException(String str) {
        super(str);
        causeSevere = false;
    }

    /**
     * <p>Constructor for EnergyException.</p>
     *
     * @param causeSevere a boolean.
     */
    public EnergyException(boolean causeSevere) {
        super();
        this.causeSevere = causeSevere;
    }

    /**
     * <p>Constructor for EnergyException.</p>
     *
     * @param str         a {@link java.lang.String} object.
     * @param causeSevere a boolean.
     */
    public EnergyException(String str, boolean causeSevere) {
        super(str);
        this.causeSevere = causeSevere;
    }

    /**
     * <p>doCauseSevere.</p>
     *
     * @return a boolean.
     */
    public boolean doCauseSevere() {
        return causeSevere;
    }
}
