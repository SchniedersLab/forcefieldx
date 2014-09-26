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
package ffx.ui;

/**
 * The UIFileCloser class wraps the closing of an FFXSystem within a thread.
 *
 * @author Michael J. Schnieders
 *
 */
public final class UIFileCloser implements Runnable {

    FFXSystem ffxSystem;

    /**
     * <p>
     * Constructor for FileCloser.</p>
     *
     * @param sys a {@link ffx.ui.FFXSystem} object.
     */
    public UIFileCloser(FFXSystem sys) {
        ffxSystem = sys;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        ffxSystem.destroy();
    }
}
