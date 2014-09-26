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
package ffx.potential.parsers;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.potential.bonded.MolecularAssembly;

/**
 * The FileOpener interface specifies Runnable objects which can return one or
 * more MolecularAssemblies. Implementing classes should not be constructed
 * except by a class implementing PotentialsFunctions; one should pass that
 * class the File to be opened, which constructs an implementation of
 * FileOpener.
 *
 * To some extent, this interface is legacy code of when I was trying to
 * implement the open() methods from MainPanel instead of openWait(); it should
 * be possible to simply wrap the methods into the classes implementing
 * PotentialsFunctions.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface FileOpener extends Runnable {

    @Override
    public void run();

    public MolecularAssembly getAssembly();

    public MolecularAssembly[] getAllAssemblies();

    public CompositeConfiguration getProperties();

    public CompositeConfiguration[] getAllProperties();
}
