/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential.parsers;

import org.apache.commons.configuration2.CompositeConfiguration;

import ffx.potential.MolecularAssembly;

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
 *
 */
public interface FileOpener extends Runnable {

    /** {@inheritDoc} */
    @Override
    void run();

    /**
     * <p>getAssembly.</p>
     *
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    MolecularAssembly getAssembly();

    /**
     * <p>getAllAssemblies.</p>
     *
     * @return an array of {@link ffx.potential.MolecularAssembly} objects.
     */
    MolecularAssembly[] getAllAssemblies();

    /**
     * <p>getProperties.</p>
     *
     * @return a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
     */
    CompositeConfiguration getProperties();

    /**
     * <p>getAllProperties.</p>
     *
     * @return an array of {@link org.apache.commons.configuration2.CompositeConfiguration} objects.
     */
    CompositeConfiguration[] getAllProperties();
}
