package ffx.xray.groovy

import ffx.numerics.Potential
import ffx.xray.DiffractionRefinementData
import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.crystal.Crystal
import ffx.crystal.ReflectionList
import ffx.crystal.Resolution
import ffx.potential.MolecularAssembly
import ffx.xray.parsers.CIFFilter
import ffx.xray.parsers.MTZWriter
import ffx.xray.parsers.MTZWriter.MTZType

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The CIF2MTZ script saves a CIF file to MTZ format.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.CIF2MTZ &lt;filename&gt;
 */
@Command(description = " Convert a CIF file to MTZ format.", name = "ffxc xray.CIFMTZ")
class CIF2MTZ extends AlgorithmsScript {

    /**
     * A CIF filename.
     */
    @Parameters(arity = "2", paramLabel = "file", description = "A PDB file and a CIF diffraction file.")
    private ArrayList<String> filenames = null

    private MolecularAssembly[] systems;
    private DiffractionRefinementData refinementdata;

    /**
     * Execute the script.
     */
    @Override
    CIF2MTZ run() {

        if (!init()) {
            return this
        }

        String pdb = filenames.get(0)
        String cif = filenames.get(1)

        logger.info("\n Running CIF2MTZ on " + cif)

        // Use PotentialsFunctions methods instead of Groovy method closures to do work.
        systems = algorithmFunctions.open(pdb)

        CIFFilter ciffilter = new CIFFilter()
        ReflectionList reflectionlist = ciffilter.getReflectionList(new File(cif), systems[0].getProperties())

        if (reflectionlist == null) {
            println(" Using crystal information from the PDB file to generate MTZ file.")

            Crystal crystal = systems[0].getCrystal().getUnitCell()
            double res = ciffilter.getResolution(new File(cif), crystal)
            if (res < 0.0) {
                println(" Resolution could not be determined from the PDB and CIF files.")
                return this
            }

            Resolution resolution = new Resolution(res)
            reflectionlist = new ReflectionList(crystal, resolution, systems[0].getProperties())
        }

        refinementdata = new DiffractionRefinementData(systems[0].getProperties(), reflectionlist)
        ciffilter.readFile(new File(cif), reflectionlist, refinementdata, systems[0].getProperties())

        MTZWriter mtzwriter = new MTZWriter(reflectionlist, refinementdata, FilenameUtils.removeExtension(cif) + ".mtz", MTZType.DATAONLY)
        mtzwriter.write()

        return this
    }

    @Override
    public List<Potential> getPotentials() {
        if (systems == null) {
            return new ArrayList<Potential>();
        } else {
            return Arrays.stream(systems).filter { a -> a != null }.
                    map { a -> a.getPotentialEnergy() }.
                    filter { e -> e != null }.
                    collect(Collectors.toList());
        }
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
