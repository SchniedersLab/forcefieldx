//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.xray;

import java.io.File;
import java.util.List;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.configuration2.CompositeConfiguration;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.utils.PotentialsFileOpener;
import ffx.utilities.Keyword;
import ffx.xray.parsers.MTZFilter;

/**
 * <p>TimerTest class.</p>
 *
 * @author Timothy D. Fenn
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class TimerTest {

    private static final Logger logger = Logger.getLogger(TimerTest.class.getName());

    /**
     * <p>main.</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        // Parameters collection from original Timer script
        String pdbname = System.getProperty("pdbFile", "1N7S.pdb");
        String mtzname = System.getProperty("mtzFile", null);
        final boolean ciOnly = false;

        boolean ci = System.getProperty("ffx.ci", "false").equalsIgnoreCase("true");
        if (!ci && ciOnly) {
            return;
        }

        // Load the structure
        MolecularAssembly molecularAssembly;
        File structure = new File(pdbname);
        PotentialsFileOpener opener = new PotentialsFileOpener(structure);
        opener.run();
        molecularAssembly = opener.getAssembly();
        File mtzFile = new File(mtzname);

        // Load any properties associated with it
        CompositeConfiguration properties = Keyword.loadProperties(structure);

        // Read in Fo/sigFo/FreeR
        MTZFilter mtzFilter = new MTZFilter();
        Crystal crystal = Crystal.checkProperties(properties);
        Resolution resolution = Resolution.checkProperties(properties);
        ReflectionList reflectionList;
        if (crystal == null || resolution == null) {
            reflectionList = mtzFilter.getReflectionList(mtzFile);
        } else {
            reflectionList = new ReflectionList(crystal, resolution);
        }

        DiffractionRefinementData refinementData = new DiffractionRefinementData(properties, reflectionList);
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();

        // Associate molecular assembly with the structure, set up forcefield
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbFile = new PDBFilter(structure, molecularAssembly, forceField, properties);
        pdbFile.readFile();
        pdbFile.applyAtomProperties();
        molecularAssembly.finalize(true, forceField);
        ForceFieldEnergy.energyFactory(molecularAssembly, pdbFile.getCoordRestraints());

        List<Atom> atomList = molecularAssembly.getAtomList();
        Atom[] atomArray = atomList.toArray(new Atom[0]);

        // set up FFT and run it
        ParallelTeam parallelTeam = new ParallelTeam();
        CrystalReciprocalSpace crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, false);
        crs.computeDensity(refinementData.fc);
        refinementData.setCrystalReciprocalSpace_fc(crs);
        crs = new CrystalReciprocalSpace(reflectionList,
                atomArray, parallelTeam, parallelTeam, true);
        crs.computeDensity(refinementData.fs);
        refinementData.setCrystalReciprocalSpace_fs(crs);

        ScaleBulkMinimize scaleBulkMinimize = new ScaleBulkMinimize(reflectionList, refinementData, crs, parallelTeam);
        scaleBulkMinimize.minimize(6, 1.0e-4);

        SigmaAMinimize sigmaAMinimize = new SigmaAMinimize(reflectionList, refinementData, parallelTeam);
        sigmaAMinimize.minimize(7, 2.0e-2);

        SplineMinimize splineMinimize = new SplineMinimize(reflectionList, refinementData, refinementData.spline, SplineEnergy.Type.FOFC);
        splineMinimize.minimize(7, 1e-5);

        CrystalStats crystalStats = new CrystalStats(reflectionList, refinementData);

        scaleBulkMinimize = new ScaleBulkMinimize(reflectionList, refinementData, refinementData.crs_fs, parallelTeam);
        ScaleBulkEnergy scaleBulkEnergy = scaleBulkMinimize.getScaleBulkEnergy();
        int n = scaleBulkMinimize.getNumberOfVariables();
        double[] x = new double[n];
        double[] g = new double[n];
        scaleBulkMinimize.getCoordinates(x);
        scaleBulkEnergy.energyAndGradient(x, g);

        logger.info("SCATTER TEST");
        for (int i = 0; i < 30; i++) {
            long time = -System.nanoTime();
            scaleBulkEnergy.energyAndGradient(x, g);
            time += System.nanoTime();
            logger.info(format(" Time %12.8f", time * 1.0e-9));
        }
    }
}
