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
package test

import ffx.potential.bonded.MSNode

import java.util.stream.Collectors
import java.util.Arrays;

import org.apache.commons.io.FilenameUtils

import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

/**
 * The PACCOM script calculates intracrystollographic distances.
 *
 * created by Okimasa OKADA 2017/3/31
 * revised by Okimasa OKADA 2019/2/25
 * ported to FFX by Aaron Nessler and Micheal Schnieders 2020
 * <br>
 * Usage:
 * <br>
 * ffxc test.PACCOM &lt;filename&gt;
 */
@Command(description = " Compare crystal packings based on intermolecular distances.", name = "ffxc test.PACCOM")
class PACCOM extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null
    private MolecularAssembly[] assemblies;

    public double[][] coordinates = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    PACCOM run() {

        // The following should become flags at some point:
        // PSUEDO: PARAMETERS INITIALIZED IN ORIGINAL PACCOM

        int n1 = 2; // 1st atom to fit
        int n2 = 6; // 2nd atom to fit
        int n3 = 7; // 3rd atom to fit
        int n_mol_cut = 32; //Number of molecules for RMSD calc
        //int n_mol = 32; //Number of molecules for fitting
        int par1 = 8; // parameter for second round fitting (default is probably OK)
        int par2 = 6; // parameter for second round fitting (default is probably OK)
        int par3 = 8; // parameter for second round fitting (default is probably OK)

        // PSUEDO: READ IN STRUCTURES

        if (!init()) {
            return
        }
        // Ensure file exists/create active molecular assembly
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.openAll(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            assemblies = [activeAssembly]
        }
        // Get file
        String modelFilename = activeAssembly.getFile().getAbsolutePath()
        logger.info("\n Calculating distance for: " + modelFilename)
        // PSUEDO: Now that the molecular assembly is obtained, we will center each system to the origin and
        //  record distances from each molecule to the origin for each crystal
        double[] origin = [0.0,0.0,0.0]
        double[][] distance1 = new double[assemblies.size()][assemblies[0].getMolecules()[0].getAtomList().size()];
        int sysIndex = 0;
        // Loop over each system.
        for (int i = 0; i < assemblies.length; i++) {
            def system = assemblies[i]
            system.moveCenter(origin)
            Crystal crystal = system.getCrystal().getUnitCell()
            List<MSNode> molecules = system.getMolecules();
            List<Atom> atoms = system.getAtomList()

            // Sum over molecules to determine center most in system

            coordinates = new double[atoms.size()][3]
            double[] cart = new double[3]
            double[] frac = new double[3]

            //TODO: Currently assumes cartesian coordinates input (Masa's did the same). Potential flag to support fracitonal?
            int molIndex = 0;
            for (MSNode molecule in molecules) {
                double sx = 0.0;
                double sy = 0.0;
                double sz = 0.0;
                for (Atom atom in atoms) {
                    atom.getXYZ(cart)
                    crystal.toFractionalCoordinates(cart, frac)
                    atom.moveTo(frac)

                    sx+=frac[0]; // Summation of x points
                    sy+=frac[1]; // Summation of y points
                    sz+=frac[2]; // Summation of z points
                }
                sx=sx/atoms.size();
                sy=sy/atoms.size();
                sz=sz/atoms.size();
                //Record each average molecular distance for each crystal system
                distance1[sysIndex][molIndex++] = Math.sqrt(Math.pow(sx, 2) +  Math.pow(sy, 2) + Math.pow(sz, 2));
            }
            sysIndex++;
        }

        //TODO: PSEUDO: SORT BY DISTANCE TO CENTER OF ALL ATOMS TO DETERMINE MOLECULE TO MOVE TO ORIGIN

        for (double [] molDist in distance1){
            Arrays.sort(molDist);
        }

        //TODO: PSUEDO: MOVE CENTER-MOST MOLECULE OF FIRST LIST TO ORIGIN. ORIENTED n1 AT ORIGIN
        //TODO: PSUEDO: RETAIN CLOSEST MOLECULES (# DETERMINED BY n_mol_cut) CREATES "SHELL"
        //TODO: PSUEDO: MOVE n1 ATOM OF mTH MOLECULE OF 2ND LIST TO ORIGIN. m defined at line 325 of comp_rmsd_33_ffx_L.f
        //TODO: PSUEDO: ORIENT MOLECULES SO CENTER-MOST MOLECULE ON XY PLANE
        // PSUEDO: AT THIS POINT THE 1ST TWO MOLECULES IN THE TWO LISTS ARE SUPERPOSED (origin, x-axis, and x-y plane)
        //TODO: PSUEDO: nm (USUALLY 3) MOLECULES NEAR ORIGIN ARE USED. 1ST MOLECULE IS REPLACED BY mTH MOLECULE
        //TODO: PSUEDO: MOLECULES IN 2ND LIST ARRANGED TO MATCH 1ST LIST BASED ON DISTANCE (CLOSEST nm MOLECULES BY n1 ATOM).
        // PSUEDO: AT THIS POINT FIRST nm MOLECULES SHOULD BE CLOSE TO EACH OTHER.
        //TODO: PSUEDO: USE n1 ATOMS OF FIRST nm MOLECULES IN EACH LIST TO SUPERIMPOSE MOLECULES OF LISTS.
        //TODO: PSUEDO: PAIR MOLECULES OF EACH LIST FOR RMSD
        //TODO: PSUEDO: CALL RMSD METHOD DEVELOPED BY KEN DILL

// Rest comes from Cart2Frac which this script was based off of. Will need to edit as necessary...

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }

        String fileName = FilenameUtils.getName(modelFilename)
        String ext = FilenameUtils.getExtension(fileName)
        fileName = FilenameUtils.removeExtension(fileName)

        String dirName = FilenameUtils.getFullPath(saveDir.getAbsolutePath())

        if (ext.toUpperCase().contains("XYZ")) {
            potentialFunctions.saveAsXYZ(assemblies[0], new File(dirName + fileName + ".xyz"))
        } else {
            potentialFunctions.saveAsPDB(assemblies, new File(dirName + fileName + ".pdb"))
        }

        return this
    }

    @Override
    public List<Potential> getPotentials() {
        if (assemblies == null) {
            return new ArrayList<Potential>();
        } else {
            return Arrays.stream(assemblies).
                    filter { a -> a != null }.
                    map { a -> a.getPotentialEnergy() }.
                    filter { e -> e != null }.
                    collect(Collectors.toList());
        }
    }
}
