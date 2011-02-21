/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.xray;

import ffx.crystal.CCP4MapWriter;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.configuration.CompositeConfiguration;

/**
 *
 * @author Tim Fenn
 */
public class RealSpaceData implements DataContainer {

    private static final Logger logger = Logger.getLogger(RealSpaceData.class.getName());
    protected final MolecularAssembly assembly[];
    protected final RealSpaceFile dataname[];
    protected final String modelname;
    protected final int n;
    protected final Crystal crystal[];
    protected final RealSpaceRefinementData refinementdata[];
    protected final RefinementModel refinementmodel;
    // settings
    public final double xweight;

    public RealSpaceData(MolecularAssembly assembly[],
            CompositeConfiguration properties, RealSpaceFile... datafile) {

        this.assembly = assembly;
        this.modelname = assembly[0].getFile().getName();
        this.dataname = datafile;
        this.n = datafile.length;
        crystal = new Crystal[n];
        refinementdata = new RealSpaceRefinementData[n];

        xweight = properties.getDouble("xweight", 1.0);

        for (int i = 0; i < n; i++) {
            crystal[i] = datafile[i].realspacefilter.getCrystal(datafile[i].filename, properties);

            if (crystal[i] == null) {
                logger.severe("CCP4 map file does not contain full crystal information!");
            }
        }

        for (int i = 0; i < n; i++) {
            refinementdata[i] = new RealSpaceRefinementData(properties);
            datafile[i].realspacefilter.readFile(datafile[i].filename,
                    refinementdata[i], properties);
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  X-ray refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementmodel = new RefinementModel(assembly);




        CCP4MapWriter tst = new CCP4MapWriter(refinementdata[0].ext[0], refinementdata[0].ext[1], refinementdata[0].ext[2],
                crystal[0], "/tmp/foo.map");
        tst.setStride(1);
        tst.write(refinementdata[0].data);
    }

    @Override
    public Atom[] getAtomArray() {
        return refinementmodel.atomarray;
    }

    @Override
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return refinementmodel.altresidues;
    }

    @Override
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return refinementmodel.altmolecules;
    }

    @Override
    public MolecularAssembly[] getMolecularAssembly() {
        return assembly;
    }

    @Override
    public RefinementModel getRefinementModel() {
        return refinementmodel;
    }

    @Override
    public double getWeight() {
        return xweight;
    }

    @Override
    public String printOptimizationHeader() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String printOptimizationUpdate() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String printEnergyUpdate() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
