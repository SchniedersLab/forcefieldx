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
import ffx.numerics.TriCubicSpline;
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

    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties,
                new RealSpaceFile(assembly));
    }

    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties, RealSpaceFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, datafile);
    }

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
            sb.append("  Real space refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementmodel = new RefinementModel(assembly);

        /*
        CCP4MapWriter tst = new CCP4MapWriter(refinementdata[0].ori[0], refinementdata[0].ori[1],
        refinementdata[0].ori[2], refinementdata[0].ext[0], refinementdata[0].ext[1],
        refinementdata[0].ext[2], refinementdata[0].ni[0], refinementdata[0].ni[1],
        refinementdata[0].ni[2], crystal[0], "/tmp/foo.map");
        tst.setStride(1);
        tst.write(refinementdata[0].data);
         */
    }

    /**
     * compute real space target value, also fills in atomic derivatives
     * @return target value sum over all data sets
     */
    public double computeRealSpaceTarget() {
        double xyz[] = new double[3];
        double uvw[] = new double[3];
        double grad[] = new double[3];
        double scalar[][][] = new double[4][4][4];

        double sum;
        for (int i = 0; i < n; i++) {
            sum = 0.0;
            TriCubicSpline spline = new TriCubicSpline(refinementdata[i].ni[0],
                    refinementdata[i].ni[1], refinementdata[i].ni[2]);
            for (Atom a : refinementmodel.atomarray) {
                a.getXYZ(xyz);
                a.setXYZGradient(0.0, 0.0, 0.0);
                uvw[0] = xyz[0] / crystal[i].a;
                uvw[1] = xyz[1] / crystal[i].b;
                uvw[2] = xyz[2] / crystal[i].c;

                // Logic to find atom in 3d scalar field box
                final double frx = refinementdata[i].ni[0] * uvw[0];
                final int ifrx = ((int) Math.floor(frx)) - refinementdata[i].ori[0];
                final double dfrx = frx - Math.floor(frx);

                final double fry = refinementdata[i].ni[1] * uvw[1];
                final int ifry = ((int) Math.floor(fry)) - refinementdata[i].ori[1];
                final double dfry = fry - Math.floor(fry);

                final double frz = refinementdata[i].ni[2] * uvw[2];
                final int ifrz = ((int) Math.floor(frz)) - refinementdata[i].ori[2];
                final double dfrz = frz - Math.floor(frz);

                if (ifrx - 1 < 0 || ifrx + 2 > refinementdata[i].ext[0]
                        || ifry - 1 < 0 || ifry + 2 > refinementdata[i].ext[1]
                        || ifrz - 1 < 0 || ifrz + 2 > refinementdata[i].ext[2]) {
                    logger.warning("atom: " + a.toString() + " is outside the scalar field box! Ignoring...");
                    continue;
                }

                // fill in scalar 4x4 array for interpolation
                int uii, vii, wii;
                for (int ui = ifrx - 1; ui < ifrx + 3; ui++) {
                    uii = ui - (ifrx - 1);
                    for (int vi = ifry - 1; vi < ifry + 3; vi++) {
                        vii = vi - (ifry - 1);
                        for (int wi = ifrz - 1; wi < ifrz + 3; wi++) {
                            wii = wi - (ifrz - 1);

                            scalar[uii][vii][wii] = refinementdata[i].getDataIndex(ui, vi, wi);
                        }
                    }
                }

                // scale and interpolate
                double scale = -1.0 * dataname[i].weight * a.getAtomType().atomicWeight;
                double val = spline.spline(dfrx, dfry, dfrz, scalar, grad);
                sum += scale * val;
                grad[0] /= crystal[i].a;
                grad[1] /= crystal[i].b;
                grad[2] /= crystal[i].c;
                a.addToXYZGradient(scale * grad[0], scale * grad[1], scale * grad[2]);
            }
            refinementdata[i].densityscore = sum;
        }

        sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += refinementdata[i].densityscore;
        }
        return sum;
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
        return "Density score";
    }

    @Override
    public String printOptimizationUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%6.2f ", refinementdata[i].densityscore));
        }
        return sb.toString();
    }

    @Override
    public String printEnergyUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("     dataset %d (weight: %5.1f): chemical energy: %8.2f density score: %8.2f\n",
                    i + 1,
                    dataname[i].weight,
                    assembly[0].getPotentialEnergy().getTotal(),
                    dataname[i].weight * refinementdata[i].densityscore));
        }
        return sb.toString();
    }
}
