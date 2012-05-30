/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.numerics.TriCubicSpline;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>RealSpaceData class.</p>
 *
 * @author Tim Fenn
 * @version $Id: $
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
    protected double lambda = 1.0;
    // settings
    public double xweight;

    /**
     * construct a real space data assembly, assumes a real space map with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param diffractiondata {@link ffx.xray.DiffractionData diffraction data}
     */
    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties, DiffractionData diffractiondata) {
        this(new MolecularAssembly[]{assembly}, properties, diffractiondata);
    }

    /**
     * construct a real space data assembly, assumes a real space map with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param diffractiondata {@link ffx.xray.DiffractionData diffraction data}
     */
    public RealSpaceData(MolecularAssembly assembly[],
            CompositeConfiguration properties, DiffractionData diffractiondata) {
        this.assembly = assembly;
        this.modelname = assembly[0].getFile().getName();
        this.dataname = null;
        this.n = 1;
        crystal = new Crystal[n];
        crystal[0] = diffractiondata.crystal[0];
        refinementdata = new RealSpaceRefinementData[n];
        refinementdata[0] = new RealSpaceRefinementData(properties);
        refinementdata[0].periodic = true;

        xweight = properties.getDouble("xweight", 1.0);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  Real space refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        if (!diffractiondata.scaled[0]) {
            diffractiondata.scaleBulkFit();
            diffractiondata.printStats();
            // energy();
        }

        // get Fo-Fc difference density
        diffractiondata.crs_fc[0].computeAtomicGradients(diffractiondata.refinementdata[0].fofc1,
                diffractiondata.refinementdata[0].freer,
                diffractiondata.refinementdata[0].rfreeflag,
                RefinementMode.COORDINATES);

        refinementdata[0].ori[0] = 0;
        refinementdata[0].ori[1] = 0;
        refinementdata[0].ori[2] = 0;
        int extx = (int) diffractiondata.crs_fc[0].getXDim();
        int exty = (int) diffractiondata.crs_fc[0].getYDim();
        int extz = (int) diffractiondata.crs_fc[0].getZDim();
        refinementdata[0].ext[0] = extx;
        refinementdata[0].ext[1] = exty;
        refinementdata[0].ext[2] = extz;
        refinementdata[0].ni[0] = extx;
        refinementdata[0].ni[1] = exty;
        refinementdata[0].ni[2] = extz;
        refinementdata[0].data = new double[extx * exty * extz];
        for (int k = 0; k < extz; k++) {
            for (int j = 0; j < exty; j++) {
                for (int i = 0; i < extx; i++) {
                    int index1 = (i + extx * (j + exty * k));
                    int index2 = 2 * (i + extx * (j + exty * k));
                    refinementdata[0].data[index1] = diffractiondata.crs_fc[0].densityGrid[index2];
                }
            }
        }

        // now set up the refinement model
        refinementmodel = new RefinementModel(assembly);
    }

    /**
     * construct a real space data assembly, assumes a real space map with a
     * weight of 1.0 using the same name as the molecular assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     */
    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties,
                new RealSpaceFile(assembly));
    }

    /**
     * construct a real space data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link RealSpaceFile} to be refined against
     */
    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties, RealSpaceFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, datafile);
    }

    /**
     * construct a real space data assembly
     *
     * @param assembly {@link ffx.potential.bonded.MolecularAssembly molecular assembly}
     * object array (typically containing alternate conformer assemblies), used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link RealSpaceFile} to be refined against
     */
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

            if (refinementdata[i].ori[0] == 0
                    && refinementdata[i].ori[1] == 0
                    && refinementdata[i].ori[2] == 0
                    && refinementdata[i].ext[0] == refinementdata[i].ni[0]
                    && refinementdata[i].ext[1] == refinementdata[i].ni[1]
                    && refinementdata[i].ext[2] == refinementdata[i].ni[2]) {
                refinementdata[i].periodic = true;
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append("\nRefinement data settings:\n");
            sb.append("  Real space refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementmodel = new RefinementModel(assembly);
    }

    /**
     * compute real space target value, also fills in atomic derivatives
     *
     * @return target value sum over all data sets
     * @param dLambda a boolean.
     */
    public double computeRealSpaceTarget(boolean dLambda) {
        double xyz[] = new double[3];
        double uvw[] = new double[3];
        double grad[] = new double[3];
        double scalar[][][] = new double[4][4][4];

        double sum;
        for (int i = 0; i < n; i++) {
            sum = 0.0;
            TriCubicSpline spline = new TriCubicSpline();
            for (Atom a : refinementmodel.atomarray) {
                if (dLambda && !a.applyLambda()) {
                    continue;
                }

                final double lambdai = dLambda ? 1.0 : a.applyLambda() ? lambda : 1.0;
                a.getXYZ(xyz);
                a.setXYZGradient(0.0, 0.0, 0.0);
                crystal[i].toFractionalCoordinates(xyz, uvw);

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

                if (!refinementdata[i].periodic) {
                    if (ifrx - 1 < 0 || ifrx + 2 > refinementdata[i].ext[0]
                            || ifry - 1 < 0 || ifry + 2 > refinementdata[i].ext[1]
                            || ifrz - 1 < 0 || ifrz + 2 > refinementdata[i].ext[2]) {
                        logger.warning("atom: " + a.toString() + " is outside the scalar field box! Ignoring...");
                        continue;
                    }
                }

                // fill in scalar 4x4 array for interpolation
                int uii, vii, wii;
                for (int ui = ifrx - 1; ui < ifrx + 3; ui++) {
                    uii = ui - (ifrx - 1);
                    int pui = Crystal.mod(ui, refinementdata[i].ext[0]);
                    for (int vi = ifry - 1; vi < ifry + 3; vi++) {
                        vii = vi - (ifry - 1);
                        int pvi = Crystal.mod(vi, refinementdata[i].ext[1]);
                        for (int wi = ifrz - 1; wi < ifrz + 3; wi++) {
                            wii = wi - (ifrz - 1);
                            int pwi = Crystal.mod(wi, refinementdata[i].ext[2]);

                            if (refinementdata[i].periodic) {
                                scalar[uii][vii][wii] = refinementdata[i].getDataIndex(pui, pvi, pwi);
                            } else {
                                scalar[uii][vii][wii] = refinementdata[i].getDataIndex(ui, vi, wi);
                            }
                        }
                    }
                }

                // scale and interpolate
                double scale;
                if (dataname == null) {
                    scale = -1.0 * lambdai * a.getAtomType().atomicWeight;
                } else {
                    scale = -1.0 * lambdai * dataname[i].weight * a.getAtomType().atomicWeight;
                }
                double val = spline.spline(dfrx, dfry, dfrz, scalar, grad);
                sum += scale * val;

                grad[0] = grad[0] * refinementdata[i].ni[0];
                grad[1] = grad[1] * refinementdata[i].ni[1];
                grad[2] = grad[2] * refinementdata[i].ni[2];
                // transpose of toFractional
                xyz[0] = grad[0] * crystal[i].A00 + grad[1] * crystal[i].A01 + grad[2] * crystal[i].A02;
                xyz[1] = grad[0] * crystal[i].A10 + grad[1] * crystal[i].A11 + grad[2] * crystal[i].A12;
                xyz[2] = grad[0] * crystal[i].A20 + grad[1] * crystal[i].A21 + grad[2] * crystal[i].A22;
                a.addToXYZGradient(scale * xyz[0], scale * xyz[1], scale * xyz[2]);
            }
            refinementdata[i].densityscore = sum;
        }

        sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += refinementdata[i].densityscore;
        }
        return sum;
    }

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        this.lambda = lambda;
    }

    /** {@inheritDoc} */
    @Override
    public Atom[] getAtomArray() {
        return refinementmodel.atomarray;
    }

    /** {@inheritDoc} */
    @Override
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return refinementmodel.altresidues;
    }

    /** {@inheritDoc} */
    @Override
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return refinementmodel.altmolecules;
    }

    /** {@inheritDoc} */
    @Override
    public MolecularAssembly[] getMolecularAssembly() {
        return assembly;
    }

    /** {@inheritDoc} */
    @Override
    public RefinementModel getRefinementModel() {
        return refinementmodel;
    }

    /** {@inheritDoc} */
    @Override
    public double getWeight() {
        return xweight;
    }
    
    /** {@inheritDoc} */
    @Override
    public void setWeight(double weight) {
        this.xweight = weight;
    }

    /** {@inheritDoc} */
    @Override
    public String printOptimizationHeader() {
        return "Density score";
    }

    /** {@inheritDoc} */
    @Override
    public String printOptimizationUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%6.2f ", refinementdata[i].densityscore));
        }
        return sb.toString();
    }

    /** {@inheritDoc} */
    @Override
    public String printEnergyUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("     dataset %d (weight: %5.1f): chemical energy: %8.2f density score: %8.2f\n",
                    i + 1,
                    dataname[i].weight,
                    assembly[0].getPotentialEnergy().getTotalEnergy(),
                    dataname[i].weight * refinementdata[i].densityscore));
        }
        return sb.toString();
    }
}
