/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.xray;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.crystal.Crystal;
import ffx.numerics.TriCubicSpline;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Residue;
import ffx.xray.RefinementMinimize.RefinementMode;

/**
 * <p>
 * RealSpaceData class.</p>
 *
 * @author Timothy D. Fenn
 *
 */
public class RealSpaceData implements DataContainer {

    private static final Logger logger = Logger.getLogger(RealSpaceData.class.getName());

    /**
     * The collection of MolecularAssembly instances that model the real space
     * data.
     */
    protected final MolecularAssembly[] molecularAssemblies;
    /**
     * The real space data files.
     */
    protected final RealSpaceFile[] realSpaceFile;
    /**
     * The name of the model.
     */
    protected final String modelName;
    /**
     * The number of real space data sources to consider.
     */
    protected final int n;
    /**
     * The collection of crystals that describe the PBC and symmetry of each
     * data source.
     */
    protected final Crystal crystal[];
    /**
     * The real space refinement data for each data source.
     */
    protected final RealSpaceRefinementData[] refinementData;
    /**
     * The refinement model.
     */
    protected final RefinementModel refinementModel;
    /**
     * The total real space energy.
     */
    private double realSpaceEnergy = 0.0;
    /**
     * The partial derivative of the real space energy with respect to lambda.
     */
    private double realSpacedUdL = 0.0;
    /**
     * The real space gradient.
     */
    private double realSpaceGradient[];
    /**
     * The partial derivative of the real space gradient with respect to lambda.
     */
    private double realSpacedUdXdL[];
    /**
     * The weight of the data relative to the weight of the force field.
     */
    public double xweight;
    /**
     * The current value of the state variable lambda.
     */
    protected double lambda = 1.0;
    /**
     * A flag to indicate if the lambda term is active.
     */
    private boolean lambdaTerm = false;

    /**
     * Construct a real space data molecularAssemblies, assumes a real space map
     * with a weight of 1.0 using the same name as the molecular
     * molecularAssemblies.
     *
     * @param molecularAssembly
     * {@link ffx.potential.MolecularAssembly molecular molecularAssemblies}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param diffractionData {@link ffx.xray.DiffractionData diffraction data}
     */
    public RealSpaceData(MolecularAssembly molecularAssembly,
            CompositeConfiguration properties, DiffractionData diffractionData) {
        this(new MolecularAssembly[]{molecularAssembly}, properties, diffractionData);
    }

    /**
     * Construct a real space data molecularAssemblies, assumes a real space map
     * with a weight of 1.0 using the same name as the molecularAssemblies.
     *
     * @param molecularAssemblies
     * {@link ffx.potential.MolecularAssembly molecular molecularAssemblies}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param diffractionData {@link ffx.xray.DiffractionData diffraction data}
     */
    public RealSpaceData(MolecularAssembly molecularAssemblies[],
            CompositeConfiguration properties, DiffractionData diffractionData) {
        this.molecularAssemblies = molecularAssemblies;
        this.modelName = molecularAssemblies[0].getFile().getName();
        this.realSpaceFile = null;
        this.n = 1;
        crystal = new Crystal[n];
        crystal[0] = diffractionData.crystal[0];
        refinementData = new RealSpaceRefinementData[n];
        refinementData[0] = new RealSpaceRefinementData();
        refinementData[0].periodic = true;

        xweight = properties.getDouble("xweight", 1.0);
        lambdaTerm = properties.getBoolean("lambdaterm", false);

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Refinement Settings\n");
            sb.append("  Real space refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        if (!diffractionData.scaled[0]) {
            diffractionData.scaleBulkFit();
            diffractionData.printStats();
        }

        // get Fo-Fc difference density
        diffractionData.crs_fc[0].computeAtomicGradients(diffractionData.refinementData[0].fofc1,
                diffractionData.refinementData[0].freer,
                diffractionData.refinementData[0].rfreeflag,
                RefinementMode.COORDINATES);

        refinementData[0].ori[0] = 0;
        refinementData[0].ori[1] = 0;
        refinementData[0].ori[2] = 0;
        int extx = (int) diffractionData.crs_fc[0].getXDim();
        int exty = (int) diffractionData.crs_fc[0].getYDim();
        int extz = (int) diffractionData.crs_fc[0].getZDim();
        refinementData[0].ext[0] = extx;
        refinementData[0].ext[1] = exty;
        refinementData[0].ext[2] = extz;
        refinementData[0].ni[0] = extx;
        refinementData[0].ni[1] = exty;
        refinementData[0].ni[2] = extz;
        refinementData[0].data = new double[extx * exty * extz];
        for (int k = 0; k < extz; k++) {
            for (int j = 0; j < exty; j++) {
                for (int i = 0; i < extx; i++) {
                    int index1 = (i + extx * (j + exty * k));
                    int index2 = 2 * (i + extx * (j + exty * k));
                    refinementData[0].data[index1] = diffractionData.crs_fc[0].densityGrid[index2];
                }
            }
        }

        // now set up the refinement model
        refinementModel = new RefinementModel(molecularAssemblies);
    }

    /**
     * Construct a real space data molecularAssemblies, assumes a real space map
     * with a weight of 1.0 using the same name as the molecular
     * molecularAssemblies.
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular molecularAssemblies}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     */
    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties) {
        this(new MolecularAssembly[]{assembly}, properties,
                new RealSpaceFile(assembly));
    }

    /**
     * Construct a real space data molecularAssemblies.
     *
     * @param assembly
     * {@link ffx.potential.MolecularAssembly molecular molecularAssemblies}
     * object, used as the atomic model for comparison against the data
     * @param properties system properties file
     * @param datafile one or more {@link RealSpaceFile} to be refined against
     */
    public RealSpaceData(MolecularAssembly assembly,
            CompositeConfiguration properties, RealSpaceFile... datafile) {
        this(new MolecularAssembly[]{assembly}, properties, datafile);
    }

    /**
     * Construct a real space data molecularAssemblies.
     *
     * @param molecularAssemblies
     * {@link ffx.potential.MolecularAssembly molecular molecularAssemblies}
     * object array (typically containing alternate conformer assemblies), used
     * as the atomic model for comparison against the data
     * @param properties system properties file
     * @param dataFile one or more {@link RealSpaceFile} to be refined against
     */
    public RealSpaceData(MolecularAssembly molecularAssemblies[],
            CompositeConfiguration properties, RealSpaceFile... dataFile) {

        this.molecularAssemblies = molecularAssemblies;
        this.modelName = molecularAssemblies[0].getFile().getName();
        this.realSpaceFile = dataFile;
        this.n = dataFile.length;
        crystal = new Crystal[n];
        refinementData = new RealSpaceRefinementData[n];

        xweight = properties.getDouble("xweight", 1.0);
        lambdaTerm = properties.getBoolean("lambdaterm", false);

        for (int i = 0; i < n; i++) {
            crystal[i] = dataFile[i].realSpaceFileFilter.getCrystal(dataFile[i].filename, properties);

            if (crystal[i] == null) {
                logger.severe("CCP4 map file does not contain full crystal information!");
            }
        }

        for (int i = 0; i < n; i++) {
            refinementData[i] = new RealSpaceRefinementData();
            dataFile[i].realSpaceFileFilter.readFile(dataFile[i].filename,
                    refinementData[i], properties);

            if (refinementData[i].ori[0] == 0
                    && refinementData[i].ori[1] == 0
                    && refinementData[i].ori[2] == 0
                    && refinementData[i].ext[0] == refinementData[i].ni[0]
                    && refinementData[i].ext[1] == refinementData[i].ni[1]
                    && refinementData[i].ext[2] == refinementData[i].ni[2]) {
                refinementData[i].periodic = true;
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            sb.append(" Refinement Settings\n");
            sb.append("  Real space refinement weight (xweight): " + xweight + "\n");
            logger.info(sb.toString());
        }

        // now set up the refinement model
        refinementModel = new RefinementModel(molecularAssemblies);
    }

    /**
     * compute real space target value, also fills in atomic derivatives
     *
     * @return target value sum over all data sets.
     */
    public double computeRealSpaceTarget() {
        double xyz[] = new double[3];
        double uvw[] = new double[3];
        double grad[] = new double[3];
        double scalar[][][] = new double[4][4][4];

        // Zero out the realSpaceTarget energy.
        realSpaceEnergy = 0.0;
        // Zero out the realSpacedUdL energy.
        realSpacedUdL = 0.0;
        // Initialize gradient to zero; allocate space if necessary.
        int nUsedAtoms = refinementModel.usedAtoms.length;
        int nActive = refinementModel.activeAtoms.length;

        int nGrad = nActive * 3;
        if (realSpaceGradient == null || realSpaceGradient.length < nGrad) {
            realSpaceGradient = new double[nGrad];
        } else {
            Arrays.fill(realSpaceGradient, 0.0);
        }
        // Initialize dUdXdL to zero; allocate space if necessary.
        if (realSpacedUdXdL == null || realSpacedUdXdL.length < nGrad) {
            realSpacedUdXdL = new double[nGrad];
        } else {
            Arrays.fill(realSpacedUdXdL, 0.0);
        }

        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            TriCubicSpline spline = new TriCubicSpline();
            int index = 0;
            for (Atom a : refinementModel.usedAtoms) {
                double lambdai = 1.0;
                double dUdL = 0.0;
                if (lambdaTerm && a.applyLambda()) {
                    lambdai = lambda;
                    dUdL = 1.0;
                    a.setLambdaXYZGradient(0.0, 0.0, 0.0);
                }
                a.getXYZ(xyz);
                a.setXYZGradient(0.0, 0.0, 0.0);
                crystal[i].toFractionalCoordinates(xyz, uvw);

                // Logic to find atom in 3d scalar field box
                final double frx = refinementData[i].ni[0] * uvw[0];
                final int ifrx = ((int) Math.floor(frx)) - refinementData[i].ori[0];
                final double dfrx = frx - Math.floor(frx);

                final double fry = refinementData[i].ni[1] * uvw[1];
                final int ifry = ((int) Math.floor(fry)) - refinementData[i].ori[1];
                final double dfry = fry - Math.floor(fry);

                final double frz = refinementData[i].ni[2] * uvw[2];
                final int ifrz = ((int) Math.floor(frz)) - refinementData[i].ori[2];
                final double dfrz = frz - Math.floor(frz);

                if (!refinementData[i].periodic) {
                    if (ifrx - 1 < 0 || ifrx + 2 > refinementData[i].ext[0]
                            || ifry - 1 < 0 || ifry + 2 > refinementData[i].ext[1]
                            || ifrz - 1 < 0 || ifrz + 2 > refinementData[i].ext[2]) {
                        logger.warning("atom: " + a.toString() + " is outside the scalar field box! Ignoring...");
                        continue;
                    }
                }

                // fill in scalar 4x4 array for interpolation
                int uii, vii, wii;
                for (int ui = ifrx - 1; ui < ifrx + 3; ui++) {
                    uii = ui - (ifrx - 1);
                    int pui = Crystal.mod(ui, refinementData[i].ext[0]);
                    for (int vi = ifry - 1; vi < ifry + 3; vi++) {
                        vii = vi - (ifry - 1);
                        int pvi = Crystal.mod(vi, refinementData[i].ext[1]);
                        for (int wi = ifrz - 1; wi < ifrz + 3; wi++) {
                            wii = wi - (ifrz - 1);
                            int pwi = Crystal.mod(wi, refinementData[i].ext[2]);

                            if (refinementData[i].periodic) {
                                scalar[uii][vii][wii] = refinementData[i].getDataIndex(pui, pvi, pwi);
                            } else {
                                scalar[uii][vii][wii] = refinementData[i].getDataIndex(ui, vi, wi);
                            }
                        }
                    }
                }

                // scale and interpolate
                double scale;
                double scaledUdL;
                if (realSpaceFile == null) {
                    scale = -1.0 * lambdai * a.getAtomType().atomicWeight;
                    scaledUdL = -1.0 * dUdL * a.getAtomType().atomicWeight;
                } else {
                    scale = -1.0 * lambdai * realSpaceFile[i].weight * a.getAtomType().atomicWeight;
                    scaledUdL = -1.0 * dUdL * realSpaceFile[i].weight * a.getAtomType().atomicWeight;
                }
                double val = spline.spline(dfrx, dfry, dfrz, scalar, grad);
                sum += scale * val;
                realSpacedUdL += scaledUdL * val;

                if (a.isActive()) {
                    grad[0] = grad[0] * refinementData[i].ni[0];
                    grad[1] = grad[1] * refinementData[i].ni[1];
                    grad[2] = grad[2] * refinementData[i].ni[2];
                    // transpose of toFractional
                    xyz[0] = grad[0] * crystal[i].A00 + grad[1] * crystal[i].A01 + grad[2] * crystal[i].A02;
                    xyz[1] = grad[0] * crystal[i].A10 + grad[1] * crystal[i].A11 + grad[2] * crystal[i].A12;
                    xyz[2] = grad[0] * crystal[i].A20 + grad[1] * crystal[i].A21 + grad[2] * crystal[i].A22;

                    int gradIndex = index * 3;
                    int gx = gradIndex;
                    int gy = gradIndex + 1;
                    int gz = gradIndex + 2;
                    realSpaceGradient[gx] += scale * xyz[0];
                    realSpaceGradient[gy] += scale * xyz[1];
                    realSpaceGradient[gz] += scale * xyz[2];
                    realSpacedUdXdL[gx] += scaledUdL * xyz[0];
                    realSpacedUdXdL[gy] += scaledUdL * xyz[1];
                    realSpacedUdXdL[gz] += scaledUdL * xyz[2];
                    index++;
                }

            }
            refinementData[i].densityscore = sum;
        }

        for (int i = 0; i < n; i++) {
            realSpaceEnergy += refinementData[i].densityscore;
        }

        int index = 0;
        for (Atom a : refinementModel.activeAtoms) {
            int gradIndex = index * 3;
            int gx = gradIndex;
            int gy = gradIndex + 1;
            int gz = gradIndex + 2;
            a.addToXYZGradient(realSpaceGradient[gx], realSpaceGradient[gy], realSpaceGradient[gz]);
            index++;
        }

        return realSpaceEnergy;
    }

    public double getRealSpaceEnergy() {
        return realSpaceEnergy;
    }

    public double getdEdL() {
        return realSpacedUdL;
    }

    public double[] getRealSpaceGradient(double gradient[]) {
        int nAtoms = refinementModel.activeAtoms.length;
        if (gradient == null || gradient.length < nAtoms * 3) {
            gradient = new double[nAtoms * 3];
        }
        for (int i = 0; i < nAtoms * 3; i++) {
            gradient[i] += realSpaceGradient[i];
        }
        return gradient;
    }

    public double[] getdEdXdL(double gradient[]) {
        int nAtoms = refinementModel.activeAtoms.length;
        if (gradient == null || gradient.length < nAtoms * 3) {
            gradient = new double[nAtoms * 3];
        }
        for (int i = 0; i < nAtoms * 3; i++) {
            gradient[i] += realSpacedUdXdL[i];
        }
        return gradient;
    }

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        this.lambda = lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Atom[] getAtomArray() {
        return refinementModel.usedAtoms;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Atom[] getActiveAtomArray() {
        return refinementModel.activeAtoms;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<ArrayList<Residue>> getAltResidues() {
        return refinementModel.altResidues;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<ArrayList<Molecule>> getAltMolecules() {
        return refinementModel.altMolecules;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public MolecularAssembly[] getMolecularAssemblies() {
        return molecularAssemblies;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public RefinementModel getRefinementModel() {
        return refinementModel;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getWeight() {
        return xweight;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setWeight(double weight) {
        this.xweight = weight;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printOptimizationHeader() {
        return "Density score";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printOptimizationUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("%6.2f ", refinementData[i].densityscore));
        }
        return sb.toString();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String printEnergyUpdate() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(String.format("     dataset %d (weight: %5.1f): chemical energy: %8.2f density score: %8.2f\n",
                    i + 1,
                    realSpaceFile[i].weight,
                    molecularAssemblies[0].getPotentialEnergy().getTotalEnergy(),
                    realSpaceFile[i].weight * refinementData[i].densityscore));
        }
        return sb.toString();
    }
}
