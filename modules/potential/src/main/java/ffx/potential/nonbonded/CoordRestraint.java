/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential.nonbonded;

import java.util.logging.Logger;

import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.pow;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;

import static ffx.numerics.VectorMath.rsq;

/**
 * Restrain atoms to their initial coordinates.
 *
 * @author Michael J. Schnieders
 */
public class CoordRestraint implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(CoordRestraint.class.getName());

    private final Atom atoms[];
    private final double initialCoordinates[][];
    private int nAtoms = 0;
    /**
     * Force constant in Kcal/mole/Angstrom.
     */
    private final double forceConstant;
    private final double a1[] = new double[3];
    private final double dx[] = new double[3];
    private double lambda = 1.0;
    private final double lambdaExp = 1.0;
    private double lambdaPow = pow(lambda, lambdaExp);
    private double dLambdaPow = lambdaExp * pow(lambda, lambdaExp - 1.0);
    private double d2LambdaPow = 0;
    private final double lambdaWindow = 1.0;
    private double dEdL = 0.0;
    private double d2EdL2 = 0.0;
    private final double lambdaGradient[];
    private boolean lambdaTerm = false;
    private int atom1Index;
    private int atom2Index;
    private int atom3Index;
    private AtomType atom1Type;
    private AtomType atom2Type;
    private AtomType atom3Type;

//    private AtomType definePlane1 = null;
//    private AtomType definePlane2 = null;
//    private AtomType definePlane3 = null;
    /**
     * This CoordRestraint is based on the unit cell parameters and symmetry
     * operators of the supplied crystal.
     *
     * @param atoms the Atom array to base this CoordRestraint on.
     * @param forceField the ForceField to apply.
     */
    public CoordRestraint(Atom[] atoms, ForceField forceField) {
        this.atoms = atoms;
        nAtoms = atoms.length;

        //lambdaTerm = false;
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);

        if (lambdaTerm) {
            lambdaGradient = new double[nAtoms * 3];
        } else {
            lambdaGradient = null;
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }

        forceConstant = forceField.getDouble(ForceField.ForceFieldDouble.RESTRAINT_K, 10.0);

        logger.info("\n Coordinate Restraint:");

        initialCoordinates = new double[3][nAtoms];
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            initialCoordinates[0][i] = a.getX();
            initialCoordinates[1][i] = a.getY();
            initialCoordinates[2][i] = a.getZ();
            a.print();
        }
//        if (System.getProperty("atomrestraint1") != null) {
//            setPlaneAtoms();
//        }
    }

    public double residual(boolean gradient, boolean print) {
        //pdbAtomRestraints uses atom indexes to restrain specific atoms
        String atomIndexRestraints = System.getProperty("pdbAtomRestraints");
        if (atomIndexRestraints != null) {
            String tokens[] = atomIndexRestraints.split(",");
        //Pick up atom index for reference when looking at muliple molecules
            atom1Index = Integer.parseInt(tokens[0]);
            atom2Index = Integer.parseInt(tokens[1]);
            atom3Index = Integer.parseInt(tokens[2]);
        }
        
        //xyzAtomRestraints uses atom types to restrain specific atoms. This can result in more
        //atoms being restrained than desired since atom types are not unique to each atom.
        String atomTypeRestraints = System.getProperty("xyzAtomRestraints");
        if (atomTypeRestraints != null) {
            String tokens[] = atomTypeRestraints.split(",");
        //Pick up atom type for reference when looking at muliple molecules
            atom1Type = atoms[Integer.parseInt(tokens[0])].getAtomType();
            atom2Type = atoms[Integer.parseInt(tokens[1])].getAtomType();
            atom3Type = atoms[Integer.parseInt(tokens[2])].getAtomType();
        }
        if (lambdaTerm) {
            dEdL = 0.0;
            d2EdL2 = 0.0;
            fill(lambdaGradient, 0.0);
        }
        int atomIndex = 0;
        //Assuming that the first molecule in pdb is labeled as 1.
        int moleculeNumber = 1;
        double residual = 0.0;
        double fx2 = forceConstant * 2.0;
        for (int i = 0; i < nAtoms; i++) {
            // Current atomic coordinates.
            Atom atom = atoms[i];
            logger.info("Number of atoms in atoms[] = " + atom.getResidueNumber());
            if (atomTypeRestraints != null) {
                if (atom.getAtomType() == atom1Type || atom.getAtomType() == atom2Type || atom.getAtomType() == atom3Type) {
                    atom.getXYZ(a1);
                    // Compute their separation.
                    dx[0] = a1[0] - initialCoordinates[0][i];
                    dx[1] = a1[1] - initialCoordinates[1][i];
                    dx[2] = a1[2] - initialCoordinates[2][i];
                    logger.info(String.format("restrained atom = %d", atom.getXYZIndex()));
                }
                else {
                    continue;
                }
                double r2 = rsq(dx);
                // Apply the minimum image convention.
                // double r2 = crystal.image(dx);
                //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
                residual += r2;
                if (gradient || lambdaTerm) {
                    final double dedx = dx[0] * fx2;
                    final double dedy = dx[1] * fx2;
                    final double dedz = dx[2] * fx2;
                    if (gradient) {
                        atom.addToXYZGradient(lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
                    }
                    if (lambdaTerm) {
                        int j3 = i * 3;
                        lambdaGradient[j3] = dLambdaPow * dedx;
                        lambdaGradient[j3 + 1] = dLambdaPow * dedy;
                        lambdaGradient[j3 + 2] = dLambdaPow * dedz;
                    }
                }
            } else if (atomIndexRestraints != null) {
                atomIndex += 1;
                // Following if statement works only for .pdb files since they give residueNumber
                if (atom.getResidueNumber() != moleculeNumber) {
                    atomIndex = 1;
                    moleculeNumber = atom.getResidueNumber();
                }
                if (atomIndex == atom1Index || atomIndex == atom2Index || atomIndex == atom3Index) {
                    atom.getXYZ(a1);
                    // Compute their separation.
                    dx[0] = a1[0] - initialCoordinates[0][i];
                    dx[1] = a1[1] - initialCoordinates[1][i];
                    dx[2] = a1[2] - initialCoordinates[2][i];
                    logger.info(String.format("restrained atom = %d", atom.getXYZIndex()));
                }
                else {
                    continue;
                }
                double r2 = rsq(dx);
                // Apply the minimum image convention.
                // double r2 = crystal.image(dx);
                //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
                residual += r2;
                if (gradient || lambdaTerm) {
                    final double dedx = dx[0] * fx2;
                    final double dedy = dx[1] * fx2;
                    final double dedz = dx[2] * fx2;
                    if (gradient) {
                        atom.addToXYZGradient(lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
                    }
                    if (lambdaTerm) {
                        int j3 = i * 3;
                        lambdaGradient[j3] = dLambdaPow * dedx;
                        lambdaGradient[j3 + 1] = dLambdaPow * dedy;
                        lambdaGradient[j3 + 2] = dLambdaPow * dedz;
                    }
                }
            } else {

                if (atom.isHydrogen()) {
                    continue;
                }

                atom.getXYZ(a1);
                // Compute their separation.
                dx[0] = a1[0] - initialCoordinates[0][i];
                dx[1] = a1[1] - initialCoordinates[1][i];
                dx[2] = a1[2] - initialCoordinates[2][i];
                double r2 = rsq(dx);
                // Apply the minimum image convention.
                // double r2 = crystal.image(dx);
                //logger.info(String.format(" %d %16.8f", j, Math.sqrt(r2)));
                residual += r2;
                if (gradient || lambdaTerm) {
                    final double dedx = dx[0] * fx2;
                    final double dedy = dx[1] * fx2;
                    final double dedz = dx[2] * fx2;
                    if (gradient) {
                        atom.addToXYZGradient(lambdaPow * dedx, lambdaPow * dedy, lambdaPow * dedz);
                    }
                    if (lambdaTerm) {
                        int j3 = i * 3;
                        lambdaGradient[j3] = dLambdaPow * dedx;
                        lambdaGradient[j3 + 1] = dLambdaPow * dedy;
                        lambdaGradient[j3 + 2] = dLambdaPow * dedz;
                    }
                }
            }
        }
        if (lambdaTerm) {
            dEdL = dLambdaPow * forceConstant * residual;
            d2EdL2 = d2LambdaPow * forceConstant * residual;
        }

        //logger.info(String.format(" Restraint Energy %16.8f", forceConstant * residual * lambdaPow));
        return forceConstant * residual * lambdaPow;
    }

//    public final void setPlaneAtoms() {
//        int planeAtom1 = Integer.parseInt(System.getProperty("atomrestraint1"));
//        int planeAtom2 = Integer.parseInt(System.getProperty("atomrestraint2"));
//        int planeAtom3 = Integer.parseInt(System.getProperty("atomrestraint3"));
//        for (int i = 0; i < nAtoms; i++) {
//            // Current atomic coordinates.
//            Atom atom = atoms[i];
//            if (atom.getXYZIndex() == planeAtom1) {
//                definePlane1 = atom.getAtomType();
//            } else if (atom.getXYZIndex() == planeAtom2) {
//                definePlane2 = atom.getAtomType();
//                
//                //Should the following else if statetment be == planeAtom3??
//
//            } else if (atom.getXYZIndex() == planeAtom3) {
//                definePlane3 = atom.getAtomType();
//            }
//        }
//    }
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;
            if (this.lambda <= lambdaWindow) {
                double dldgl = 1.0 / lambdaWindow;
                double l = dldgl * this.lambda;
                double l2 = l * l;
                double l3 = l2 * l;
                double l4 = l2 * l2;
                double l5 = l4 * l;
                double c3 = 10.0;
                double c4 = -15.0;
                double c5 = 6.0;
                double threeC3 = 3.0 * c3;
                double sixC3 = 6.0 * c3;
                double fourC4 = 4.0 * c4;
                double twelveC4 = 12.0 * c4;
                double fiveC5 = 5.0 * c5;
                double twentyC5 = 20.0 * c5;
                lambdaPow = c3 * l3 + c4 * l4 + c5 * l5;
                dLambdaPow = (threeC3 * l2 + fourC4 * l3 + fiveC5 * l4) * dldgl;
                d2LambdaPow = (sixC3 * l + twelveC4 * l2 + twentyC5 * l3) * dldgl * dldgl;
            } else {
                lambdaPow = 1.0;
                dLambdaPow = 0.0;
                d2LambdaPow = 0.0;
            }
        } else {
            this.lambda = 1.0;
            lambdaPow = 1.0;
            dLambdaPow = 0.0;
            d2LambdaPow = 0.0;
        }
    }

    public void setLambdaTerm(boolean lambdaTerm) {
        this.lambdaTerm = lambdaTerm;
        setLambda(lambda);
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getdEdL() {
        if (lambdaTerm) {
            return dEdL;
        } else {
            return 0.0;
        }
    }

    @Override
    public double getd2EdL2() {
        if (lambdaTerm) {
            return d2EdL2;
        } else {
            return 0.0;
        }
    }

    @Override
    public void getdEdXdL(double[] gradient) {
        if (lambdaTerm) {
            int n3 = nAtoms * 3;
            for (int i = 0; i < n3; i++) {
                gradient[i] += lambdaGradient[i];
            }
        }
    }
}
