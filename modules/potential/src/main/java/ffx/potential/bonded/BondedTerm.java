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
package ffx.potential.bonded;

import java.util.List;
import java.util.logging.Logger;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import ffx.potential.bonded.Atom.Resolution;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The BondedTerm class is extended by all Valence Geometry classes (bond,
 * angle, dihedral, torsion, etc.).
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public abstract class BondedTerm extends MSNode {

    private static final Logger logger = Logger.getLogger(BondedTerm.class.getName());
    /**
     * This method sets the Term's id and key by concatenating the respective id
     * and keys of the Atoms that are used in forming the term. Order can be
     * reversed for help in assigning force field parameters for the Term.
     */
    private static StringBuilder idtemp = new StringBuilder();

    protected String id;
    public Atom atoms[]; // Atoms that are used to form this term
    public Bond bonds[]; // Bonds that are used to form this term
    protected double value; // Value of the term
    protected double energy; // Energy of the term
    protected double esvLambda = 1.0;       // Lambda value of ESV, if present
    /**
     * d(lambda_switching_function)/dL.
     * For the linear switch E=L*E1+(1-L)*E0, chain is +/- unity.
     * For other switches E=S(L)*E1+(1-S(L))*E0, set chain to dSdL.
     */
    protected double dedesvChain = 1.0; // Handles sign flip d.t. d[(1-L)*E]/dL

    /**
     * Default Constructor
     */
    public BondedTerm() {
        super("", 1);
        setAllowsChildren(false);
    }

    /**
     * Constructor which sets the Term's id.
     *
     * @param i a {@link java.lang.String} object.
     */
    public BondedTerm(String i) {
        this();
        id = i;
    }

    /**
     * Checks if all atoms in this BondedTerm are of the given resolution.
     *
     * @param resolution
     *
     * @return true if all atoms in this term are at the same resolution.
     */
    public boolean isResolution(Resolution resolution) {
        for (Atom atom : atoms) {
            if (atom.getResolution() != resolution) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks if at least one atom in this BondedTerm is of the given
     * resolution.
     *
     * @param resolution
     *
     * @return true if at least one atom in this term is of the specified
     * resolution.
     */
    public boolean containsResolution(Resolution resolution) {
        for (Atom atom : atoms) {
            if (atom.getResolution() == resolution) {
                return true;
            }
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        super.destroy();
        id = null;
        value = 0;
        return true;
    }

    /**
     * <p>
     * containsHydrogen</p>
     *
     * @return a boolean.
     */
    public boolean containsHydrogen() {
        for (Atom atom : atoms) {
            if (atom.isHydrogen()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Check if any atom of this BondedTerm has the Lambda flag set.
     *
     * @return True if Lambda is applied to one of the BondedTerm atoms.
     */
    public boolean applyLambda() {
        for (int i = 0; i < atoms.length; i++) {
            if (atoms[i].applyLambda()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Check if all atoms of this BondedTerm have the Lambda flag set.
     *
     * @return True if Lambda is applied to all of the BondedTerm atoms.
     */
    public boolean applyAllLambda() {
        for (int i = 0; i < atoms.length; i++) {
            if (!atoms[i].applyLambda()) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * Overridden method that returns true if object is equals to this, is of
     * the same Class and has the same id.
     */
    @Override
    public final boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        BondedTerm other = (BondedTerm) object;
        return getID().equals(other.getID());
    }

    /**
     * Get the constituent Atom specified by index.
     *
     * @param index a int.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getAtom(int index) {
        if (index >= 0 && index < atoms.length) {
            return atoms[index];
        }
        return null;
    }

    /**
     * Returns all of the Atoms contained in this BondedTerm, regardless of
     * whether they are child nodes in the tree structure. Returns a new array,
     * not a reference to the original array.
     *
     * @return Atoms in this BondedTerm
     */
    public Atom[] getAtomArray() {
        return getAtomArray(true);
    }

    /**
     * Returns all of the Atoms contained in this BondedTerm, regardless of
     * whether they are child nodes in the tree structure.
     *
     * @param returnCopy If true, return a new copy of the Atom array.
     * @return Atoms in this BondedTerm
     */
    public Atom[] getAtomArray(boolean returnCopy) {
        if (returnCopy) {
            int nAtoms = atoms.length;
            Atom[] retAtoms = new Atom[nAtoms];
            System.arraycopy(atoms, 0, retAtoms, 0, nAtoms);
            return retAtoms;
        } else {
            return atoms;
        }
    }

    /**
     * <p>
     * containsAtom</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @return a boolean.
     */
    public boolean containsAtom(Atom atom) {
        for (Atom a : atoms) {
            if (a.equals(atom)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Get the constituent Bond specified by index.
     *
     * @param index a int.
     * @return a {@link ffx.potential.bonded.Bond} object.
     */
    public Bond getBond(int index) {
        if (index >= 0 && index < atoms.length) {
            return bonds[index];
        }
        return null;
    }

    /**
     * Get the Term's id.
     *
     * @return a {@link java.lang.String} object.
     */
    public String getID() {
        return id;
    }

    /**
     * Get the Term's value.
     *
     * @return a double.
     */
    public double getValue() {
        return value;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final int hashCode() {
        return hash(SEED, getID().hashCode());
    }

    /**
     * {@inheritDoc}
     *
     * Prints the toString method to stdout
     */
    @Override
    public void print() {
        logger.info(toString());
    }

    /**
     * Add a constituent Atom to the Term.
     *
     * @param a an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public void setAtoms(Atom a[]) {
        atoms = a;
    }

    /**
     * Add constituent Bonds to the Term.
     *
     * @param b an array of {@link ffx.potential.bonded.Bond} objects.
     */
    public void setBonds(Bond b[]) {
        bonds = b;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
            Material mat) {
        if (atoms == null) {
            return;
        }
        for (Atom atom : atoms) {
            atom.setColor(newColorModel, color, mat);
        }
    }

    /**
     * Sets the Term's id.
     *
     * @param i a {@link java.lang.String} object.
     */
    public void setID(String i) {
        id = i;
    }

    /**
     * <p>
     * setID_Key</p>
     *
     * @param reverse a boolean.
     */
    public final void setID_Key(boolean reverse) {
        Atom a;
        // Reuse the string buffers
        if (idtemp.length() > 0) {
            idtemp.delete(0, idtemp.length());
        }
        if (atoms != null) {
            if (!reverse) {
                for (int i = 0; i < atoms.length; i++) {
                    a = atoms[i];
                    if (i != 0) {
                        idtemp.append("  ").append(a.toShortString());
                    } else {
                        idtemp.append(a.toShortString());
                    }
                }
            } else { // Reverse Order
                for (int i = 0; i < atoms.length; i++) {
                    a = atoms[i];
                    if (i != 0) {
                        idtemp.append("  ").append(a.toShortString());
                    } else {
                        idtemp.append(a.toShortString());
                    }
                }
            }
            id = idtemp.substring(0).toString().intern();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSelected(boolean b) {
        super.setSelected(b);
        if (atoms == null) {
            return;
        }
        for (Atom a : atoms) {
            a.setSelected(b);
        }
        if (!(this instanceof Bond)) {
            if (bonds == null) {
                return;
            }
            for (Bond bond : bonds) {
                bond.setSelected(b);
            }
        }
    }

    /**
     * Sets the Term's value.
     *
     * @param v a double.
     */
    public void setValue(double v) {
        value = v;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel newViewModel,
            List<BranchGroup> newShapes) {
        if (atoms == null) {
            return;
        }
        for (Atom atom : atoms) {
            atom.setView(newViewModel, newShapes);
        }
        if (bonds == null) {
            return;
        }
        for (Bond bond : bonds) {
            bond.setView(newViewModel, newShapes);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.2f,%7.2f)", id, value, energy);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract void update();
    
    /**
     * Under a linear switching function, E=L*E1+(1-L)*E0, chainRule is +1 or -1
     * for lambda and (1-lambda) terms, respectively. Other switches should 
     * set this to d(switch)/d(lambda) as well.
     */
    public void setEsvLambda(double lambda, double chainRule) {
        esvLambda = lambda;
        dedesvChain = chainRule;
    }
    
    /**
     * Derivative with respect to attached ExtendedVariable lambda, if any.
     */
    public final double getdEdEsv() {
        return dedesvChain * energy / esvLambda;
    }
}
