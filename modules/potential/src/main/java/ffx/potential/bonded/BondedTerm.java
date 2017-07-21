/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.logging.Logger;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

import edu.rit.pj.reduction.SharedDouble;

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
@SuppressWarnings({"serial", "CloneableImplementsClone"})
public abstract class BondedTerm extends MSNode implements BondedEnergy, Comparable<BondedTerm> {

    private static final Logger logger = Logger.getLogger(BondedTerm.class.getName());
    /**
     * This method sets the Term's id and key by concatenating the respective id
     * and keys of the Atoms that are used in forming the term. Order can be
     * reversed for help in assigning force field parameters for the Term.
     */
    private static StringBuilder idtemp = new StringBuilder();

    protected String id;
    protected Atom atoms[]; // Atoms that are used to form this term
    protected Bond bonds[]; // Bonds that are used to form this term
    protected double value; // Value of the term
    protected double energy; // Energy of the term
    
    protected boolean esvTerm = false;
    /**
     * Lambda value of attached ESV, if present.
     */
    protected double esvLambda = 1.0;
    /**
     * d(lambda_switching_function)/dL.
     * For the linear switch E=L*E1+(1-L)*E0, chain is +/- unity.
     * For other switches E=S(L)*E1+(1-S(L))*E0, set chain to dSdL.
     */
    protected double dedesvChain = 1.0; // Handles sign flip d.t. d[(1-L)*E]/dL
    /**
     * Target for extended variable derivatives.
     */
    protected double esvDerivLocal = 0.0;
    /**
     * Reference to the ESV derivative reduction variable.
     */
    private SharedDouble esvDerivShared = null;
    /**
     * If set, derivative components are filed by source type.
     */
    private HashMap<Class<? extends BondedTerm>,SharedDouble> decompositionMap = null;
    private boolean decomposeEsvDeriv = false;

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
    
    @Override
    public int compareTo(BondedTerm t) {
        return Objects.compare(this, t, bondedComparator);
    }
    public static BondedComparator bondedComparator = new BondedComparator();
    public static class BondedComparator implements Comparator<BondedTerm> {
        private BondedComparator() {}   // singleton
        private static final List<Class<? extends BondedTerm>> naturalOrder =
                new ArrayList<Class<? extends BondedTerm>>() {{
                        add(Bond.class);
                        add(Angle.class);
                        add(StretchBend.class);
                        add(OutOfPlaneBend.class);
                        add(Torsion.class);
                        add(PiOrbitalTorsion.class);
                }};
        /**
         * Sort using position in the naturalOrder list; fallback to alphabetical.
         */
        @Override
        public int compare(BondedTerm me, BondedTerm other) {
            final Class<? extends BondedTerm> oc = other.getClass();
            final Class<? extends BondedTerm> myc = me.getClass();
            int myidx = naturalOrder.indexOf(me.getClass());
            int uridx = naturalOrder.indexOf(other.getClass());
            if (myidx >= 0 && uridx >= 0) {
                return Integer.compare(myidx, uridx);
            } else {
                return String.CASE_INSENSITIVE_ORDER.compare(myc.toString(), oc.toString());
            }
        }
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
    
    public boolean isExtendedSystemMember() {
        return esvTerm;
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
     * Get the Term's value.
     *
     * @return a double.
     */
    public double getValue() {
        return value;
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
     * {@inheritDoc}
     */
    @Override
    public final int hashCode() {
        return hash(SEED, getID().hashCode());
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
        if (atoms == null) {
			return;
        }
        // Reuse the string buffers
        idtemp.delete(0, idtemp.length());
        for (int i = 0; i < atoms.length; i++) {
            Atom a = (reverse) ? atoms[atoms.length - i] : atoms[i];
            if (i != 0) {
                idtemp.append("  ");
            }
            idtemp.append(a.describe(Atom.Descriptions.XyzIndex_Name));
        }
        id = idtemp.toString().intern();
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
     * Under a linear switching function, E=L*E1+(1-L)*E0, chainRule is +1 or -1
     * for lambda and (1-lambda) terms, respectively. Other switches should 
     * set this to d(switch)/d(lambda) as well.
     */
    public void attachExtendedVariable(double lambda, double chainRule, 
            SharedDouble esvBondedDeriv,
            HashMap<Class<? extends BondedTerm>,SharedDouble> decomposition) {
        esvTerm = true;
        esvLambda = lambda;
        dedesvChain = chainRule;
        esvDerivShared = esvBondedDeriv;
        esvDerivLocal = 0.0;        
        if (decomposition != null) {
            decompositionMap = decomposition;
            decomposeEsvDeriv = true;
        } else {
            decompositionMap = null;
            decomposeEsvDeriv = false;
        }
    }
    public void attachExtendedVariable(double lambda, double chainRule,
            SharedDouble esvBondedDeriv) {
        attachExtendedVariable(lambda, chainRule, esvBondedDeriv, null);
    }
    
    public void detachExtendedVariable() {
        esvTerm = false;
        esvLambda = 1.0;
        dedesvChain = 0.0;
        esvDerivLocal = 0.0;
        esvDerivShared = null;
        decompositionMap = null;
        decomposeEsvDeriv = false;
    }
    
    /**
     * Derivative with respect to attached ExtendedVariable lambda, if any.
     * Double.isFinite() check protects against dEdEsv=(energy*chain/lambda) for lambda=0.0
     */
    protected final void setEsvDeriv(double dEdEsv) {
        if (esvTerm) {
            esvDerivLocal = dEdEsv;
        }
    }

    public void reduceEsvDeriv() {
        if (esvTerm) {
//            logf(" :: %.6f from %s", esvDerivLocal, this.toString());
            esvDerivShared.addAndGet(esvDerivLocal);
            if (decomposeEsvDeriv) {
                Class<? extends BondedTerm> source = this.getClass();
                SharedDouble dub = decompositionMap.get(source);
                if (dub == null) {
                    decompositionMap.put(source, new SharedDouble(esvDerivLocal));
                } else {
                    dub.addAndGet(esvDerivLocal);
                }
            }
            esvDerivLocal = 0.0;
        }
    }
}
