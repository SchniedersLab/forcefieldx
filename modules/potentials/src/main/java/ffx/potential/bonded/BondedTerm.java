/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
package ffx.potential.bonded;

import java.util.List;
import java.util.logging.Logger;

import javax.media.j3d.BranchGroup;
import javax.media.j3d.Material;
import javax.vecmath.Color3f;

/**
 * The BondedTerm class is extended by all Valence Geometry classes (bond,
 * angle, dihedral, torsion, etc.).
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class BondedTerm extends MSNode {

    private static final Logger logger = Logger.getLogger(BondedTerm.class.getName());
    /**
     * This method sets the Term's id and key by concatenating the respective id
     * and keys of the Atoms that are used in forming the term. Order can be
     * reversed for help in finding corresponding Molecular Mechanics file
     * entries for the Term.
     */
    private static StringBuffer idtemp = new StringBuffer();
    protected String id;
    public Atom atoms[]; // Atoms that are used to form this term
    public Bond bonds[]; // Bonds that are used to form this term
    protected double value; // Value of the term
    protected double energy; // Energy of the term

    /**
     * Default Constructor
     */
    public BondedTerm() {
        super("", 1);
        setAllowsChildren(false);
    }

    /**
     * Constructor which sets the Term's id.
     */
    public BondedTerm(String i) {
        this();
        id = i;
    }

    @Override
    public boolean destroy() {
        super.destroy();
        id = null;
        value = 0;
        return true;
    }

    /**
     * Overidden method that returns true if object is equals to this, is of
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
     */
    public Atom getAtom(int index) {
        if (index >= 0 && index < atoms.length) {
            return atoms[index];
        }
        return null;
    }

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
     */
    public Bond getBond(int index) {
        if (index >= 0 && index < atoms.length) {
            return bonds[index];
        }
        return null;
    }

    /**
     * Get the Term's id.
     */
    public String getID() {
        return new String(id);
    }

    /**
     * Get the Term's value.
     */
    public double getValue() {
        return value;
    }

    @Override
    public final int hashCode() {
        return HashCodeUtil.hash(HashCodeUtil.BONDTERMSEED, getID().hashCode());
    }

    /**
     * Prints the toString method to stdout
     */
    @Override
    public void print() {
        logger.info(toString());
    }

    /**
     * Add a constituent Atom to the Term.
     */
    public void setAtoms(Atom a[]) {
        atoms = a;
    }

    /**
     * Add constituent Bonds to the Term.
     */
    public void setBonds(Bond b[]) {
        bonds = b;
    }

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
     */
    public void setID(String i) {
        id = new String(i);
    }

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
     */
    public void setValue(double v) {
        value = v;
    }

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
     * Overidden toString Method returns the Term's id.
     */
    @Override
    public String toString() {
        return String.format("%s  (%7.2f,%7.2f)", id, value, energy);
    }

    @Override
    public abstract void update();
}
