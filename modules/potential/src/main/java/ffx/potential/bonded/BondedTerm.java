// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
// ******************************************************************************
package ffx.potential.bonded;

import ffx.numerics.Constraint;
import ffx.potential.bonded.Atom.Resolution;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Material;
import org.jogamp.vecmath.Color3f;

import java.io.Serial;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.logging.Logger;

/**
 * The BondedTerm class is extended by all Valence Geometry classes (bond, angle, dihedral, torsion,
 * etc.).
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class BondedTerm extends MSNode implements BondedEnergy, Comparable<BondedTerm> {

  @Serial
  private static final long serialVersionUID = 1L;

  private static final Logger logger = Logger.getLogger(BondedTerm.class.getName());
  /**
   * This method sets the Term's id and key by concatenating the respective id and keys of the Atoms
   * that are used in forming the term. Order can be reversed for help in assigning force field
   * parameters for the Term.
   */
  private static final StringBuilder idTemp = new StringBuilder();
  /**
   * Constant <code>bondedComparator</code>
   */
  private static final BondedComparator bondedComparator = new BondedComparator();
  /**
   * ID of this BondedTerm.
   */
  protected String id;
  /**
   * Atoms that are used to form this term.
   */
  protected Atom[] atoms;
  /**
   * Bonds that are used to form this term.
   */
  protected Bond[] bonds;
  /**
   * Value of the term (e.g. bond length, angle, dihedral angle, etc).
   */
  protected double value;
  /**
   * Energy of the term (kcal/mol).
   */
  protected double energy;
  /**
   * Flag indicating if this term is constrained.
   */
  private boolean isConstrained = false;

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
   * Check if any atom of this BondedTerm has the Lambda flag set.
   *
   * @return True if Lambda is applied to one of the BondedTerm atoms.
   */
  public boolean applyLambda() {
    for (Atom atom : atoms) {
      if (atom.applyLambda()) {
        return true;
      }
    }
    return false;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int compareTo(BondedTerm t) {
    return Objects.compare(this, t, bondedComparator);
  }

  /**
   * containsHydrogen
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
   * This method returns true if any atom is marked as being part of a neural network.
   *
   * @return True if any atom is marked as being part of a neural network.
   */
  public boolean isNeuralNetwork() {
    for (Atom atom : atoms) {
      if (atom.isNeuralNetwork()) {
        return true;
      }
    }
    return false;
  }

  /**
   * This method removes terms from a list that are marked as being part of a neural network.
   *
   * @param list The list to check.
   * @param <T>  The Bonded term class.
   */
  public static <T extends BondedTerm> void removeNeuralNetworkTerms(List<T> list) {
    // Remove terms handled by a neutral network.
    List<T> removeList = new ArrayList<>();
    for (T term : list) {
      if (term.isNeuralNetwork()) {
        removeList.add(term);
      }
    }
    list.removeAll(removeList);
  }

  /**
   * Checks if at least one atom in this BondedTerm is of the given resolution.
   *
   * @param resolution a {@link ffx.potential.bonded.Atom.Resolution} object.
   * @return true if at least one atom in this term is of the specified resolution.
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
   * {@inheritDoc}
   *
   * <p>Overridden method that returns true if object is equals to this, is of the same Class and
   * has the same id.
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
   * @param index The index of the Atom to return.
   * @return a {@link ffx.potential.bonded.Atom} object.
   */
  public Atom getAtom(int index) {
    if (index >= 0 && index < atoms.length) {
      return atoms[index];
    }
    return null;
  }

  /**
   * Returns a reference to the Atoms contained in this BondedTerm, regardless of whether they are child nodes in
   * the tree structure.
   *
   * @return Atoms in this BondedTerm
   */
  public Atom[] getAtoms() {
    return atoms;
  }


  /**
   * Returns all Atoms contained in this BondedTerm, regardless of whether they are child nodes in
   * the tree structure. Returns a new array, not a reference to the original array.
   *
   * @return Atoms in this BondedTerm
   */
  public Atom[] getAtomArray() {
    return getAtomArray(true);
  }

  /**
   * Returns all Atoms contained in this BondedTerm, regardless of whether they are child nodes in
   * the tree structure.
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
   * Get the constituent Bond specified by index.
   *
   * @param index The index of the Bond to return.
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
   * This method returns true if any atom is being used.
   *
   * @return True if any atom is being used.
   */
  public boolean getUse() {
    for (Atom atom : atoms) {
      if (atom.getUse()) {
        return true;
      }
    }
    return false;
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
   * Get the Term's value.
   *
   * @return a double.
   */
  public double getValue() {
    return value;
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
  public int hashCode() {
    return Objects.hash(getID());
  }

  /**
   * Check if this BondedTerm is constrained.
   *
   * @return If constrained.
   */
  public boolean isConstrained() {
    return isConstrained;
  }

  /**
   * Check if this BondedTerm is lambda-sensitive (e.g., a softcore dihedral).
   *
   * @return True if Lambda affects the energy of this term.
   */
  public boolean isLambdaScaled() {
    return false;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Prints the toString method to stdout
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
  public void setAtoms(Atom[] a) {
    atoms = a;
  }

  /**
   * Add constituent Bonds to the Term.
   *
   * @param b an array of {@link ffx.potential.bonded.Bond} objects.
   */
  public void setBonds(Bond[] b) {
    bonds = b;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
    if (atoms == null) {
      return;
    }
    for (Atom atom : atoms) {
      atom.setColor(newColorModel, color, mat);
    }
  }

  /**
   * Sets the Constraint on this bond (clearing it if null). May recursively set the Constraint on
   * component terms (i.e. an Angle will call setConstraint on its component Bonds).
   *
   * @param c Constraint or null to clear.
   */
  public void setConstraint(Constraint c) {
    isConstrained = c != null;
  }

  /**
   * setID_Key
   *
   * @param reverse a boolean.
   */
  public final void setID_Key(boolean reverse) {
    if (atoms == null) {
      return;
    }
    // Reuse the string buffers
    idTemp.delete(0, idTemp.length());
    for (int i = 0; i < atoms.length; i++) {
      Atom a = (reverse) ? atoms[atoms.length - 1 - i] : atoms[i];
      if (i != 0) {
        idTemp.append("  ");
      }
      idTemp.append(a.describe(Atom.Descriptions.XyzIndex_Name));
    }
    id = idTemp.toString().intern();
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
   * {@inheritDoc}
   */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
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
   * <p>Overridden toString Method returns the Term's id.
   */
  @Override
  public String toString() {
    return String.format("%s  (%7.2f,%7.2f)", id, value, energy);
  }

  /**
   * Check if all atoms of this BondedTerm have the Lambda flag set.
   *
   * @return True if Lambda is applied to all BondedTerm atoms.
   */
  boolean applyAllLambda() {
    for (Atom atom : atoms) {
      if (!atom.applyLambda()) {
        return false;
      }
    }
    return true;
  }

  /**
   * containsAtom
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @return a boolean.
   */
  boolean containsAtom(Atom atom) {
    for (Atom a : atoms) {
      if (a.equals(atom)) {
        return true;
      }
    }
    return false;
  }

  public static class BondedComparator implements Comparator<BondedTerm> {

    private static final List<Class<? extends BondedTerm>> naturalOrder =
        new ArrayList<>() {
          {
            add(Bond.class);
            add(Angle.class);
            add(StretchBend.class);
            add(UreyBradley.class);
            add(OutOfPlaneBend.class);
            add(Torsion.class);
            add(ImproperTorsion.class);
            add(PiOrbitalTorsion.class);
            add(StretchTorsion.class);
            add(AngleTorsion.class);
            add(TorsionTorsion.class);
          }
        };

    private BondedComparator() {
    } // singleton

    /**
     * Sort using position in the naturalOrder list; fallback to alphabetical.
     */
    @Override
    public int compare(BondedTerm bondedTerm1, BondedTerm bondedTerm2) {
      final Class<? extends BondedTerm> class1 = bondedTerm2.getClass();
      final Class<? extends BondedTerm> class2 = bondedTerm1.getClass();
      int order1 = naturalOrder.indexOf(bondedTerm1.getClass());
      int order2 = naturalOrder.indexOf(bondedTerm2.getClass());
      if (order1 >= 0 && order2 >= 0) {
        return Integer.compare(order1, order2);
      } else {
        return String.CASE_INSENSITIVE_ORDER.compare(class2.toString(), class1.toString());
      }
    }
  }
}
