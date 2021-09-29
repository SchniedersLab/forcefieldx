// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

import static ffx.potential.bonded.AminoAcidUtils.AA1toAA3;
import static ffx.potential.bonded.NucleicAcidUtils.NA1toNA3;
import static java.lang.System.arraycopy;

import ffx.potential.bonded.AminoAcidUtils.AminoAcid1;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.NucleicAcidUtils.NucleicAcid1;
import ffx.potential.bonded.NucleicAcidUtils.NucleicAcid3;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TitrationUtils;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jogamp.java3d.Canvas3D;
import org.jogamp.java3d.J3DGraphics2D;
import org.jogamp.java3d.Material;
import org.jogamp.java3d.Node;
import org.jogamp.vecmath.Color3f;
import org.jogamp.vecmath.Point2d;
import org.jogamp.vecmath.Point3d;

/**
 * The Residue class represents individual amino acids or nucleic acid bases.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class Residue extends MSGroup implements Comparable<Residue> {

  private static final Logger logger = Logger.getLogger(Residue.class.getName());
  /** Compare residues first on seg ID, then residue number, then residue type, then name. */
  private static final Comparator<Residue> resComparator =
      Comparator.comparing(Residue::getSegID)
          .thenComparingInt(Residue::getResidueNumber)
          .thenComparing(Residue::getResidueType)
          .thenComparing(Residue::getName);

  /** Constant <code>NA3Color</code> */
  private static final HashMap<NucleicAcid3, Color3f> NA3Color = new HashMap<>();
  /** Constant <code>AA3Color</code> */
  private static final HashMap<AminoAcid3, Color3f> AA3Color = new HashMap<>();
  /** Constant <code>SSTypeColor</code> */
  private static final HashMap<SSType, Color3f> SSTypeColor = new HashMap<>();
  /** Constant <code>Ramachandran="new String[17]"</code> */
  public static String[] Ramachandran = new String[17];

  private static final Point3d point3d = new Point3d();
  private static final Point2d point2d = new Point2d();

  static {
    NA3Color.put(NucleicAcid3.ADE, RendererCache.RED);
    NA3Color.put(NucleicAcid3.CYT, RendererCache.MAGENTA);
    NA3Color.put(NucleicAcid3.GUA, RendererCache.BLUE);
    NA3Color.put(NucleicAcid3.URI, RendererCache.YELLOW);
    NA3Color.put(NucleicAcid3.DAD, RendererCache.RED);
    NA3Color.put(NucleicAcid3.DCY, RendererCache.MAGENTA);
    NA3Color.put(NucleicAcid3.DGU, RendererCache.BLUE);
    NA3Color.put(NucleicAcid3.DTY, RendererCache.ORANGE);
    NA3Color.put(NucleicAcid3.MP1, RendererCache.GREEN);
    NA3Color.put(NucleicAcid3.DP2, RendererCache.GREEN);
    NA3Color.put(NucleicAcid3.TP3, RendererCache.GREEN);
    NA3Color.put(NucleicAcid3.UNK, RendererCache.CYAN);

    AA3Color.put(AminoAcid3.ALA, RendererCache.GRAY);
    AA3Color.put(AminoAcid3.ARG, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.ASN, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.ASP, RendererCache.RED);
    AA3Color.put(AminoAcid3.CYS, RendererCache.YELLOW);
    AA3Color.put(AminoAcid3.GLN, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.GLU, RendererCache.RED);
    AA3Color.put(AminoAcid3.GLH, RendererCache.RED);
    AA3Color.put(AminoAcid3.GLD, RendererCache.RED);
    AA3Color.put(AminoAcid3.GLY, RendererCache.GRAY);
    AA3Color.put(AminoAcid3.ILE, RendererCache.GRAY);
    AA3Color.put(AminoAcid3.LEU, RendererCache.GRAY);
    AA3Color.put(AminoAcid3.LYS, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.MET, RendererCache.YELLOW);
    AA3Color.put(AminoAcid3.PHE, RendererCache.GREEN);
    AA3Color.put(AminoAcid3.PRO, RendererCache.ORANGE);
    AA3Color.put(AminoAcid3.SER, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.THR, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.TRP, RendererCache.GREEN);
    AA3Color.put(AminoAcid3.TYR, RendererCache.GREEN);
    AA3Color.put(AminoAcid3.VAL, RendererCache.GRAY);
    AA3Color.put(AminoAcid3.HIS, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.HIE, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.HID, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.ORN, RendererCache.ORANGE);
    AA3Color.put(AminoAcid3.AIB, RendererCache.ORANGE);
    AA3Color.put(AminoAcid3.PCA, RendererCache.ORANGE);
    AA3Color.put(AminoAcid3.FOR, RendererCache.RED);
    AA3Color.put(AminoAcid3.ACE, RendererCache.RED);
    AA3Color.put(AminoAcid3.NH2, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.NME, RendererCache.BLUE);
    AA3Color.put(AminoAcid3.UNK, RendererCache.MAGENTA);

    SSTypeColor.put(SSType.NONE, RendererCache.WHITE);
    SSTypeColor.put(SSType.SHEET, RendererCache.PINK);
    SSTypeColor.put(SSType.HELIX, RendererCache.BLUE);
    SSTypeColor.put(SSType.TURN, RendererCache.YELLOW);

    Ramachandran[0] = "Default (Extended)       [-135.0  135.0]";
    Ramachandran[1] = "Alpha Helix (R)          [ -57.0  -47.0]";
    Ramachandran[2] = "Alpha Helix (L)          [  57.0   47.0]";
    Ramachandran[3] = "3-10 Helix               [ -49.0  -26.0]";
    Ramachandran[4] = "Pi Helix                 [ -57.0  -70.0]";
    Ramachandran[5] = "Polyproline II Helix     [ -79.0  149.0]";
    Ramachandran[6] = "Parallel Beta Strand     [-119.0  113.0]";
    Ramachandran[7] = "Antiparallel Beta Strand [-139.0  135.0]";
    Ramachandran[8] = "Beta-Hairpin 2' (i+1)    [  90.0 -170.0]";
    Ramachandran[9] = "Beta-Hairpin 2' (i+2)    [ -80.0  -10.0]";
    Ramachandran[10] = "Beta-Hairpin 1' (i+1)    [  57.0   47.0]";
    Ramachandran[11] = "Beta-Hairpin 1' (i+2)    [  57.0   47.0]";
    Ramachandran[12] = "Beta-Hairpin 1  (i+1)    [ -57.0  -47.0]";
    Ramachandran[13] = "Beta-Hairpin 1  (i+2)    [ -57.0  -47.0]";
    Ramachandran[14] = "Beta-Hairpin 1  (i+3)    [  90.0 -170.0]";
    Ramachandran[15] = "Beta-Hairpin 3' (i+1)    [  57.0   47.0]";
    Ramachandran[16] = "Beta-Hairpin 3' (i+2)    [ -80.0  -10.0]";
  }

  /** Residue type. */
  ResidueType residueType;
  /** The residue number of this residue in a chain. */
  private int resNumber;
  /** Possibly redundant PDB chain ID. */
  private Character chainID;
  /** Unique segID. */
  private String segID;
  /** The rotamers for this residue. */
  private Rotamer[] rotamers = null;
  /** The current rotamer in use. */
  private Rotamer currentRotamer = null;
  /** Short string describing this residue. */
  private String shortString = null;
  /** 3-letter amino acid code. */
  private AminoAcid3 aa;
  /** 3-letter nucleic acid code. */
  private NucleicAcid3 na;
  /**
   * If this is set, then ASP, GLU, LYS and HIS Rotamers will be titratable.
   */
  private TitrationUtils titrationUtils = null;

  /**
   * These arrays store default coordinates for certain atoms in nucleic acid Residues. C1', O4', and
   * C4' are the critical sugar atoms off which every other atom is drawn when applyRotamer is
   * called. The backbone corrections, however, move these atoms, so they must be reverted to their
   * original coordinates each time applyRotamer is called.
   *
   * <p>O3' North and South coordinates are technically non-essential, as they could be derived from
   * C1', O4', C4', and a given sugar pucker, however, it is much less computationally expensive to
   * calculate them once and then store them.
   *
   * <p>TODO: Add O3' coordinates for the DNA C3'-exo configuration.
   */
  private double[] O3sNorthCoords, O3sSouthCoords, C1sCoords, O4sCoords, C4sCoords;

  /**
   * Default Constructor where num is this Residue's position in the Polymer.
   *
   * @param num a int.
   * @param rt a {@link ResidueType} object.
   */
  public Residue(int num, ResidueType rt) {
    super();
    resNumber = num;
    residueType = rt;
    assignResidueType();
  }

  /**
   * Constructor for Residue.
   *
   * @param name a {@link java.lang.String} object.
   * @param rt a {@link ResidueType} object.
   */
  public Residue(String name, ResidueType rt) {
    super(name);
    residueType = rt;
    assignResidueType();
  }

  /**
   * Name is the residue's 3 letter abbreviation and num is its position in the Polymer.
   *
   * @param name a {@link java.lang.String} object.
   * @param num a int.
   * @param rt a {@link ResidueType} object.
   */
  public Residue(String name, int num, ResidueType rt) {
    this(name, rt);
    resNumber = num;
  }

  /**
   * Name is the residue's 3 letter abbreviation and num is its position in the Polymer.
   *
   * @param name a {@link java.lang.String} object.
   * @param resNumber a int.
   * @param rt a {@link ResidueType} object.
   * @param chainID a {@link java.lang.Character} object.
   * @param segID a {@link java.lang.String} object.
   */
  public Residue(String name, int resNumber, ResidueType rt, Character chainID, String segID) {
    this(name, rt);
    this.resNumber = resNumber;
    this.chainID = chainID;
    this.segID = segID;
  }

  /**
   * As above, with atoms being a FNode with this Residue's atoms as child nodes
   *
   * @param name a {@link java.lang.String} object.
   * @param num a int.
   * @param atoms a {@link ffx.potential.bonded.MSNode} object.
   * @param rt a {@link ResidueType} object.
   * @param forceField the ForceField to use when created bonded terms.
   */
  public Residue(String name, int num, MSNode atoms, ResidueType rt, ForceField forceField) {
    super(name, atoms);
    resNumber = num;
    residueType = rt;
    assignResidueType();
    finalize(true, forceField);
  }

  @Override
  public void setName(String name) {
    super.setName(name);
    if (residueType != null) {
      assignResidueType();
    }
    for (Atom atom : getAtomList()) {
      atom.setResName(name);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Allows adding Atoms to the Residue.
   */
  @Override
  public MSNode addMSNode(MSNode o) {
    Atom currentAtom = null;
    if (o instanceof Atom) {
      Atom newAtom = (Atom) o;
      Character newAlt = newAtom.getAltLoc();
      MSNode atoms = getAtomNode();
      currentAtom = (Atom) atoms.contains(newAtom);
      if (currentAtom == null) {
        currentAtom = newAtom;
        currentAtom.setResName(getName());
        currentAtom.setResidueNumber(resNumber);
        atoms.add(newAtom);
        setFinalized(false);
      } else {
        // Allow overwriting of the root alternate conformer (' ' or 'A').
        Character currentAlt = currentAtom.getAltLoc();
        if (currentAlt.equals(' ') || currentAlt.equals('A')) {
          if (!newAlt.equals(' ') && !newAlt.equals('A')) {
            newAtom.setXyzIndex(currentAtom.getXyzIndex());
            atoms.remove(currentAtom);
            currentAtom = newAtom;

            currentAtom.setResName(getName());
            currentAtom.setResidueNumber(resNumber);

            atoms.add(currentAtom);
            setFinalized(false);
          }
        }
      }
    } else {
      logger.warning("Only an Atom can be added to a Residue.");
    }
    return currentAtom;
  }

  @Override
  public int compareTo(Residue o) {
    if (this.equals(o)) {
      return 0;
    }
    return resComparator.compare(this, o);
  }

  /** {@inheritDoc} */
  @Override
  public void drawLabel(Canvas3D canvas, J3DGraphics2D g2d, Node node) {
    if (RendererCache.labelResidues) {
      double[] d = getCenter();
      point3d.x = d[0];
      point3d.y = d[1];
      point3d.z = d[2];
      RendererCache.getScreenCoordinate(canvas, node, point3d, point2d);
      g2d.drawString(getName(), (float) point2d.x, (float) point2d.y);
    }
    if (RendererCache.labelAtoms) {
      super.drawLabel(canvas, g2d, node);
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    Residue residue = (Residue) o;
    return Objects.equals(segID, residue.segID)
        && Objects.equals(getResidueNumber(), residue.getResidueNumber())
        && residueType == residue.residueType
        && Objects.equals(getName(), residue.getName());
  }

  /**
   * {@inheritDoc}
   *
   * <p>The Finalize method should be called once all atoms have been added to the Residue. Geometry
   * objects (Bonds, Angles, etc) are then formed, followed by a determination of under-constrained
   * (Dangeling) atoms.
   */
  @Override
  public void finalize(boolean finalizeGeometry, ForceField forceField) {
    setFinalized(false);
    getAtomNode().setName("Atoms (" + getAtomList().size() + ")");
    if (finalizeGeometry) {
      assignBondedTerms(forceField);
      removeLeaves();
    }
    setCenter(getMultiScaleCenter(false));
    setFinalized(true);
  }

  public void setTitrationUtils(TitrationUtils titrationUtils) {
    this.titrationUtils = titrationUtils;
  }

  public TitrationUtils getTitrationUtils() {
    return titrationUtils;
  }

  /**
   * getAminoAcid3.
   *
   * @return a {@link AminoAcid3} object.
   */
  public AminoAcid3 getAminoAcid3() {
    if (this.residueType != ResidueType.AA) {
      throw new IllegalArgumentException(
          String.format(" This residue is " + "not an amino acid: %s", this));
    }
    return aa;
  }

  /**
   * Returns a list of backbone atoms; for our purposes, nucleic acid backbone atoms are those of the
   * nucleobase.
   *
   * @return List of backbone (or nucleobase) atoms.
   */
  public List<Atom> getBackboneAtoms() {
    List<Atom> atoms = getAtomList();
    List<Atom> ret;
    switch (residueType) {
      case NA:
        ret = new ArrayList<>(atoms);
        for (Atom atom : atoms) {
          String name = atom.getName().toUpperCase();
          if (name.contains("'")
              || name.equals("P")
              || name.startsWith("OP")
              || name.equals("H5T")
              || name.equals("H3T")) {
            ret.remove(atom);
          }
        }
        return ret;
      case AA:
        ret = new ArrayList<>();
        tryAddAtom(ret, "N");
        tryAddAtom(ret, "CA");
        tryAddAtom(ret, "C");
        tryAddAtom(ret, "O");
        tryAddAtom(ret, "OXT"); // C-terminal residues
        tryAddAtom(ret, "OT2"); // Probably alternate name for OXT.
        tryAddAtom(ret, "H1"); // N-terminal residues
        tryAddAtom(ret, "H2");
        tryAddAtom(ret, "H3");
        tryAddAtom(ret, "H");
        tryAddAtom(ret, "HA");
        tryAddAtom(ret, "HA2"); // Glycines
        tryAddAtom(ret, "HA3");
        return ret;
      default:
        return new ArrayList<>(1); // Return empty list.
    }
  }

  /**
   * Returns this Residues Parent Polymer name.
   *
   * @return a {@link java.lang.Character} object.
   */
  public Character getChainID() {
    return chainID;
  }

  /**
   * Setter for the field <code>chainID</code>.
   *
   * @param c a {@link java.lang.Character} object.
   */
  public void setChainID(Character c) {
    chainID = c;

    for (Atom atom : getAtomList()) {
      atom.setChainID(c);
    }
  }

  /**
   * Returns the Residue bonded to this Residue at this Residue's 3' or C-terminal end. Any use of
   * this method to add Residues to a sliding window or similar MUST not add that residue if that
   * residue has no Rotamers, as several algorithms (such as the distance matrix) assume that all
   * Residues being optimized have Rotamers.
   *
   * @return The next Residue.
   */
  public Residue getNextResidue() {
    switch (residueType) {
      case AA: {
        Atom carbon = (Atom) getAtomNode("C");
        if (carbon == null) {
          return null;
        }
        List<Bond> bonds = carbon.getBonds();
        for (Bond b : bonds) {
          Atom other = b.get1_2(carbon);
          if (other.getName().equalsIgnoreCase("N")) {
            return (Residue) other.getParent().getParent();
          }
        }
        break;
      }
      case NA:
        Atom oxygen = (Atom) getAtomNode("O3'");
        if (oxygen == null) {
          return null;
        }
        List<Bond> bonds = oxygen.getBonds();
        for (Bond b : bonds) {
          Atom other = b.get1_2(oxygen);
          if (other.getName().equalsIgnoreCase("P")) {
            return (Residue) other.getParent().getParent();
          }
        }
        break;
      default:
        return null;
    }
    return null;
    // Will generally indicate that you passed in a chain-terminal residue.
  }

  /**
   * getNucleicAcid3.
   *
   * @return a {@link NucleicAcid3} object.
   */
  public NucleicAcid3 getNucleicAcid3() {
    if (this.residueType != ResidueType.NA) {
      throw new IllegalArgumentException(
          String.format(" This residue is not a nucleic acid: %s", this));
    }
    return na;
  }

  /**
   * Returns the NucleicAcid3 corresponding to this Residue, with additional robust checking for 1-
   * or 2-letter names.
   *
   * @param matchShortName Try to match 1- or 2-letter names (e.g. A to ADE).
   * @return a {@link NucleicAcid3} object.
   */
  public NucleicAcid3 getNucleicAcid3(boolean matchShortName) {
    NucleicAcid3 na3 = getNucleicAcid3();
    if (na3 == NucleicAcidUtils.NucleicAcid3.UNK && matchShortName) {
      switch (getName()) {
        case "A":
          return NucleicAcidUtils.NucleicAcid3.ADE;
        case "C":
          return NucleicAcidUtils.NucleicAcid3.CYT;
        case "G":
          return NucleicAcidUtils.NucleicAcid3.GUA;
        case "T":
          return NucleicAcidUtils.NucleicAcid3.THY;
        case "U":
          return NucleicAcidUtils.NucleicAcid3.URI;
        case "DA":
          return NucleicAcidUtils.NucleicAcid3.DAD;
        case "DC":
          return NucleicAcidUtils.NucleicAcid3.DCY;
        case "DG":
          return NucleicAcidUtils.NucleicAcid3.DGU;
        case "DT":
          return NucleicAcidUtils.NucleicAcid3.DTY;
        case "DU":
          throw new IllegalArgumentException(
              " No NucleicAcid3 enum exists for DU (presumed to be deoxy-uracil)!");
      }
    }
    return na3;
  }

  /**
   * Returns the Residue bonded to this Residue at this Residue's 5' or N-terminal end. Any use of
   * this method to add Residues to a sliding window or similar MUST not add that residue if that
   * residue has no Rotamers, as several algorithms (such as the distance matrix) assume that all
   * Residues being optimized have Rotamers.
   *
   * @return The previous Residue.
   */
  public Residue getPreviousResidue() {
    switch (residueType) {
      case AA: {
        Atom nitrogen = (Atom) getAtomNode("N");
        if (nitrogen == null) {
          return null;
        }
        List<Bond> bonds = nitrogen.getBonds();
        for (Bond b : bonds) {
          Atom other = b.get1_2(nitrogen);
          if (other.getName().equalsIgnoreCase("C")) {
            return (Residue) other.getParent().getParent();
          }
        }
        break;
      }
      case NA:
        Atom phosphate = (Atom) getAtomNode("P");
        if (phosphate == null) {
          return null;
        }
        List<Bond> bonds = phosphate.getBonds();
        for (Bond b : bonds) {
          Atom other = b.get1_2(phosphate);
          if (other.getName().equalsIgnoreCase("O3'")) {
            return (Residue) other.getParent().getParent();
          }
        }
        break;
      default:
        return null;
    }
    return null;
    // Will generally indicate that you passed in a chain-starting residue.
  }

  /**
   * Returns a reference Atom for a Residue, primarily intended for rough distance calculations. This
   * atom should be roughly centrally located within the residue, and be invariant.
   *
   * @return A reference Atom.
   */
  public Atom getReferenceAtom() {
    Atom atom = null;
    switch (this.getResidueType()) {
      case AA:
        atom = (Atom) this.getAtomNode("CA");
        break;
      case NA:
        // If pyrimidine, atom will be N1.  Else, if purine,
        // N1 will return null, so grab N9.
        atom = (Atom) this.getAtomNode("N1");
        if (atom == null) {
          atom = (Atom) this.getAtomNode("N9");
        }
        break;
      default:
        break;
    }
    if (atom == null) {
      atom = (Atom) this.getAtomNode(0);
    }
    return atom;
  }

  /**
   * Returns this Residue's sequence number.
   *
   * @return a int.
   */
  public int getResidueNumber() {
    return resNumber;
  }

  /**
   * Getter for the field <code>residueType</code>.
   *
   * @return a {@link ResidueType} object.
   */
  public ResidueType getResidueType() {
    return residueType;
  }

  /**
   * Get the current rotamer.
   *
   * @return a {@link ffx.potential.bonded.Rotamer} object.
   */
  public Rotamer getRotamer() {
    return currentRotamer;
  }

  /**
   * Set the current rotamer.
   *
   * @param rotamer a {@link ffx.potential.bonded.Rotamer} object.
   */
  public void setRotamer(Rotamer rotamer) {
    this.currentRotamer = rotamer;
  }

  /**
   * Return all currently set rotamers.
   *
   * @return All current rotamers for this residue.
   */
  public Rotamer[] getRotamers() {
    return rotamers;
  }

  /**
   * Resets the rotamers for this residue, potentially incorporating the original coordinates if
   * RotamerLibrary's original coordinates rotamer flag has been set.
   * <p>
   * Any rotamers that were set previously are deleted.
   *
   * @param library Rotamer library to use
   * @return An array of Rotamer.
   */
  public Rotamer[] setRotamers(RotamerLibrary library) {

    // Obtain rotamers for this residue from the RotamerLibrary.
    Rotamer[] libRotamers = library.getRotamers(this, titrationUtils);
    rotamers = libRotamers;

    // Add the original coordinates as a rotamer if requested
    // (unless the library has no rotamers for this residue).
    if (library.getUsingOrigCoordsRotamer() && rotamers != null) {
      // Define the current coordinates as a new rotamer.
      Rotamer[] originalRotamers = Rotamer.defaultRotamerFactory(this, titrationUtils);
      int nOrig = originalRotamers.length;
      // Add the original rotamers to those from the library and cache the result.
      int libRots = libRotamers.length;
      rotamers = new Rotamer[libRots + nOrig];
      // First rotamers are the original conformation.
      arraycopy(originalRotamers, 0, rotamers, 0, nOrig);
      // Copy in the library.
      arraycopy(libRotamers, 0, rotamers, nOrig, libRots);
    }

    return rotamers;
  }

  /**
   * Getter for the field <code>segID</code>.
   *
   * @return a {@link java.lang.String} object.
   */
  public String getSegID() {
    return segID;
  }

  /**
   * Returns a list of side chain atoms; for our purposes, nucleic acid side chain atoms are the
   * sugar and the phosphate.
   *
   * @return List of side chain (or nucleic backbone) atoms.
   */
  public List<Atom> getSideChainAtoms() {
    List<Atom> atoms = getAtomList();
    List<Atom> ret;
    switch (residueType) {
      case NA:
        ret = new ArrayList<>();
        for (Atom atom : atoms) {
          String name = atom.getName().toUpperCase();
          /*
           * Very conveniently for our purposes, the entire sugar is
           * denoted with ' at the end. Note that we add the side chain
           * atoms for NAs, instead of removing the backbone.
           */
          if (name.contains("'") || name.equals("P") || name.startsWith("OP")) {
            ret.add(atom);
          }
        }
        return ret;
      case AA:
        ret = new ArrayList<>(atoms);
        for (Atom atom : atoms) {
          String name = atom.getName().toUpperCase();
          if (name.equals("N")
              || name.equals("H")
              || name.equals("H1")
              || name.equals("H2")
              || name.equals("H3")
              || name.equals("CA")
              || name.startsWith("HA")
              || name.equals("C")
              || name.equals("O")
              || name.equals("OXT")
              || name.equals("OT2")) {
            ret.remove(atom);
          }
        }
        return ret;
      default:
        return null;
    }
  }

  /**
   * Returns a list of atoms liable to change during dead-end elimination repacking. For ordinary
   * amino acids: side chain atoms. For ordinary nucleic acids: sugar/phosphate backbone atoms.
   * MultiResidue over-rides this to return all atoms (as backbone atom types are nonconstant).
   *
   * @return Atoms changeable during DEE.
   */
  public List<Atom> getVariableAtoms() {
    return getSideChainAtoms();
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    return Objects.hash(segID, getResidueNumber(), residueType, getName());
  }

  /**
   * Initializes this (presumably nucleic acid) Residue's C1s, O4s, C4s, O3sNorth, and O3sSouth
   * default coordinates based on default PDB atom locations; to preserve rotamer independence, this
   * must be called before any NA rotamers are applied.
   */
  public void initializeDefaultAtomicCoordinates() {
    if (residueType != ResidueType.NA) {
      return;
    }
    boolean isDeoxy;
    try {
      switch (NucleicAcidUtils.NucleicAcid3.valueOf(this.getName())) {
        case DAD:
        case DCY:
        case DGU:
        case DTY:
          isDeoxy = true;
          break;
        case CYT:
        case ADE:
        case THY:
        case URI:
        case GUA:
        default:
          isDeoxy = false;
          break;
      }
      C1sCoords = new double[3];
      ((Atom) getAtomNode("C1'")).getXYZ(C1sCoords);
      O4sCoords = new double[3];
      ((Atom) getAtomNode("O4'")).getXYZ(O4sCoords);
      C4sCoords = new double[3];
      ((Atom) getAtomNode("C4'")).getXYZ(C4sCoords);

      /*
       With the place flag set false, applySugarPucker returns
       hypothetical O3' coordinates based on default atom positions and
       the supplied sugar pucker.
      */
      O3sNorthCoords =
          RotamerLibrary.applySugarPucker(
              this, RotamerLibrary.NucleicSugarPucker.C3_ENDO, isDeoxy, false);
      O3sSouthCoords =
          RotamerLibrary.applySugarPucker(
              this, RotamerLibrary.NucleicSugarPucker.C2_ENDO, isDeoxy, false);
    } catch (Exception e) {
      logger.log(Level.WARNING, toString(), e);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Prints "Residue Number: x" to stdout.
   */
  @Override
  public void print() {
    logger.info(" " + this);
    for (Atom a : getAtomNode().getAtomList()) {
      a.print();
    }
  }

  /**
   * revertState.
   *
   * @param state a {@link ffx.potential.bonded.ResidueState} object.
   */
  public void revertState(ResidueState state) {
    List<Atom> atomList = getAtomList();
    for (Atom atom : atomList) {
      atom.moveTo(state.getAtomCoords(atom));
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
    // If Color by Residue, pass this Residue's Color
    if (newColorModel == RendererCache.ColorModel.RESIDUE) {
      switch (residueType) {
        case AA:
          color = AA3Color.get(aa);
          break;
        case NA:
          color = NA3Color.get(na);
          break;
        default:
          color = null;
      }
      if (color == null) {
        return;
      }
      mat = RendererCache.materialFactory(color);
    } else if (newColorModel == RendererCache.ColorModel.STRUCTURE) {
      color = SSTypeColor.get(SSType.NONE);
      mat = RendererCache.materialFactory(color);
    }
    super.setColor(newColorModel, color, mat);
  }

  /**
   * setNumber
   *
   * @param n a int.
   */
  public void setNumber(int n) {
    resNumber = n;
    for (Atom atom : getAtomList()) {
      atom.setResidueNumber(n);
    }
  }

  /**
   * storeCoordinateArray.
   *
   * @return an array of {@link double} objects.
   */
  public double[][] storeCoordinateArray() {
    List<Atom> atomList = getAtomList();
    int nAtoms = atomList.size();
    double[][] coords = new double[nAtoms][3];
    int i = 0;
    for (Atom atom : atomList) {
      atom.getXYZ(coords[i++]);
    }
    return coords;
  }

  /**
   * storeState.
   *
   * @return a {@link ffx.potential.bonded.ResidueState} object.
   */
  public ResidueState storeState() {
    return new ResidueState(this, this);
  }

  /**
   * Formats this residue with some optional inclusions.
   *
   * <p>[residue type]-[chain ID]ResNumber-Name.
   *
   * @param addResType Include the residue type
   * @param addChainID Include the chain ID.
   * @return A descriptive string.
   */
  public String toFormattedString(boolean addResType, boolean addChainID) {
    StringBuilder sb = new StringBuilder();
    if (addResType) {
      sb.append(residueType.toString()).append("-");
    }
    if (addChainID) {
      sb.append(chainID);
    }
    sb.append(this);
    return sb.toString();
  }

  /**
   * A descriptive string based on a given rotamer.
   *
   * <p>[chain ID]ResNumber-RotamerName.
   *
   * @param rotamer The rotamer to use.
   * @return A descriptive string.
   */
  public String toString(Rotamer rotamer) {
    return "" + chainID + resNumber + "-" + rotamer.getName();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    if (shortString == null) {
      shortString = "" + resNumber + "-" + getName();
    }
    return shortString;
  }

  /**
   * Add an array of rotamers to this Residue's cached array of rotamers.
   *
   * @param rotamers The rotamers to add.
   */
  void addRotamers(Rotamer[] rotamers) {
    for (Rotamer rotamer : rotamers) {
      addRotamer(rotamer);
    }
  }

  /**
   * Add a rotamer to this Residue's cached array of rotamers.
   *
   * @param rotamer The rotamer to add.
   */
  void addRotamer(Rotamer rotamer) {
    if (rotamers != null) {
      Rotamer[] libRotamers = rotamers;
      int nRots = libRotamers.length;
      rotamers = new Rotamer[nRots + 1];
      arraycopy(libRotamers, 0, rotamers, 0, nRots);
      rotamers[rotamers.length - 1] = rotamer;
    } else {
      rotamers = new Rotamer[1];
      rotamers[0] = rotamer;
    }
  }

  /**
   * deleteAtom
   *
   * @param atomToDelete a {@link ffx.potential.bonded.Atom} object.
   */
  void deleteAtom(Atom atomToDelete) {
    MSNode atoms = getAtomNode();
    if (atoms.contains(atomToDelete) != null) {
      logger.info(
          " The following atom is being deleted from the model:\n" + atomToDelete);
      atoms.remove(atomToDelete);
    }
  }

  /**
   * Returns the position of this (presumably nucleic acid) Residue's default O3' coordinates given a
   * North pucker.
   *
   * @return a new double[] with default XYZ coordinates for O3' in a North pucker.
   */
  double[] getO3sNorth() {
    double[] ret = new double[3];
    arraycopy(O3sNorthCoords, 0, ret, 0, ret.length);
    return ret;
  }

  /**
   * Returns the position of this (presumably nucleic acid) Residue's default O3' coordinates given a
   * South pucker.
   *
   * @return a new double[] with default XYZ coordinates for O3' in a South pucker.
   */
  double[] getO3sSouth() {
    double[] ret = new double[3];
    arraycopy(O3sSouthCoords, 0, ret, 0, ret.length);
    return ret;
  }

  /**
   * Returns the position of this (presumably nucleic acid) Residue's original C1' coordinates.
   *
   * @return a new double[] with original XYZ coordinates for C1'.
   */
  double[] getC1sCoords() {
    double[] ret = new double[3];
    arraycopy(C1sCoords, 0, ret, 0, ret.length);
    return ret;
  }

  /**
   * Returns the position of this (presumably nucleic acid) Residue's original O4' coordinates.
   *
   * @return a new double[] with original XYZ coordinates for O4'.
   */
  double[] getO4sCoords() {
    double[] ret = new double[3];
    arraycopy(O4sCoords, 0, ret, 0, ret.length);
    return ret;
  }

  /**
   * Returns the position of this (presumably nucleic acid) Residue's original C4' coordinates.
   *
   * @return a new double[] with original XYZ coordinates for C4'.
   */
  double[] getC4sCoords() {
    double[] ret = new double[3];
    arraycopy(C4sCoords, 0, ret, 0, ret.length);
    return ret;
  }

  /**
   * Uses a name to add an Atom to a List<Atom> if the Atom exists for this residue.
   *
   * @param atList List to add to.
   * @param name Atom to add.
   * @return If atom exists.
   */
  private boolean tryAddAtom(List<Atom> atList, String name) {
    try {
      Atom at = (Atom) getAtomNode(name);
      if (at != null) {
        atList.add(at);
        return true;
      } else {
        return false;
      }
    } catch (Exception ex) {
      return false;
    }
  }

  private void assignResidueType() {
    String name = getName().toUpperCase();
    switch (residueType) {
      case AA:
        aa = null;
        try {
          if (name.length() >= 2) {
            aa = AminoAcidUtils.AminoAcid3.valueOf(name);
          } else if (name.length() == 1) {
            AminoAcid1 aa1 = AminoAcidUtils.AminoAcid1.valueOf(name);
            aa = AA1toAA3.get(aa1);
          }
        } catch (Exception e) {
          logger.fine(String.format("Exception assigning AA3 for residue: %s", name));
          aa = AminoAcidUtils.AminoAcid3.UNK;
        }
        break;
      case NA:
        na = null;
        try {
          if (name.length() >= 2) {
            na = NucleicAcidUtils.NucleicAcid3.parse(name);
          } else if (name.length() == 1) {
            NucleicAcid1 na1 = NucleicAcidUtils.NucleicAcid1.valueOf(name);
            na = NA1toNA3.get(na1);
          }
        } catch (Exception e) {
          na = NucleicAcidUtils.NucleicAcid3.UNK;
        }
        break;
    }
  }

  public enum SSType {
    NONE,
    HELIX,
    SHEET,
    TURN
  }

  /**
   * Residue type [NA, AA, UNK].
   */
  public enum ResidueType {
    NA,
    AA,
    UNK
  }
}
