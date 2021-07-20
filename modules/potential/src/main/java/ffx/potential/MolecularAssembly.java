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
package ffx.potential;

import static ffx.numerics.math.DoubleMath.length;
import static ffx.numerics.math.DoubleMath.sub;
import static ffx.potential.extended.ExtUtils.prop;
import static java.lang.String.format;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Indexing;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSGroup;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.RendererCache;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.parameters.ForceField;
import ffx.utilities.StringUtils;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.jogamp.java3d.Appearance;
import org.jogamp.java3d.BoundingSphere;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.ColoringAttributes;
import org.jogamp.java3d.GeometryArray;
import org.jogamp.java3d.Group;
import org.jogamp.java3d.LineArray;
import org.jogamp.java3d.LineAttributes;
import org.jogamp.java3d.Link;
import org.jogamp.java3d.Material;
import org.jogamp.java3d.Node;
import org.jogamp.java3d.RenderingAttributes;
import org.jogamp.java3d.Shape3D;
import org.jogamp.java3d.SharedGroup;
import org.jogamp.java3d.Switch;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.utils.picking.PickTool;
import org.jogamp.vecmath.Color3f;
import org.jogamp.vecmath.Matrix3d;
import org.jogamp.vecmath.Point3d;
import org.jogamp.vecmath.Vector3d;

/**
 * The MolecularAssembly class is a collection of Polymers, Hetero Molecules, Ions and Water
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MolecularAssembly extends MSGroup {

  /** Constant <code>atomIndexing</code> */
  public static final Indexing atomIndexing =
      prop(Indexing.class, "sys.atomIndexing", Indexing.XYZ);

  private static final Logger logger = Logger.getLogger(MolecularAssembly.class.getName());
  private static final double[] a = new double[3];
  /** Persistent index parallel to xyzIndex. */
  public static int persistentAtomIndexer = 1;

  private final List<String> headerLines = new ArrayList<>();
  // Data Nodes
  private final MSNode ions = new MSNode("Ions");
  private final HashMap<String, Molecule> ionHashMap = new HashMap<>();
  private final MSNode water = new MSNode("Waters");
  private final HashMap<String, Molecule> waterHashMap = new HashMap<>();
  private final MSNode molecules = new MSNode("Hetero Molecules");
  private final HashMap<String, Molecule> moleculeHashMap = new HashMap<>();
  private final List<BranchGroup> myNewShapes = new ArrayList<>();
  protected ForceField forceField;
  // MolecularAssembly member variables
  private File file;
  private ForceFieldEnergy potentialEnergy;
  private CompositeConfiguration properties;
  private Vector3d offset;
  private int cycles = 1;
  private int currentCycle = 1;
  // 3D Graphics Nodes - There is a diagram explaining the MolecularAssembly
  // Scenegraph below
  private BranchGroup branchGroup;
  private TransformGroup originToRot;
  private Transform3D originToRotT3D;
  private Vector3d originToRotV3D;
  private TransformGroup rotToCOM;
  private Transform3D rotToCOMT3D;
  private Vector3d rotToCOMV3D;
  private BranchGroup base;
  private Switch switchGroup;
  private Shape3D wire;
  private BranchGroup vrml;
  private TransformGroup vrmlTG;
  private Transform3D vrmlTd;
  private BranchGroup childNodes;
  private Atom[] atomLookUp;
  private LineAttributes lineAttributes;
  private boolean visible = false;
  private FractionalMode fractionalMode = FractionalMode.MOLECULE;
  private double[][] fractionalCoordinates;

  /**
   * Constructor for MolecularAssembly.
   *
   * @param name a {@link java.lang.String} object.
   */
  public MolecularAssembly(String name) {
    super(name);
    getAtomNode().setName("MacroMolecules");
    add(molecules);
    add(ions);
    add(water);
  }

  /**
   * Constructor for MolecularAssembly.
   *
   * @param name a {@link java.lang.String} object.
   * @param Polymers a {@link ffx.potential.bonded.MSNode} object.
   */
  public MolecularAssembly(String name, MSNode Polymers) {
    super(name, Polymers);
  }

  /**
   * Constructor for MolecularAssembly.
   *
   * @param name a {@link java.lang.String} object.
   * @param Polymers a {@link ffx.potential.bonded.MSNode} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public MolecularAssembly(String name, MSNode Polymers, CompositeConfiguration properties) {
    this(name, Polymers);
    this.properties = properties;
  }

  /**
   * Adds a header line to this MolecularAssembly (particularly for PDB formats)
   *
   * @param line Line to add.
   */
  public void addHeaderLine(String line) {
    headerLines.add(line);
  }

  /** {@inheritDoc} */
  @Override
  public MSNode addMSNode(MSNode o) {
    List<MSNode> Polymers = getAtomNodeList();
    if (o instanceof Atom) {
      Atom atom = (Atom) o;
      if (atom.isModRes()) {
        return getResidue(atom, true, Residue.ResidueType.AA);
      } else if (!atom.isHetero()) {
        return getResidue(atom, true);
      } else {
        return getMolecule(atom, true);
      }
    } else if (o instanceof Residue) {
      Residue residue = (Residue) o;
      Character chainID = residue.getChainID();
      String segID = residue.getSegID();
      int index = Polymers.indexOf(new Polymer(chainID, segID));
      // See if the polymer already exists.
      if (index != -1) {
        Polymer c = (Polymer) Polymers.get(index);
        setFinalized(false);
        return c.addMSNode(residue);
      } else {
        Polymer newc = new Polymer(chainID, segID);
        getAtomNode().add(newc);
        setFinalized(false);
        return newc.addMSNode(residue);
      }
    } else if (o instanceof Polymer) {
      Polymer c = (Polymer) o;
      int index = Polymers.indexOf(c);
      if (index == -1) {
        getAtomNode().add(c);
        setFinalized(false);
        return c;
      } else {
        return Polymers.get(index);
      }
    } else if (o instanceof Molecule) {
      Molecule m = (Molecule) o;
      String key = m.getKey();
      if (m.getAtomNode().getChildCount() == 1) {
        ions.add(m);
        if (ionHashMap.containsKey(key)) {
          logger.info(" Ion map already contains " + m);
        } else {
          ionHashMap.put(m.getKey(), m);
        }
        return m;
      } else if (Utilities.isWaterOxygen((Atom) m.getAtomNode().getChildAt(0))) {
        water.add(m);
        if (waterHashMap.containsKey(key)) {
          logger.info(" Water map already contains " + m);
        } else {
          waterHashMap.put(m.getKey(), m);
        }
        return m;
      } else {
        molecules.add(m);
        if (moleculeHashMap.containsKey(key)) {
          logger.info(" Molecule map already contains " + m);
        } else {
          moleculeHashMap.put(m.getKey(), m);
        }
        return m;
      }
    } else {
      String message = "Programming error in MolecularAssembly addNode";
      logger.log(Level.SEVERE, message);
      return o;
    }
  }

  /**
   * Applies a randomly drawn density to a molecular system's crystal.
   *
   * @param ucDensity Target density.
   */
  public void applyRandomDensity(double ucDensity) {
    if (ucDensity > 0) {
      logger.info(
          format("\n Applying random unit cell axes with target density of %6.3f\n", ucDensity));
      // The replicates crystal is needed here (not the unit cell).
      Crystal crystal = getCrystal();
      if (!crystal.aperiodic()) {
        double mass = getMass();
        crystal.randomParameters(ucDensity, mass);
        potentialEnergy.setCrystal(crystal);
      } else {
        logger.fine(String.format(" Potential %s is an aperiodic system!", potentialEnergy));
      }
    }
  }

  /**
   * Randomizes position in the unit cell of each molecule by applying a Cartesian SymOp with a
   * random translation.
   *
   * @param x SymOp with translation range -x/2 .. x/2 (0 for random placement in the unit cell,
   *     negative for no SymOp)
   */
  public void applyRandomSymOp(double x) {
    // The Unit Cell is needed here (not the ReplicatesCrystal).
    Crystal crystal = getCrystal().getUnitCell();
    if (crystal.aperiodic() || x < 0) {
      return;
    }
    double[] xyz = new double[3];
    List<MSNode> molecules = getMolecules();
    int moleculeNum = 1;
    for (MSNode msNode : molecules) {
      List<Atom> atoms = msNode.getAtomList();
      SymOp symOp;
      if (x == 0) {
        double[] translation = crystal.getRandomCartTranslation();
        symOp = SymOp.randomSymOpFactory(translation);
      } else {
        symOp = SymOp.randomSymOpFactory(x);
      }
      logger.info(
          format(
              "\n Applying random Cartesian SymOp to molecule %d:\n%s",
              moleculeNum, symOp));
      for (Atom atom : atoms) {
        atom.getXYZ(xyz);
        crystal.applyCartesianSymOp(xyz, xyz, symOp);
        atom.setXYZ(xyz);
      }
      moleculeNum++;
    }
  }

  /** center */
  public void center() {
    double[] center = getMultiScaleCenter(false);
    offset = new Vector3d(center);
    if (vrml != null) {
      vrmlTd.set(offset);
      vrmlTG.setTransform(vrmlTd);
    }
    offset.negate();
    originToRotV3D.set(offset);
    originToRotT3D.setTranslation(originToRotV3D);
    originToRot.setTransform(originToRotT3D);
    rotToCOMT3D.setIdentity();
    rotToCOM.setTransform(rotToCOMT3D);
    offset.negate();
    rotateAbout(offset);
    originToRotT3D.get(offset);
  }

  /**
   * centerAt
   *
   * @param d an array of double.
   */
  public void centerAt(double[] d) {
    double[] Rc = {0, 0, 0};
    double[] c = new double[3];
    ListIterator<Atom> li;
    int i, num = getAtomList().size();
    for (li = getAtomList().listIterator(); li.hasNext(); ) {
      (li.next()).getXYZ(a);
      Rc[0] += a[0];
      Rc[1] += a[1];
      Rc[2] += a[2];
    }
    for (i = 0; i < 3; i++) {
      Rc[i] /= num;
    }
    sub(d, Rc, c);
    for (li = getAtomList().listIterator(); li.hasNext(); ) {
      (li.next()).move(c);
    }
  }

  /**
   * centerView
   *
   * @param rot a boolean.
   * @param trans a boolean.
   */
  public void centerView(boolean rot, boolean trans) {
    originToRot.getTransform(originToRotT3D);
    if (rot) {
      Matrix3d m3d = new Matrix3d();
      m3d.setIdentity();
      originToRotT3D.setRotation(m3d);
      // rotToCOMT3D.setRotation(m3d);
    }
    if (trans) {
      originToRotV3D.set(offset);
      originToRotT3D.set(originToRotV3D);
    }
    originToRot.setTransform(originToRotT3D);
    // rotToCOM.setTransform(rotToCOMT3D);
  }

  /** Compute fractional coordinates. */
  public void computeFractionalCoordinates() {

    // Count up the number of fractional coordinate entities.
    fractionalCount();

    Crystal unitCell = getCrystal().getUnitCell();
    double[] com = new double[3];

    switch (fractionalMode) {
      case MOLECULE:
        int iMolecule = 0;
        Polymer[] polymers = getChains();
        if (polymers != null && polymers.length > 0) {
          // Find the center of mass
          for (Polymer polymer : polymers) {
            List<Atom> list = polymer.getAtomList();
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
              double m = atom.getMass();
              com[0] += atom.getX() * m;
              com[1] += atom.getY() * m;
              com[2] += atom.getZ() * m;
              totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
          }
        }

        // Loop over each molecule
        List<MSNode> molecules = getMolecules();
        for (MSNode molecule : molecules) {
          List<Atom> list = molecule.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
        }

        // Loop over each water
        List<MSNode> waters = getWaters();
        for (MSNode water : waters) {
          List<Atom> list = water.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
        }

        // Loop over each ion
        List<MSNode> ions = getIons();
        for (MSNode ion : ions) {
          List<Atom> list = ion.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
        }
        break;
      case ATOM:
        Atom[] atoms = getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
          atoms[i].getXYZ(com);
          unitCell.toFractionalCoordinates(com, fractionalCoordinates[i]);
        }
        break;
      case OFF:
        break;
    }
  }

  /** createBox */
  public void createBox() {
    int vertices = 8;
    LineArray la =
        new LineArray(
            4 * vertices,
            GeometryArray.COORDINATES | GeometryArray.COLOR_4 | GeometryArray.NORMALS);
    la.setCapability(LineArray.ALLOW_COORDINATE_WRITE);
    la.setCapability(LineArray.ALLOW_COORDINATE_READ);
    la.setCapability(LineArray.ALLOW_COLOR_WRITE);
    la.setCapability(LineArray.ALLOW_COUNT_READ);
    la.setCapability(LineArray.ALLOW_INTERSECT);
    la.setCapability(LineArray.ALLOW_FORMAT_READ);

    ColoringAttributes cola =
        new ColoringAttributes(new Color3f(), ColoringAttributes.SHADE_GOURAUD);
    Appearance app = new Appearance();
    lineAttributes = new LineAttributes();
    lineAttributes.setLineWidth(RendererCache.bondwidth);
    lineAttributes.setCapability(LineAttributes.ALLOW_WIDTH_WRITE);
    lineAttributes.setLineAntialiasingEnable(true);
    app.setLineAttributes(lineAttributes);
    app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_READ);
    app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_WRITE);
    RenderingAttributes ra = new RenderingAttributes();
    ra.setAlphaTestValue(0.1f);
    ra.setAlphaTestFunction(RenderingAttributes.GREATER);
    ra.setDepthBufferEnable(true);
    ra.setDepthBufferWriteEnable(true);
    app.setRenderingAttributes(ra);
    app.setColoringAttributes(cola);
    Shape3D wireframe = new Shape3D(la, app);
    // PickTool.setCapabilities(wire, PickTool.INTERSECT_COORD);
    wireframe.setUserData(this);
    wireframe.setBounds(new BoundingSphere(new Point3d(0, 0, 0), 10.0));
    try {
      wireframe.setBoundsAutoCompute(false);
    } catch (Exception e) {
      e.printStackTrace();
    }
    wireframe.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
    wireframe.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
    // return wire;
  }

  /**
   * deleteMolecule
   *
   * @param molecule a {@link ffx.potential.bonded.Molecule} object.
   */
  public void deleteMolecule(Molecule molecule) {
    List<MSNode> list = ions.getChildList();
    for (MSNode node : list) {
      Molecule m = (Molecule) node;
      if (molecule == m) {
        ions.remove(m);
        ionHashMap.remove(m.getKey());
        return;
      }
    }
    list = water.getChildList();
    for (MSNode node : list) {
      Molecule m = (Molecule) node;
      if (molecule == m) {
        water.remove(m);
        waterHashMap.remove(m.getKey());
        return;
      }
    }
    list = molecules.getChildList();
    for (MSNode node : list) {
      Molecule m = (Molecule) node;
      if (molecule == m) {
        molecules.remove(m);
        moleculeHashMap.remove(m.getKey());
        return;
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean destroy() {
    try {
      if (potentialEnergy == null) {
        finishDestruction();
        return true;
      } else {
        // potentialEnergy.destroy() will loop around to call finishDestruction().
        // This is a poor construction, but is necessary due to prior decisions.
        return potentialEnergy.destroy();
      }
    } catch (Exception ex) {
      logger.warning(String.format(" Exception in destroying a MolecularAssembly: %s", ex));
      logger.info(Utilities.stackTraceToString(ex));
      return false;
    }
  }

  /** detach */
  public void detach() {
    synchronized (this) {
      if (branchGroup != null && branchGroup.isLive()) {
        branchGroup.detach();
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public void finalize(boolean finalizeGroups, ForceField forceField) {
    setFinalized(false);
    if (finalizeGroups) {
      bondTime = 0;
      angleTime = 0;
      stretchBendTime = 0;
      ureyBradleyTime = 0;
      outOfPlaneBendTime = 0;
      torsionTime = 0;
      piOrbitalTorsionTime = 0;
      torsionTorsionTime = 0;
      List<MSNode> Polymers = getAtomNodeList();
      for (MSNode msNode : Polymers) {
        MSGroup group = (MSGroup) msNode;
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(" Finalizing bonded terms for polymer " + group.toString());
        }
        try {
          group.finalize(true, forceField);
        } catch (Exception e) {
          String message = "Fatal exception finalizing " + group.toString();
          logger.log(Level.SEVERE, message, e);
          System.exit(-1);
        }
        if (logger.isLoggable(Level.FINE)) {
          Runtime runtime = Runtime.getRuntime();
          long occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
          long MB = 1024 * 1024;
          logger.fine(
              "\n In-Use Memory   (Mb): "
                  + occupiedMemory / MB
                  + "\n Free Memory     (Mb): "
                  + runtime.freeMemory() / MB
                  + "\n Total Memory    (Mb): "
                  + runtime.totalMemory() / MB);
        }
      }
      for (MSNode m : molecules.getChildList()) {
        Molecule molecule = (Molecule) m;
        molecule.finalize(true, forceField);
      }
      for (MSNode m : water.getChildList()) {
        Molecule molecule = (Molecule) m;
        molecule.finalize(true, forceField);
      }
      for (MSNode m : ions.getChildList()) {
        Molecule molecule = (Molecule) m;
        molecule.finalize(true, forceField);
      }
      if (logger.isLoggable(Level.FINE)) {
        StringBuilder sb = new StringBuilder("\n Time to create bonded energy terms\n\n");
        sb.append(format(" Bond Streching     %10.3f\n", bondTime * 1.0e-9));
        sb.append(format(" Angle Bending      %10.3f\n", angleTime * 1.0e-9));
        sb.append(format(" Stretch-Bend       %10.3f\n", stretchBendTime * 1.0e-9));
        sb.append(format(" Urey-Bradley       %10.3f\n", ureyBradleyTime * 1.0e-9));
        sb.append(format(" Out-of-Plane Bend  %10.3f\n", outOfPlaneBendTime * 1.0e-9));
        sb.append(format(" Torsionanl Angle   %10.3f\n", torsionTime * 1.0e-9));
        sb.append(format(" Pi-Orbital Torsion %10.3f\n", piOrbitalTorsionTime * 1.0e-9));
        sb.append(format(" Torsion-Torsion    %10.3f\n", torsionTorsionTime * 1.0e-9));
        logger.fine(sb.toString());
      }
    }
    if (!GraphicsEnvironment.isHeadless()) {
      createScene(!finalizeGroups);
      center();
    }
    removeLeaves();

    // Apply the Heavy Hydrogen flag.
    boolean heavyHydrogen = forceField.getBoolean("HEAVY_HYDROGENS", false);
    heavyHydrogen = forceField.getBoolean("HEAVY_HYDROGEN", heavyHydrogen);
    if (heavyHydrogen) {
      applyHeavyHydrogen();
    }

    setFinalized(true);
  }

  /**
   * findAtom
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @return a {@link ffx.potential.bonded.Atom} object.
   */
  public Atom findAtom(Atom atom) {
    if (!atom.isHetero() || atom.isModRes()) {
      Polymer polymer = getPolymer(atom.getChainID(), atom.getSegID(), false);
      if (polymer != null) {
        Residue res = polymer.getResidue(atom.getResidueName(), atom.getResidueNumber(), false);
        if (res != null) {
          MSNode node = res.getAtomNode();
          return (Atom) node.contains(atom);
        }
      }
      return null;
    } else {
      String resName = atom.getResidueName();
      int resNum = atom.getResidueNumber();
      String segID = atom.getSegID();
      String key = resNum + resName + segID;
      Molecule m = ionHashMap.get(key);
      if (m == null) {
        m = waterHashMap.get(key);
        if (m == null) {
          m = moleculeHashMap.get(key);
          if (m == null) {
            return null;
          }
        }
      }
      return (Atom) m.contains(atom);
    }
  }

  /**
   * Count the number of fractional coordinate entities in the system. If fractionalMode is
   * MOLECULE, then the count is equal to the number of molecules. If fractionalMode is ATOM, then
   * the count is the number of atoms. Otherwise the count is zero.
   *
   * @return The number of fractional coordinate entities.
   */
  public int fractionalCount() {
    int count = 0;
    switch (fractionalMode) {
      case MOLECULE:
        // Move polymers togethers.
        Polymer[] polymers = getChains();
        if (polymers != null && polymers.length > 0) {
          count += polymers.length;
        }
        List<MSNode> molecules = getMolecules();
        if (molecules != null) {
          count += molecules.size();
        }

        List<MSNode> waters = getWaters();
        if (waters != null) {
          count += waters.size();
        }

        List<MSNode> ions = getIons();
        if (ions != null) {
          count += ions.size();
        }
        break;
      case ATOM:
        count = getAtomArray().length;
        break;
      case OFF:
        count = 0;
        break;
    }

    if (fractionalCoordinates == null || fractionalCoordinates.length != count) {
      fractionalCoordinates = new double[count][3];
    }

    return count;
  }

  /**
   * getActiveAtomArray
   *
   * @return an array of active {@link ffx.potential.bonded.Atom} objects.
   */
  public Atom[] getActiveAtomArray() {
    List<Atom> atoms = getAtomList();
    List<Atom> activeAtoms = new ArrayList<>();
    for (Atom a : atoms) {
      if (a.isActive()) {
        activeAtoms.add(a);
      }
    }
    Atom[] atomArray = activeAtoms.toArray(new Atom[0]);
    Arrays.sort(atomArray);
    return atomArray;
  }

  /**
   * Gets all bonded entities in this MolecularAssembly, where an entity can be a polymer, molecule,
   * monoatomic ion, or monoatomic gas (i.e. noble gas atoms).
   *
   * @return All bonded groups of atoms and all singleton atoms.
   */
  public List<MSNode> getAllBondedEntities() {
    // Initial construction via HashSet to eliminate duplicates.
    Set<MSNode> allBondedNodes = new HashSet<>(getIons());
    allBondedNodes.addAll(getMolecules());
    allBondedNodes.addAll(getWaters());
    Polymer[] polys = getChains();
    if (polys != null && polys.length > 0) {
      allBondedNodes.addAll(Arrays.asList(polys));
    }
    return new ArrayList<>(allBondedNodes);
  }

  /**
   * getAltLocations
   *
   * @return an array of {@link java.lang.String} objects.
   */
  public String[] getAltLocations() {
    return null;
  }

  /**
   * getAtomArray
   *
   * @return an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public Atom[] getAtomArray() {
    List<Atom> atoms = getAtomList();
    Atom[] atomArray = atoms.toArray(new Atom[0]);
    for (int i = 0; i < atoms.size(); i++) {
      atomArray[i].setXyzIndex(i + 1);
    }

    return atomArray;
  }

  /**
   * getAtomFromWireVertex
   *
   * @param i a int.
   * @return a {@link ffx.potential.bonded.Atom} object.
   */
  public Atom getAtomFromWireVertex(int i) {
    if (atomLookUp != null && atomLookUp.length > i) {
      return atomLookUp[i];
    }
    return null;
  }

  /**
   * getBackBoneAtoms
   *
   * @return a {@link java.util.List} object.
   */
  public List<Atom> getBackBoneAtoms() {
    List<Atom> backbone = new ArrayList<>();
    List<Residue> residues = getResidueList();
    for (Residue residue : residues) {
      backbone.addAll(residue.getBackboneAtoms());
    }
    return backbone;
  }

  /**
   * Getter for the field <code>branchGroup</code>.
   *
   * @return a {@link org.jogamp.java3d.BranchGroup} object.
   */
  public BranchGroup getBranchGroup() {
    return branchGroup;
  }

  /**
   * getChain
   *
   * @param name a {@link java.lang.String} object.
   * @return a {@link ffx.potential.bonded.Polymer} object.
   */
  public Polymer getChain(String name) {
    for (MSNode node : getAtomNodeList()) {
      if (node instanceof Polymer) {
        String chainName = node.getName();
        if (chainName.equalsIgnoreCase(name)) {
          return (Polymer) node;
        } else if (name.contentEquals(" ")) {
          return (Polymer) node;
        }
        /* TODO: Right now if a molecular assembly has no chain ID, the first node is returned.
        This functionality will not work for cases where multiple chains exist and one of those chains has no name.*/
      }
    }
    return null;
  }

  /**
   * getChainNames
   *
   * @return an array of {@link java.lang.String} objects.
   */
  public String[] getChainNames() {
    List<String> temp = new ArrayList<>();
    for (MSNode node : getAtomNodeList()) {
      if (node instanceof Polymer) {
        temp.add(node.getName());
      }
    }
    if (temp.isEmpty()) {
      return null;
    }

    String[] names = new String[temp.size()];
    for (int i = 0; i < temp.size(); i++) {
      names[i] = temp.get(i);
    }

    return names;
  }

  /**
   * getChains
   *
   * @return an array of {@link ffx.potential.bonded.Polymer} objects.
   */
  public Polymer[] getChains() {
    List<Polymer> polymers = new ArrayList<>();
    for (MSNode node : getAtomNodeList()) {
      if (node instanceof Polymer) {
        polymers.add((Polymer) node);
      }
    }
    if (polymers.isEmpty()) {
      return null;
    }
    return polymers.toArray(new Polymer[0]);
  }

  /**
   * Sums up charge of the system, checking nonstandard residues for non-unitary charges.
   *
   * @param alwaysLog Log non-unitary charge warnings for all nodes
   * @return System charge
   */
  public double getCharge(boolean alwaysLog) {
    double totalCharge = 0;
    for (MSNode node : getNodeList()) {
      double charge = 0;
      boolean isNonstandard = false;
      for (Atom atom : node.getAtomList()) {
        charge += atom.getCharge(forceField);
        if (atom.isModRes()) {
          isNonstandard = true;
        }
      }
      if ((alwaysLog || isNonstandard) && (Math.abs(Math.round(charge) - charge) > 1.0E-5)) {
        logger.warning(
            String.format(" Node %s has non-unitary charge %12.8f", node, charge));
      }
      totalCharge += charge;
    }
    return totalCharge;
  }

  /**
   * getCrystal
   *
   * @return a {@link ffx.crystal.Crystal} object.
   */
  public Crystal getCrystal() {
    if (potentialEnergy == null) {
      return null;
    }
    return potentialEnergy.getCrystal();
  }

  /**
   * Set the Crystal for the Potential of this MolecularAssembly.
   *
   * @param crystal Crystal instance.
   */
  public void setCrystal(Crystal crystal) {
    if (potentialEnergy != null) {
      potentialEnergy.setCrystal(crystal);
    }
  }

  /**
   * Getter for the field <code>currentCycle</code>.
   *
   * @return a int.
   */
  public int getCurrentCycle() {
    return currentCycle;
  }

  /**
   * Setter for the field <code>currentCycle</code>.
   *
   * @param c a int.
   */
  public void setCurrentCycle(int c) {
    if (c <= cycles && c > 0) {
      currentCycle = c;
      for (Atom atom : getAtomList()) {
        atom.setCurrentCycle(currentCycle);
      }
    }
  }

  /**
   * Getter for the field <code>cycles</code>.
   *
   * @return a int.
   */
  public int getCycles() {
    return cycles;
  }

  /**
   * Setter for the field <code>cycles</code>.
   *
   * @param c a int.
   */
  public void setCycles(int c) {
    cycles = c;
  }

  /** {@inheritDoc} */
  @Override
  public double getExtent() {
    double[] Rc = {0, 0, 0};
    int num = getAtomList().size();
    for (Atom atom : getAtomList()) {
      atom.getXYZ(a);
      Rc[0] += a[0];
      Rc[1] += a[1];
      Rc[2] += a[2];
    }

    for (int i = 0; i < 3; i++) {
      Rc[i] /= num;
    }

    double r, d = 0;
    double[] xyz = new double[3];
    for (Atom atom : getAtomList()) {
      atom.getXYZ(xyz);
      sub(xyz, Rc, xyz);
      r = length(xyz);
      if (d < r) {
        d = r;
      }
    }
    return d;
  }

  /**
   * Getter for the field <code>file</code>.
   *
   * @return a {@link java.io.File} object.
   */
  public File getFile() {
    return file;
  }

  /**
   * Setter for the field <code>file</code>.
   *
   * @param f a {@link java.io.File} object.
   */
  public void setFile(File f) {
    if (f == null) {
      return;
    }
    file = f;
  }

  /**
   * Getter for the field <code>forceField</code>.
   *
   * @return a {@link ffx.potential.parameters.ForceField} object.
   */
  public ForceField getForceField() {
    return forceField;
  }

  /**
   * Setter for the field <code>forceField</code>.
   *
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   */
  public void setForceField(ForceField forceField) {
    this.forceField = forceField;
  }

  /**
   * Getter for the field <code>fractionalMode</code>.
   *
   * @return a {@link ffx.potential.MolecularAssembly.FractionalMode} object.
   */
  public FractionalMode getFractionalMode() {
    return fractionalMode;
  }

  /**
   * Setter for the field <code>fractionalMode</code>.
   *
   * @param mode a {@link ffx.potential.MolecularAssembly.FractionalMode} object.
   */
  public void setFractionalMode(FractionalMode mode) {
    fractionalMode = mode;
  }

  /**
   * Gets the header lines associated with this MolecularAssembly (particularly for PDB)
   *
   * @return Header lines.
   */
  public String[] getHeaderLines() {
    String[] ret = new String[headerLines.size()];
    headerLines.toArray(ret);
    return ret;
  }

  /**
   * Getter for the field <code>ions</code>.
   *
   * @return a {@link java.util.List} object.
   */
  public List<MSNode> getIons() {
    return ions.getChildList();
  }

  /**
   * getMass.
   *
   * @return a double.
   */
  public double getMass() {
    double mass = 0;
    for (Atom atom : getAtomArray()) {
      mass += atom.getMass();
    }
    return mass;
  }

  /**
   * This method assigns a unique integer to every molecule in the MolecularAssembly beginning at 0.
   * An integer array with these values for each atom is returned.
   *
   * @return an array of molecule numbers for each atom.
   */
  public int[] getMoleculeNumbers() {
    int[] moleculeNumber = new int[getAtomList().size()];
    int current = 0;
    // Loop over polymers together
    Polymer[] polymers = getChains();
    if (polymers != null && polymers.length > 0) {
      for (Polymer polymer : polymers) {
        List<Atom> atomList = polymer.getAtomList();
        for (Atom atom : atomList) {
          moleculeNumber[atom.getXyzIndex() - 1] = current;
        }
        current++;
      }
    }

    // Loop over each molecule
    for (MSNode molecule : molecules.getChildList()) {
      List<Atom> atomList = molecule.getAtomList();
      for (Atom atom : atomList) {
        moleculeNumber[atom.getXyzIndex() - 1] = current;
        atom.setMoleculeNumber(current);
      }
      current++;
    }

    // Loop over each water
    for (MSNode wat : water.getChildList()) {
      List<Atom> atomList = wat.getAtomList();
      for (Atom atom : atomList) {
        moleculeNumber[atom.getXyzIndex() - 1] = current;
      }
      current++;
    }

    // Loop over each ion
    for (MSNode ion : ions.getChildList()) {
      List<Atom> atomList = ion.getAtomList();
      for (Atom atom : atomList) {
        moleculeNumber[atom.getXyzIndex() - 1] = current;
      }
      current++;
    }

    return moleculeNumber;
  }

  /**
   * Getter for the field <code>molecules</code>.
   *
   * @return a {@link java.util.List} object.
   */
  public List<MSNode> getMolecules() {
    return molecules.getChildList();
  }

  /**
   * getNodeList
   *
   * @return a {@link java.util.List} object.
   */
  public List<MSNode> getNodeList() {
    List<MSNode> residues = new ArrayList<>();
    ListIterator<MSNode> li, lj;
    MSNode o;
    Polymer c;
    for (li = getAtomNodeList().listIterator(); li.hasNext(); ) {
      o = li.next();
      if (o instanceof Polymer) {
        c = (Polymer) o;
        for (lj = c.getAtomNodeList().listIterator(); lj.hasNext(); ) {
          o = lj.next();
          if (o instanceof Residue) {
            residues.add(o);
          }
        }
      }
    }

    residues.addAll(ions.getChildList());
    residues.addAll(water.getChildList());
    residues.addAll(molecules.getChildList());

    return residues;
  }

  /**
   * Getter for the field <code>offset</code>.
   *
   * @return a {@link org.jogamp.vecmath.Vector3d} object.
   */
  public Vector3d getOffset() {
    if (offset == null) {
      offset = new Vector3d(0.0, 0.0, 0.0);
    }

    return offset;
  }

  /**
   * Setter for the field <code>offset</code>.
   *
   * @param o a {@link org.jogamp.vecmath.Vector3d} object.
   */
  public void setOffset(Vector3d o) {
    offset = o;
  }

  /**
   * Getter for the field <code>originToRot</code>.
   *
   * @return a {@link org.jogamp.java3d.TransformGroup} object.
   */
  public TransformGroup getOriginToRot() {
    return originToRot;
  }

  /**
   * getParallelTeam.
   *
   * @return a {@link edu.rit.pj.ParallelTeam} object.
   */
  public ParallelTeam getParallelTeam() {
    if (potentialEnergy != null) {
      return potentialEnergy.getParallelTeam();
    } else {
      return null;
    }
  }

  /**
   * getPolymer
   *
   * @param chainID a {@link java.lang.Character} object.
   * @param segID a {@link java.lang.String} object.
   * @param create a boolean.
   * @return a {@link ffx.potential.bonded.Polymer} object.
   */
  public Polymer getPolymer(Character chainID, String segID, boolean create) {
    for (MSNode node : getAtomNodeList()) {
      if (node instanceof Polymer) {
        Polymer polymer = (Polymer) node;
        if (polymer.getName().equals(segID) && polymer.getChainID().equals(chainID)) {
          return (Polymer) node;
        }
      }
    }
    if (create) {
      Polymer polymer = new Polymer(chainID, segID, true);
      addMSNode(polymer);
      return polymer;
    }

    return null;
  }

  /**
   * Getter for the field <code>potentialEnergy</code>.
   *
   * @return a {@link ffx.potential.ForceFieldEnergy} object.
   */
  public ForceFieldEnergy getPotentialEnergy() {
    return potentialEnergy;
  }

  /**
   * Getter for the field <code>properties</code>.
   *
   * @return a {@link org.apache.commons.configuration2.CompositeConfiguration} object.
   */
  public CompositeConfiguration getProperties() {
    return properties == null ? forceField.getProperties() : properties;
  }

  /**
   * getResidueList
   *
   * @return a {@link java.util.List} object.
   */
  public List<Residue> getResidueList() {
    List<Residue> residues = new ArrayList<>();
    ListIterator<MSNode> li, lj;
    MSNode o;
    Polymer c;
    for (li = getAtomNodeList().listIterator(); li.hasNext(); ) {
      o = li.next();
      if (o instanceof Polymer) {
        c = (Polymer) o;
        for (lj = c.getAtomNodeList().listIterator(); lj.hasNext(); ) {
          o = lj.next();
          if (o instanceof Residue) {
            residues.add((Residue) o);
          }
        }
      }
    }
    return residues;
  }

  /**
   * getTransformGroup
   *
   * @return a {@link org.jogamp.java3d.TransformGroup} object.
   */
  public TransformGroup getTransformGroup() {
    return originToRot;
  }

  /**
   * getWaters
   *
   * @return a {@link java.util.List} object.
   */
  public List<MSNode> getWaters() {
    return water.getChildList();
  }

  /**
   * getWireFrame
   *
   * @return a {@link org.jogamp.java3d.Node} object.
   */
  public Node getWireFrame() {
    return wire;
  }

  /**
   * isVisible
   *
   * @return a boolean.
   */
  public boolean isVisible() {
    return visible;
  }

  /**
   * loadVRML
   *
   * @return a {@link org.jogamp.java3d.BranchGroup} object.
   */
  public BranchGroup loadVRML() {
    //        try {
    //            VrmlLoader loader = new VrmlLoader();
    //            VrmlScene scene = null;
    //            if (vrmlFile != null && vrmlFile.exists()) {
    //                scene = (VrmlScene) loader.load(vrmlFile.getAbsolutePath());
    //            } else if (vrmlURL != null) {
    //                scene = (VrmlScene) loader.load(vrmlURL);
    //            } else {
    //                return null;
    //            }
    //            BranchGroup bg = scene.getSceneGroup();
    //            recurseVRML(bg);
    //            bg.setCapability(BranchGroup.ALLOW_DETACH);
    //            bg.setCapability(BranchGroup.ALLOW_BOUNDS_READ);
    //            bg.compile();
    //            return bg;
    //        } catch (Exception e) {
    //            String message = "Fatal exception loading VRML.\n";
    //            logger.log(Level.SEVERE, message, e);
    //            System.exit(-1);
    //            return null;
    //        }
    return null;
  }

  /**
   * Moves the center of all chemical entities into the primary unit cell. Somewhat experimental
   * feature; use with caution.
   */
  public void moveAllIntoUnitCell() {
    moveIntoUnitCell(getChains());
    moveIntoUnitCell(getWaters());
    moveIntoUnitCell(getIons());
    moveIntoUnitCell(getMolecules());
  }

  /**
   * moveCenter
   *
   * @param d an array of double.
   */
  public void moveCenter(double[] d) {
    for (Atom atom : getAtomList()) {
      atom.move(d);
    }
  }

  /** Move to fractional coordinates. */
  public void moveToFractionalCoordinates() {

    if (fractionalCoordinates == null) {
      return;
    }

    Crystal unitCell = getCrystal().getUnitCell();
    double[] com = new double[3];

    switch (fractionalMode) {
      case MOLECULE:
        int iMolecule = 0;
        Polymer[] polymers = getChains();
        if (polymers != null && polymers.length > 0) {
          // Find the center of mass
          for (Polymer polymer : polymers) {
            List<Atom> list = polymer.getAtomList();
            double totalMass = 0.9;
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            for (Atom atom : list) {
              double m = atom.getMass();
              com[0] += atom.getX() * m;
              com[1] += atom.getY() * m;
              com[2] += atom.getZ() * m;
              totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            // Find the new center of mass in fractional coordinates.
            unitCell.toFractionalCoordinates(com, com);
            // Find the reciprocal translation vector.
            double[] frac = fractionalCoordinates[iMolecule++];
            com[0] = frac[0] - com[0];
            com[1] = frac[1] - com[1];
            com[2] = frac[2] - com[2];
            // Convert the fractional translation vector to Cartesian coordinates.
            unitCell.toCartesianCoordinates(com, com);
            // Move all atoms.
            for (Atom atom : list) {
              atom.move(com);
            }
          }
        }

        // Loop over each molecule
        List<MSNode> molecules = getMolecules();
        for (MSNode molecule : molecules) {
          List<Atom> list = molecule.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          // Find the new center of mass in fractional coordinates.
          unitCell.toFractionalCoordinates(com, com);
          // Find the reciprocal translation vector to the previous COM.
          double[] frac = fractionalCoordinates[iMolecule++];
          com[0] = frac[0] - com[0];
          com[1] = frac[1] - com[1];
          com[2] = frac[2] - com[2];
          // Convert the fractional translation vector to Cartesian coordinates.
          unitCell.toCartesianCoordinates(com, com);
          // Move all atoms.
          for (Atom atom : list) {
            atom.move(com);
          }
        }

        // Loop over each water
        List<MSNode> waters = getWaters();
        for (MSNode water : waters) {
          List<Atom> list = water.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          // Find the new center of mass in fractional coordinates.
          unitCell.toFractionalCoordinates(com, com);
          // Find the reciprocal translation vector to the previous COM.
          double[] frac = fractionalCoordinates[iMolecule++];
          com[0] = frac[0] - com[0];
          com[1] = frac[1] - com[1];
          com[2] = frac[2] - com[2];
          // Convert the fractional translation vector to Cartesian coordinates.
          unitCell.toCartesianCoordinates(com, com);

          double r = length(com);

          // Warn if an atom is moved more than 1 Angstrom.
          if (r > 1.0) {
            int i = iMolecule - 1;
            logger.info(String.format(" %d R: %16.8f", i, r));
            logger.info(
                String.format(" %d FRAC %16.8f %16.8f %16.8f", i, frac[0], frac[1], frac[2]));
            logger.info(String.format(" %d COM  %16.8f %16.8f %16.8f", i, com[0], com[1], com[2]));
          }

          // Move all atoms.
          for (Atom atom : list) {
            atom.move(com);
          }
        }

        // Loop over each ion
        List<MSNode> ions = getIons();
        for (MSNode ion : ions) {
          List<Atom> list = ion.getAtomList();
          // Find the center of mass
          com[0] = 0.0;
          com[1] = 0.0;
          com[2] = 0.0;
          double totalMass = 0.0;
          for (Atom atom : list) {
            double m = atom.getMass();
            com[0] += atom.getX() * m;
            com[1] += atom.getY() * m;
            com[2] += atom.getZ() * m;
            totalMass += m;
          }
          com[0] /= totalMass;
          com[1] /= totalMass;
          com[2] /= totalMass;
          // Find the new center of mass in fractional coordinates.
          unitCell.toFractionalCoordinates(com, com);
          // Find the reciprocal translation vector to the previous COM.
          double[] frac = fractionalCoordinates[iMolecule++];
          com[0] = frac[0] - com[0];
          com[1] = frac[1] - com[1];
          com[2] = frac[2] - com[2];
          // Convert the fractional translation vector to Cartesian coordinates.
          unitCell.toCartesianCoordinates(com, com);
          // Move all atoms.
          for (Atom atom : list) {
            atom.move(com);
          }
        }
        break;
      case ATOM:
        Atom[] atoms = getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
          // Convert the stored factional coordinates to Cartesian coordinates in the current
          // unitcell.
          unitCell.toCartesianCoordinates(fractionalCoordinates[i], com);
          atoms[i].moveTo(com);
        }
        break;
      case OFF:
        break;
    }
  }

  /**
   * Rotate about a point in given in the System's Local Coordinates
   *
   * @param v Vector3d
   */
  public void rotateAbout(Vector3d v) {
    Vector3d newRotPoint = new Vector3d(v);
    originToRot.getTransform(originToRotT3D);
    originToRotT3D.get(originToRotV3D);
    originToRotT3D.setTranslation(new Vector3d(0, 0, 0));
    rotToCOM.getTransform(rotToCOMT3D);
    rotToCOMT3D.get(rotToCOMV3D);
    newRotPoint.add(rotToCOMV3D);
    originToRotT3D.transform(newRotPoint);
    newRotPoint.add(originToRotV3D);
    originToRotT3D.setTranslation(newRotPoint);
    rotToCOMV3D.set(v);
    rotToCOMV3D.negate();
    rotToCOMT3D.setTranslation(rotToCOMV3D);
    originToRot.setTransform(originToRotT3D);
    rotToCOM.setTransform(rotToCOMT3D);
  }

  /**
   * sceneGraphChange
   *
   * @param newShapes a {@link java.util.List} object.
   */
  public void sceneGraphChange(List<BranchGroup> newShapes) {
    if (newShapes == null) {
      newShapes = myNewShapes;
    }

    if (newShapes.isEmpty()) {
      return;
    }

    boolean reCompile = false;
    // Check for nodes (new and/or recycled) being added to this
    // MolecularAssembly
    for (ListIterator<BranchGroup> li = newShapes.listIterator(); li.hasNext(); ) {
      BranchGroup group = li.next();
      li.remove();
      // This is code for cycling between two MolecularAssemblies
      if (group.getUserData() != null) {
        logger.info(format("%s %s", group, group.getUserData().toString()));
        /*
         * Object userData = group.getUserData(); if (userData!=this) {
         * // The appearance has already been set during a recursive
         * call to setView, // although we need to turn back on
         * Pickablility. TransformGroup tg = (TransformGroup)
         * group.getChild(0); Shape3D shape = (Shape3D) tg.getChild(0);
         * shape.setPickable(true); group.setUserData(this); if
         * (!reCompile) { if (childNodes.isLive()) {
         * childNodes.detach(); } reCompile = true; }
         * childNodes.moveTo(group);
         */
      } else {
        // This is a new group since it has no userData.
        // We can not query for the identity of its parent later, so
        // we will store it as a userData reference.
        group.setUserData(this);
        if (!reCompile) {
          if (childNodes.isLive()) {
            childNodes.detach();
          }

          reCompile = true;
        }

        childNodes.addChild(group);
      }
    }
    if (reCompile) {
      childNodes.compile();
      base.addChild(childNodes);
    }
  }

  /** {@inheritDoc} */
  @Override
  public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
    for (MSNode msNode : getAtomNodeList()) {
      MSGroup group = (MSGroup) msNode;
      group.setColor(newColorModel, color, mat);
    }

    for (MSNode m : molecules.getChildList()) {
      m.setColor(newColorModel, color, mat);
    }

    for (MSNode m : water.getChildList()) {
      m.setColor(newColorModel, color, mat);
    }

    for (MSNode m : ions.getChildList()) {
      m.setColor(newColorModel, color, mat);
    }
  }

  /**
   * setPotential
   *
   * @param potentialEnergy a {@link ffx.potential.ForceFieldEnergy} object.
   */
  public void setPotential(ForceFieldEnergy potentialEnergy) {
    this.potentialEnergy = potentialEnergy;
  }

  /** {@inheritDoc} */
  @Override
  public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
    // Just Detach the whole system branch group
    if (newViewModel == RendererCache.ViewModel.DESTROY) {
      if (switchGroup != null) {
        switchGroup.setWhichChild(Switch.CHILD_NONE);
      }
      visible = false;
    } else if (newViewModel == RendererCache.ViewModel.SHOWVRML) {
      switchGroup.setWhichChild(Switch.CHILD_ALL);
    } else if (newViewModel == RendererCache.ViewModel.HIDEVRML) {
      switchGroup.setWhichChild(0);
    } else {
      setWireWidth(RendererCache.bondwidth);
      if (newViewModel == RendererCache.ViewModel.DETAIL && childNodes.isLive()) {
        childNodes.detach();
      }
      /*
       We'll collect new Scenegraph Shapes in our newShapeNode This is to avoid the case where
       setView is called from the root node and all new shapes for every MolecularAssembly would
       then be put into the same List.
       */
      super.setView(newViewModel, myNewShapes);
      List<Molecule> moleculeList = getList(Molecule.class, new ArrayList<>());
      for (Molecule m : moleculeList) {
        m.setView(newViewModel, myNewShapes);
      }
      for (MSNode m : molecules.getChildList()) {
        m.setView(newViewModel, myNewShapes);
      }
      for (MSNode m : water.getChildList()) {
        m.setView(newViewModel, myNewShapes);
      }
      for (MSNode m : ions.getChildList()) {
        m.setView(newViewModel, myNewShapes);
      }
      if (newViewModel == RendererCache.ViewModel.INVISIBLE) {
        switchGroup.setWhichChild(0);
      }
      if (newViewModel == RendererCache.ViewModel.DETAIL) {
        childNodes.compile();
        base.addChild(childNodes);
      }
    }
  }

  /**
   * setWireWidth
   *
   * @param f a float.
   */
  public void setWireWidth(float f) {
    if (wire == null) {
      return;
    }

    lineAttributes.setLineWidth(f);
  }

  /**
   * The MolecularAssembly BranchGroup has two TransformGroups between it and the "base" node where
   * geometry is attached. If the point between the two transformations is where user rotation
   * occurs. For example, if rotating about the center of mass of the system, the RotToCOM
   * transformation will be an identity transformation (ie. none). If rotation is about some atom or
   * group of atoms within the system, then the RotToCOM transformation will be a translation from
   * that point to the COM.
   *
   * @param zero boolean
   * @return BranchGroup
   */
  private BranchGroup createScene(boolean zero) {
    originToRotT3D = new Transform3D();
    originToRotV3D = new Vector3d();
    originToRot = new TransformGroup(originToRotT3D);
    branchGroup = new BranchGroup();
    rotToCOM = new TransformGroup();
    rotToCOMT3D = new Transform3D();
    rotToCOMV3D = new Vector3d();
    // Set capabilities needed for picking and moving the MolecularAssembly
    branchGroup.setCapability(BranchGroup.ALLOW_DETACH);
    originToRot.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    originToRot.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
    originToRot.setCapability(TransformGroup.ENABLE_PICK_REPORTING);
    rotToCOM.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
    rotToCOM.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    // Put the MolecularAssembly in the middle of the scene
    if (zero) {
      originToRotV3D.set(0.0, 0.0, 0.0);
      originToRotT3D.set(originToRotV3D);
      originToRot.setTransform(originToRotT3D);
    }
    wire = renderWire();
    switchGroup = new Switch(Switch.CHILD_NONE);
    switchGroup.setCapability(Switch.ALLOW_SWITCH_WRITE);
    base = new BranchGroup();
    base.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
    base.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
    childNodes = new BranchGroup();
    childNodes.setCapability(BranchGroup.ALLOW_DETACH);
    childNodes.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
    childNodes.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
    switchGroup.addChild(base);
    if (wire != null) {
      base.addChild(wire);
    }
    vrml = loadVRML();
    if (vrml != null) {
      vrmlTG = new TransformGroup();
      vrmlTd = new Transform3D();
      vrmlTG.setTransform(vrmlTd);
      vrmlTG.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
      vrmlTG.addChild(vrml);
      switchGroup.addChild(vrmlTG);
      setView(RendererCache.ViewModel.INVISIBLE, null);
    }
    switchGroup.setWhichChild(Switch.CHILD_ALL);
    rotToCOM.addChild(switchGroup);
    originToRot.addChild(rotToCOM);
    branchGroup.addChild(originToRot);
    branchGroup.compile();
    return branchGroup;
  }

  void finishDestruction() {
    detach();
    super.destroy();
  }

  /**
   * Move mass from heavy atoms to their attached hydrogens.
   *
   * <p>The mass of each hydrogen is scaled by a factor of 3. The mass of the heavy atom is reduced
   * by 2 AMU.
   */
  private void applyHeavyHydrogen() {
    List<Bond> bonds = getBondList();
    for (Bond bond : bonds) {
      Atom a1 = bond.getAtom(0);
      Atom a2 = bond.getAtom(1);
      if (a1.isHydrogen() && a2.isHeavy()) {
        double hydrogenMass = a1.getMass();
        double heavyAtomMass = a2.getMass();
        if (hydrogenMass < 1.1) {
          a2.setMass(heavyAtomMass - 2.0 * hydrogenMass);
          a1.setMass(3.0 * hydrogenMass);
        }
      } else if (a1.isHeavy() && a2.isHydrogen()) {
        double heavyAtomMass = a1.getMass();
        double hydrogenMass = a2.getMass();
        if (hydrogenMass < 1.1) {
          a1.setMass(heavyAtomMass - 2.0 * hydrogenMass);
          a2.setMass(3.0 * hydrogenMass);
        }
      }
    }
  }

  private Atom getResidue(Atom atom, boolean create) {
    return getResidue(atom, create, Residue.ResidueType.UNK);
  }

  private Atom getResidue(Atom atom, boolean create, ResidueType defaultRT) {
    Character chainID = atom.getChainID();
    String resName = atom.getResidueName();
    int resNum = atom.getResidueNumber();
    String segID = atom.getSegID();
    // Find/Create the chain
    Polymer polymer = getPolymer(chainID, segID, create);
    if (polymer == null) {
      return null;
    }
    Residue res = polymer.getResidue(resName, resNum, create, defaultRT);
    if (create && res != null) {
      return (Atom) res.addMSNode(atom);
    }
    return null;
  }

  private Atom getMolecule(Atom atom, boolean create) {
    String resName = atom.getResidueName();
    int resNum = atom.getResidueNumber();
    Character chainID = atom.getChainID();
    String segID = atom.getSegID();

    String key = resNum + resName + segID;
    Molecule m = ionHashMap.get(key);
    if (m == null) {
      m = waterHashMap.get(key);
      if (m == null) {
        m = moleculeHashMap.get(key);
      }
    }
    if (m != null) {
      return (Atom) m.addMSNode(atom);
    }

    if (create) {
      boolean isWater = false;
      boolean isIon = false;
      if (StringUtils.looksLikeWater(resName)) {
        resName = StringUtils.STANDARD_WATER_NAME;
        isWater = true;
      } else if (StringUtils.looksLikeIon(resName)) {
        resName = StringUtils.tryParseIon(resName);
        isIon = true;
      }
      atom.setResName(resName);
      m = new Molecule(resName, resNum, chainID, segID);
      m.addMSNode(atom);
      if (resName == null) {
        logger.warning(
            format(
                " Attempting to create a molecule %s with a null name on atom %s! Defaulting to creating a generic Molecule.",
                m, atom));
        molecules.add(m);
        moleculeHashMap.put(key, m);
      } else if (isWater) {
        water.add(m);
        waterHashMap.put(key, m);
      } else if (isIon) {
        ions.add(m);
        ionHashMap.put(key, m);
      } else {
        molecules.add(m);
        moleculeHashMap.put(key, m);
      }
      return atom;
    } else {
      return null;
    }
  }

  private void recurseVRML(Node node) {
    if (node instanceof Shape3D) {
      Shape3D s3d = (Shape3D) node;
      PickTool.setCapabilities(s3d, PickTool.INTERSECT_COORD);
    } else if (node instanceof SharedGroup) {
      SharedGroup sg = (SharedGroup) node;
      for (Iterator<Node> e = sg.getAllChildren(); e.hasNext(); ) {
        recurseVRML(e.next());
      }
    } else if (node instanceof BranchGroup) {
      BranchGroup bg = (BranchGroup) node;
      for (Iterator<Node> e = bg.getAllChildren(); e.hasNext(); ) {
        recurseVRML(e.next());
      }
    } else if (node instanceof TransformGroup) {
      TransformGroup vrmlTG1 = (TransformGroup) node;
      for (Iterator<Node> e = vrmlTG1.getAllChildren(); e.hasNext(); ) {
        node = e.next();
        recurseVRML(node);
      }
    } else if (node instanceof Link) {
      Link link = (Link) node;
      recurseVRML(link.getSharedGroup());
    } else if (node instanceof Group) {
      Group group = (Group) node;
      for (Iterator<Node> e = group.getAllChildren(); e.hasNext(); ) {
        recurseVRML(e.next());
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  protected void removeLeaves() {
    super.removeLeaves();
    MSNode macroNode = getAtomNode();
    if (macroNode != null) {
      if (macroNode.getChildCount() > 0) {
        getAtomNode().setName("Macromolecules " + "(" + macroNode.getChildCount() + ")");
      } else if (macroNode.getParent() == this) {
        removeChild(macroNode);
      }
    }

    if (molecules.getChildCount() == 0) {
      removeChild(molecules);
    } else {
      molecules.setName("Hetero Molecules " + "(" + molecules.getChildCount() + ")");
    }

    if (ions.getChildCount() == 0) {
      removeChild(ions);
    } else {
      ions.setName("Ions " + "(" + ions.getChildCount() + ")");
    }

    if (water.getChildCount() == 0) {
      removeChild(water);
    } else {
      water.setName("Water " + "(" + water.getChildCount() + ")");
    }
  }

  private Shape3D renderWire() {
    List<Bond> bonds = getBondList();
    int numbonds = bonds.size();
    if (numbonds < 1) {
      return null;
    }

    Vector3d bondmidpoint = new Vector3d();
    double[] mid = {0, 0, 0};
    Vector3d v1 = new Vector3d();
    Vector3d v2 = new Vector3d();
    float[] a1 = {0, 0, 0};
    float[] a2 = {0, 0, 0};
    float[] col = new float[4];

    Atom atom1, atom2;
    LineArray la =
        new LineArray(
            4 * numbonds,
            GeometryArray.COORDINATES | GeometryArray.COLOR_4 | GeometryArray.NORMALS);
    la.setCapability(LineArray.ALLOW_COORDINATE_WRITE);
    la.setCapability(LineArray.ALLOW_COORDINATE_READ);
    la.setCapability(LineArray.ALLOW_COLOR_WRITE);
    la.setCapability(LineArray.ALLOW_COUNT_READ);
    la.setCapability(LineArray.ALLOW_INTERSECT);
    la.setCapability(LineArray.ALLOW_FORMAT_READ);
    atomLookUp = new Atom[4 * numbonds];
    int i = 0;
    col[3] = 0.9f;
    for (Bond bond : bonds) {
      bond.setWire(la, i);
      atom1 = bond.getAtom(0);
      atom2 = bond.getAtom(1);
      atom1.getV3D(v1);
      atom2.getV3D(v2);
      a1[0] = (float) v1.x;
      a1[1] = (float) v1.y;
      a1[2] = (float) v1.z;
      a2[0] = (float) v2.x;
      a2[1] = (float) v2.y;
      a2[2] = (float) v2.z;
      // Find the bond center
      bondmidpoint.add(v1, v2);
      bondmidpoint.scale(0.5d);
      bondmidpoint.get(mid);

      // Atom #1
      Atom.AtomColor.get(atom1.getAtomicNumber()).get(col);
      atomLookUp[i] = atom1;
      la.setCoordinate(i, a1);
      la.setColor(i, col);
      la.setNormal(i, a2);
      i++;

      atomLookUp[i] = atom1;
      la.setCoordinate(i, mid);
      la.setColor(i, col);
      la.setNormal(i, a2);
      i++;

      // Atom #2
      Atom.AtomColor.get(atom2.getAtomicNumber()).get(col);
      atomLookUp[i] = atom2;
      la.setCoordinate(i, a2);
      la.setColor(i, col);
      la.setNormal(i, a1);
      i++;

      atomLookUp[i] = atom2;
      la.setCoordinate(i, mid);
      la.setColor(i, col);
      la.setNormal(i, a1);
      i++;
    }

    ColoringAttributes cola =
        new ColoringAttributes(new Color3f(), ColoringAttributes.SHADE_GOURAUD);
    Appearance app = new Appearance();
    lineAttributes = new LineAttributes();
    lineAttributes.setLineWidth(RendererCache.bondwidth);
    lineAttributes.setCapability(LineAttributes.ALLOW_WIDTH_WRITE);
    lineAttributes.setLineAntialiasingEnable(true);
    app.setLineAttributes(lineAttributes);
    app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_READ);
    app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_WRITE);
    RenderingAttributes ra = new RenderingAttributes();
    ra.setAlphaTestValue(0.1f);
    ra.setAlphaTestFunction(RenderingAttributes.GREATER);
    ra.setDepthBufferEnable(true);
    ra.setDepthBufferWriteEnable(true);
    app.setRenderingAttributes(ra);
    app.setColoringAttributes(cola);
    Shape3D wireframe = new Shape3D(la, app);
    // PickTool.setCapabilities(wire, PickTool.INTERSECT_COORD);
    wireframe.setUserData(this);
    wireframe.setBounds(new BoundingSphere(new Point3d(0, 0, 0), 1000.0));
    try {
      wireframe.setBoundsAutoCompute(false);
    } catch (Exception e) {
      e.printStackTrace();
    }

    wireframe.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
    wireframe.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
    wireframe.setCapability(Shape3D.ALLOW_LOCAL_TO_VWORLD_READ);
    return wireframe;
  }

  /**
   * Returns a list of all water molecules in the system, including both those properly under the
   * water node and those under the molecules node.
   *
   * @return All waters in the system
   */
  private List<Molecule> getAllWatersInclMistyped() {
    Stream<Molecule> mistyped =
        getMolecules()
            .parallelStream()
            .map((MSNode m) -> (Molecule) m)
            .filter(
                (Molecule m) -> {
                  List<Atom> atoms = m.getAtomList();
                  if (atoms.size() != 3) {
                    return false;
                  }
                  int nO = 0;
                  int nH = 0;
                  for (Atom atom : atoms) {
                    int el = atom.getAtomicNumber();
                    if (el == 1) {
                      ++nH;
                    } else if (nH == 8) {
                      ++nO;
                    }
                  }
                  return nO == 1 && nH == 2;
                });
    return Stream.concat(mistyped, getWaters().stream().map((MSNode m) -> (Molecule) m))
        .distinct()
        .collect(Collectors.toList());
  }

  /** Renames water protons to H1 and H2. */
  void renameWaterProtons() {
    for (Molecule water : getAllWatersInclMistyped()) {
      Atom H1 = water.getAtomByName("H1", false);
      Atom H2 = water.getAtomByName("H2", false);
      if (H1 != null && H2 != null) {
        continue;
      }
      for (Atom a : water.getAtomList()) {
        if (a.getAtomicNumber() == 1) {
          if (H1 == null) {
            H1 = a;
            H1.setName("H1");
          } else if (H2 == null) {
            H2 = a;
            H2.setName("H2");
          }
        }
      }
    }
  }

  private void moveIntoUnitCell(MSNode[] groups) {
    if (groups != null && groups.length > 0) {
      moveIntoUnitCell(Arrays.asList(groups));
    }
  }

  /**
   * Move the center of each listed chemical entity into the primary unit cell.
   *
   * @param groups Move each MSNode into the primary unit cell.
   */
  private void moveIntoUnitCell(List<MSNode> groups) {
    Crystal cryst = getCrystal().getUnitCell();
    if (cryst.aperiodic()) {
      return;
    }
    for (MSNode group : groups) {
      double[] com = new double[3];
      double[] xyz = new double[3];
      double[] translate = new double[3];
      List<Atom> atoms = group.getAtomList();
      double totMass = 0;
      for (Atom atom : atoms) {
        double mass = atom.getMass();
        totMass += mass;
        xyz = atom.getXYZ(xyz);
        com[0] += mass * xyz[0];
        com[1] += mass * xyz[1];
        com[2] += mass * xyz[2];
      }
      com[0] /= totMass;
      com[1] /= totMass;
      com[2] /= totMass;

      // Move the COM to the primary unit cell
      cryst.toPrimaryCell(com, translate);

      // The translation vector is difference between the new location and the current COM.
      translate[0] -= com[0];
      translate[1] -= com[1];
      translate[2] -= com[2];

      for (Atom atom : atoms) {
        atom.move(translate);
      }
    }
  }

  public enum FractionalMode {
    OFF,
    MOLECULE,
    ATOM
  }
}
