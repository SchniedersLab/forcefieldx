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
package ffx.potential.extended;

import static ffx.potential.extended.ExtConstants.RNG;
import static ffx.potential.extended.TitrationUtils.isTitratableHydrogen;
import static ffx.potential.parameters.MultipoleType.zeroM;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.reduction.SharedDouble;
import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Descriptions;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtendedSystem.ExtendedSystemConfig;
import ffx.potential.parameters.MultipoleType;
import ffx.utilities.Constants;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Logger;

/**
 * A generalized extended system variable. Treatment of ESVs: a. Bonded terms interpolate linearly
 * between end states. ESVs based on MultiResidue (e.g. TitrationESV) place a multiplier in the term
 * objects themselves. b. PME and vdW scaling and derivatives are handled inside these classes'
 * inner loops. Softcoring follows the same form used by OST lambda.
 *
 * @author Stephen LuCore
 * @since 1.0
 */
public abstract class ExtendedVariable {

  /** Constant <code>logger</code> */
  protected static final Logger logger = Logger.getLogger(ExtendedVariable.class.getName());

  public final int esvIndex;
  /** Magnitude of the discretization bias in kcal/mol. */
  private final double discrBiasMag;
  /**
   * Sigmoidal switching function. Maps lambda -> S(lambda) which has a flatter deriv near
   * zero/unity. Properties: {S(0)=0, S(1)=1, dS(0)=0, dS(1)=0, S(1-L)=1-S(L)}.
   */
  private final MultiplicativeSwitch switchingFunction;
  /* Atom lists and scaled terms                  */
  private final Residue residueForeground; // resOne*lamedh + resZero*(1-lamedh)
  private final Residue residueBackground;
  private final List<Atom> backbone;
  private final List<Atom> atomsForeground; // (L); side-chain only
  private final List<Atom>
      atomsBackground; // (1-L); side-chain only; permanently disconnected from assembly
  private final List<Atom> atomsShared; // all foreground atoms except titrating hydrogens
  private final List<Atom> atomsUnshared; // titrating (and thus foreground) atoms
  private final int moleculeNumber;
  private final SharedDouble bondedDeriv = new SharedDouble(); // bonded dUdL reduction target
  private final HashMap<Class<? extends BondedTerm>, SharedDouble>
      fgBondedDerivDecomp; // foreground dUdL by term
  private final HashMap<Class<? extends BondedTerm>, SharedDouble>
      bgBondedDerivDecomp; // background dUdL by term
  private final HashMap<Atom, Atom> fg2bg =
      new HashMap<>(); // maps multipole end points of this ESV's lambda path
  private final ExtendedSystemConfig config;
  /* Lambda and derivative variables. */
  private double lambda = 1.0; // ESVs travel on {0,1}
  private double theta; // Propagates lambda particle via "lambda=sin(theta)^2"
  private double halfThetaVelocity = 0.0; // from OST, start theta with zero velocity
  private StringBuilder SB = new StringBuilder();
  /** Discretization bias and its (chain rule) derivative. */
  private double discrBias, dDiscrBiasdL;
  /** Switched lambda value and its (chain rule) derivative. */
  private double lSwitch, dlSwitch;
  /* Bonded energy and derivative handling        */
  private List<BondedTerm> bondedFg,
      bondedBg; // valence terms for each side; mola won't see zro by default
  private MSNode termNode; // modified to contain all applicable bonded terms

  /**
   * Prefer ExtendedSystem::populate to manual ESV creation.
   *
   * @param multiRes: from TitrationUtils.titrationFactory()
   * @param initialLambda: (optional) starting position of the extended particle
   * @param esvSystem a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public ExtendedVariable(ExtendedSystem esvSystem, MultiResidue multiRes, double initialLambda) {
    this.esvIndex = esvSystem.requestIndexing();
    this.config = esvSystem.config;
    this.discrBiasMag = config.discrBias;
    this.switchingFunction = new MultiplicativeSwitch(1.0, 0.0);
    setInitialLambda(lambda);

    residueForeground = multiRes.getActive();
    termNode = residueForeground.getTermNode();
    residueBackground = multiRes.getInactive().get(0);
    moleculeNumber = residueForeground.getAtomList().get(0).getMoleculeNumber();

    backbone = new ArrayList<>();
    atomsForeground = new ArrayList<>();
    atomsBackground = new ArrayList<>();
    atomsShared = new ArrayList<>();
    atomsUnshared = new ArrayList<>();

    if (config.decomposeBonded) {
      fgBondedDerivDecomp = new HashMap<>();
      bgBondedDerivDecomp = new HashMap<>();
    } else {
      fgBondedDerivDecomp = null;
      bgBondedDerivDecomp = null;
    }

    // Fill the atom lists.
    for (String bbName : ExtConstants.backboneNames) {
      Atom bb = (Atom) residueForeground.getAtomNode(bbName);
      if (bb != null) {
        backbone.add(bb);
      }
    }
    for (Atom fg : residueForeground.getAtomList()) {
      if (!backbone.contains(fg)) {
        atomsForeground.add(fg);
        Atom bg = residueBackground.getAtomByName(fg.getName(), true);
        if (bg == null) {
          atomsUnshared.add(fg);
          /* The following check ought to be safely removable if you've
           * defined ExtendedVariables that are not TitrationESVs.      */
          assert (isTitratableHydrogen(fg));
          if (!isTitratableHydrogen(fg)) {
            logger.warning(
                format(
                    "ExtendedVariable could not identify a companion for foreground atom %s.", fg));
            throw new IllegalStateException();
          }
        } else {
          atomsShared.add(fg);
          fg2bg.put(fg, bg);
        }
      }
    }
    for (Atom a0 : residueBackground.getAtomList()) {
      if (!backbone.contains(a0)) {
        assert (!atomsForeground.contains(a0));
        assert (!isTitratableHydrogen(a0));
        if (atomsForeground.contains(a0) || isTitratableHydrogen(a0)) {
          logger.warning("Error: inappropriate background atom.");
          throw new IllegalStateException();
        }
        atomsBackground.add(a0);
        a0.setBackground();
      }
    }

    /* Assign foreground atom indices to their corresponding background atoms. */
    if (config.cloneXyzIndices) {
      for (Atom fg : fg2bg.keySet()) {
        fg2bg.get(fg).setXyzIndex(fg.getXyzIndex());
      }
    }

    // Fill bonded term list and set all esvLambda values.
    bondedFg = residueForeground.getDescendants(BondedTerm.class);
    bondedBg = residueBackground.getDescendants(BondedTerm.class);
    if (config.bonded) {
      MSNode extendedTermNode = new MSNode(format("Extended (%d)", bondedBg.size()));
      for (MSNode node : residueBackground.getTermNode().getChildList()) {
        extendedTermNode.add(node);
      }
      multiRes.getActive().getTermNode().add(extendedTermNode);
      updateBondedLambdas();
    }

    /* Background atoms don't get automatically typed by PME since they're
     * disconnected from the molecular assembly; must be done manually. */
    if (config.electrostatics) {
      esvSystem.initializeBackgroundMultipoles(atomsBackground);
      updateMultipoleTypes();
    }

    describe();
  }

  /** List all the atoms and bonded terms associated with each end state. */
  public final void describe() {
    SB.append(format(" %s\n", this.toString()));
    SB.append(format("   %-50s %-50s\n", "Shared Atoms", "(Background)"));
    for (int i = 0; i < atomsShared.size(); i++) {
      Atom ai = atomsShared.get(i);
      SB.append(
          format(
              "     %-50s %-50s\n",
              ai.describe(Atom.Descriptions.Default).trim(),
              fg2bg.get(ai).describe(Atom.Descriptions.Trim)));
    }
    SB.append(format("   Unshared Atoms\n"));
    for (Atom atom : atomsUnshared) {
      SB.append(format("%s\n", atom));
    }
    SB.append(format("   %-50s %-50s\n", "Bonded Terms", "(Background)"));
    MSNode extendedNode =
        termNode.getChildList().stream()
            .filter(node -> node.toString().contains("Extended"))
            .findAny()
            .orElse(null);
    for (MSNode term : termNode.getChildList()) {
      if (term == extendedNode) {
        continue;
      }
      MSNode background =
          (extendedNode == null)
              ? null
              : extendedNode.getChildList().stream()
                  .filter(
                      node ->
                          node.toString()
                              .startsWith(
                                  term.toString().substring(0, term.toString().indexOf("("))))
                  .findAny()
                  .orElse(null);
      String bgTermString = (background != null) ? background.toString() : "";
      SB.append(format("     %-50s %-50s\n", term.toString(), bgTermString));
    }
    logger.info(SB.toString());
  }

  /**
   * Getter for the field <code>esvIndex</code>.
   *
   * @return a int.
   */
  public final int getEsvIndex() {
    return esvIndex;
  }

  /**
   * The unswitched lambda value, ie input to S(L).
   *
   * @return a double.
   */
  public final double getLambda() {
    return lambda;
  }

  /**
   * Should only be called by ExtendedSystem since an updateListeners() call is required afterwards
   * to notify VdW and PME.
   *
   * @param lambda a double.
   */
  protected void setLambda(double lambda) {
    setLambda(lambda, true);
  }

  /**
   * getName.
   *
   * @return a {@link java.lang.String} object.
   */
  public String getName() {
    return String.format("Esv%d", esvIndex);
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return String.format("%s (%4.2f)", this.getName(), getLambda());
  }

  /**
   * Should include at least the discretization bias; add any type-specific biases (eg pH).
   *
   * @param temperature a double.
   * @param print a boolean.
   * @return a double.
   */
  protected abstract double getTotalBias(double temperature, boolean print);

  /**
   * Should include at least the discretization bias; add any type-specific biases (eg pH).
   *
   * @param temperature a double.
   * @param print a boolean.
   * @return a double.
   */
  protected abstract double getTotalBiasDeriv(double temperature, boolean print);

  /**
   * Propagate lambda using Langevin dynamics. Check that temperature goes to the value used below
   * (when set as a constant), even when sim is decoupled. Be sure it call setLambda() rather than
   * using direct access for array resizing, etc.
   *
   * @param dEdEsv a double.
   * @param dt a double.
   * @param setTemperature a double.
   */
  protected void propagate(double dEdEsv, double dt, double setTemperature) {
    if (!config.propagation) {
      return;
    }
    double rt2 = 2.0 * Constants.R * setTemperature * config.thetaFriction / dt;
    double randomForce = sqrt(rt2) * RNG.nextGaussian() / ExtConstants.forceToKcal;
    double dEdL = -dEdEsv * sin(2.0 * theta);
    halfThetaVelocity =
        (halfThetaVelocity * (2.0 * config.thetaMass - config.thetaFriction * dt)
                + ExtConstants.forceToKcalSquared * 2.0 * dt * (dEdL + randomForce))
            / (2.0 * config.thetaMass + config.thetaFriction * dt);
    theta = theta + dt * halfThetaVelocity;

    if (theta > PI) {
      theta -= 2.0 * PI;
    } else if (theta <= -PI) {
      theta += 2.0 * PI;
    }

    double sinTheta = sin(theta);
    setLambda(sinTheta * sinTheta);
  }

  private void setLambda(double lambda, boolean updateComponents) {
    this.lambda = lambda;
    this.lSwitch = (config.allowLambdaSwitch) ? switchingFunction.taper(lambda) : lambda;
    this.dlSwitch = (config.allowLambdaSwitch) ? switchingFunction.dtaper(lambda) : 1.0;
    theta = Math.asin(Math.sqrt(lambda));
    discrBias = discrBiasMag - 4 * discrBiasMag * (lambda - 0.5) * (lambda - 0.5);
    dDiscrBiasdL = -8 * discrBiasMag * (lambda - 0.5);
    if (updateComponents) {
      updateMultipoleTypes();
      updateBondedLambdas();
    }
  }

  private void setInitialLambda(double lambda) {
    setLambda(lambda, false);
  }

  /** Invoked by ExtendedSystem after lambda changes and by PME after multipole rotation. */
  protected final void updateMultipoleTypes() {
    if (!config.electrostatics) {
      return;
    }
    /* If not softcoring unshared atoms, then scale them as well
     * (with an implied zero-multipole background atom).	  */
    List<Atom> iterate = ExtUtils.joinedListView(atomsShared, atomsUnshared);
    for (Atom fg : iterate) {
      //            MultipoleType Ptype, Utype;
      Atom bg = fg2bg.get(fg);
      final MultipoleType Ptype = fg.getMultipoleType();
      final MultipoleType Utype;
      if (bg == null) {
        if (atomsUnshared.contains(fg)) {
          Utype = new MultipoleType(zeroM, Ptype.frameAtomTypes, Ptype.frameDefinition, false);
        } else {
          logger.warning("Error @ESV.updateMultipoles: bg null && !unshared.");
          Utype = null;
        }
      } else {
        Utype = bg.getMultipoleType();
      }
      MultipoleType types[] = new MultipoleType[] {Ptype, Utype};
      double mWeights[], mdotWeights[];
      if (config.allowLambdaSwitch && config.nonlinearMultipoles) {
        mWeights = new double[] {lSwitch, 1.0 - lSwitch};
        mdotWeights = new double[] {dlSwitch, -dlSwitch};
      } else {
        mWeights = new double[] {lambda, 1.0 - lambda};
        mdotWeights = new double[] {1.0, -1.0};
      }
      if (Ptype == null) {
        // Multipoles not yet defined.
        continue;
      }

      if (Utype == null) {
        SB.append(format("Error @ESV.updateMultipoleTypes: bgType null."));
        SB.append(format("   fg: %s, '%s'", fg.toString(), fg.getName()));
        SB.append(format(" Background atoms available for match: "));
        for (Atom debug : atomsBackground) {
          SB.append(format("   bg: %s, '%s'", debug.toString(), debug.getName()));
        }
        logger.warning(SB.toString());
        continue;
      }
      final double[] esvMultipole;
      final double[] esvMultipoleDot;
      final MultipoleType esvType;
      final MultipoleType esvDotType;
      try {
        esvType = MultipoleType.weightMultipoleTypes(types, mWeights, Ptype.frameAtomTypes);
        esvDotType = MultipoleType.weightMultipoleTypes(types, mdotWeights, Ptype.frameAtomTypes);
      } catch (IllegalArgumentException ex) {
        logger.warning(
            format(
                "Multipole scaling failed for fg,bg atoms: %s %s",
                fg.describe(Descriptions.Trim), bg.describe(Descriptions.Trim)));
        throw ex;
      }
      if (esvType == null || esvDotType == null) {
        logger.severe("Error @ESV.updateMultipoleTypes: M or Mdot null.");
      }
      if (isTitratableHydrogen(fg)) {
        final double scaledAlpha = fg.getPolarizeType().polarizability * getLambda();
        fg.setEsv(this, esvType, esvDotType, scaledAlpha);
      } else {
        fg.setEsv(this, esvType, esvDotType);
      }
      SB.append(
          format(
              " Assigning ESV MultipoleTypes for atom %s\n",
              fg.describe(Atom.Descriptions.Resnum_Name)));
      SB.append(format("  U: %.2f*%s\n", lambda, Ptype.toCompactBohrString()));
      SB.append(format("  P: %.2f*%s\n", 1.0 - lambda, Utype.toCompactBohrString()));
      SB.append(format("  M:      %s\n", fg.getEsvMultipole().toCompactBohrString()));
      SB.append(format("  Mdot:   %s\n", fg.getEsvMultipoleDot().toCompactBohrString()));
      if (ExtUtils.verbose) {
        logger.info(SB.toString());
      }
    }
  }

  /**
   * getLambdaSwitch.
   *
   * @return a double.
   */
  protected final double getLambdaSwitch() {
    return (config.allowLambdaSwitch) ? lSwitch : lambda; // S(L)
  }

  /**
   * getSwitchDeriv.
   *
   * @return a double.
   */
  protected final double getSwitchDeriv() {
    return (config.allowLambdaSwitch) ? dlSwitch : 1.0; // dS(L)dL
  }

  /**
   * From Shen and Huang 2016; drives ESVs to zero/unity. bias = 4B*(L-0.5)^2
   *
   * @return a double.
   */
  protected double getDiscrBias() {
    return discrBias;
  }

  /**
   * dBiasdL = -8B*(L-0.5)
   *
   * @return a double.
   */
  protected double getDiscrBiasDeriv() {
    return dDiscrBiasdL;
  }

  /**
   * getBackgroundForAtom.
   *
   * @param foreground a {@link ffx.potential.bonded.Atom} object.
   * @return a {@link ffx.potential.bonded.Atom} object.
   */
  protected Atom getBackgroundForAtom(Atom foreground) {
    return fg2bg.get(foreground);
  }

  /**
   * viewUnsharedAtoms.
   *
   * @return a {@link java.util.List} object.
   */
  protected List<Atom> viewUnsharedAtoms() {
    return Collections.unmodifiableList(atomsUnshared);
  }

  /**
   * viewSharedAtoms.
   *
   * @return a {@link java.util.List} object.
   */
  protected List<Atom> viewSharedAtoms() {
    return Collections.unmodifiableList(atomsShared);
  }

  /**
   * viewBackgroundAtoms.
   *
   * @return a {@link java.util.List} object.
   */
  protected List<Atom> viewBackgroundAtoms() {
    return Collections.unmodifiableList(atomsBackground);
  }

  private void updateBondedLambdas() {
    if (!config.bonded) {
      return;
    }
    bondedDeriv.set(0.0);
    double Sl = getLambdaSwitch();
    double dSldL = getSwitchDeriv();
    for (BondedTerm bt1 : bondedFg) {
      if (config.decomposeBonded) {
        bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv, fgBondedDerivDecomp);
        fgBondedDerivDecomp.clear();
      } else {
        bt1.attachExtendedVariable(Sl, dSldL, bondedDeriv);
      }
    }
    for (BondedTerm bt0 : bondedBg) {
      if (config.decomposeBonded) {
        bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv, bgBondedDerivDecomp);
        bgBondedDerivDecomp.clear();
      } else {
        bt0.attachExtendedVariable(1.0 - Sl, -dSldL, bondedDeriv);
      }
    }
  }

  /**
   * Getter for the field <code>bondedDeriv</code>.
   *
   * @return a double.
   */
  protected double getBondedDeriv() {
    return bondedDeriv.get();
  }

  /**
   * getBondedDerivDecomp.
   *
   * @return a {@link java.util.HashMap} object.
   */
  protected HashMap<Class<? extends BondedTerm>, SharedDouble> getBondedDerivDecomp() {
    return fgBondedDerivDecomp;
  }

  /**
   * getBackgroundBondedDerivDecomp.
   *
   * @return a {@link java.util.HashMap} object.
   */
  protected HashMap<Class<? extends BondedTerm>, SharedDouble> getBackgroundBondedDerivDecomp() {
    return bgBondedDerivDecomp;
  }
}
