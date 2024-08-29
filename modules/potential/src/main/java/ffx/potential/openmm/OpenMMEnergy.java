// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.potential.openmm;

import edu.rit.mp.CharacterBuf;
import edu.rit.pj.Comm;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Platform;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import javax.annotation.Nullable;
import java.io.File;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static java.lang.Double.isFinite;
import static java.lang.String.format;

/**
 * Compute the potential energy and derivatives using OpenMM.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("deprecation")
public class OpenMMEnergy extends ForceFieldEnergy {

  private static final Logger logger = Logger.getLogger(OpenMMEnergy.class.getName());

  /**
   * FFX Platform.
   */
  private final Platform platform;
  /**
   * OpenMM Context.
   */
  private OpenMMContext openMMContext;
  /**
   * OpenMM System.
   */
  private OpenMMSystem openMMSystem;
  /**
   * The atoms this ForceFieldEnergyOpenMM operates on.
   */
  private final Atom[] atoms;
  /**
   * Truncate the normal OpenMM Lambda Path from 0 ... 1 to Lambda_Start ... 1. This is useful for
   * conformational optimization if full removal of vdW interactions is not desired (i.e. lambdaStart
   * = ~0.2).
   */
  private double lambdaStart;
  /**
   * Use two-sided finite difference dU/dL.
   */
  private boolean twoSidedFiniteDifference = true;
  /**
   * Lambda step size for finite difference dU/dL.
   */
  private final double finiteDifferenceStepSize;

  /**
   * ForceFieldEnergyOpenMM constructor; offloads heavy-duty computation to an OpenMM Platform while
   * keeping track of information locally.
   *
   * @param molecularAssembly Assembly to construct energy for.
   * @param requestedPlatform requested OpenMM platform to be used.
   * @param nThreads          Number of threads to use in the super class ForceFieldEnergy instance.
   */
  public OpenMMEnergy(MolecularAssembly molecularAssembly, Platform requestedPlatform, int nThreads) {
    super(molecularAssembly, nThreads);

    Crystal crystal = getCrystal();
    int symOps = crystal.spaceGroup.getNumberOfSymOps();
    if (symOps > 1) {
      logger.info("");
      logger.severe(" OpenMM does not support symmetry operators.");
    }

    logger.info("\n Initializing OpenMM");

    ForceField forceField = molecularAssembly.getForceField();
    atoms = molecularAssembly.getAtomArray();

    // Load the OpenMM plugins
    this.platform = requestedPlatform;
    ffx.openmm.Platform openMMPlatform = OpenMMContext.loadPlatform(platform, forceField);

    // Create the OpenMM System.
    openMMSystem = new OpenMMSystem(this);
    openMMSystem.addForces();

    // Create the Context.
    openMMContext = new OpenMMContext(openMMPlatform, openMMSystem, atoms);

    // Expand the path [lambda-start .. 1.0] to the interval [0.0 .. 1.0].
    lambdaStart = forceField.getDouble("LAMBDA_START", 0.0);
    if (lambdaStart > 1.0) {
      lambdaStart = 1.0;
    } else if (lambdaStart < 0.0) {
      lambdaStart = 0.0;
    }

    finiteDifferenceStepSize = forceField.getDouble("FD_DLAMBDA", 0.001);
    twoSidedFiniteDifference = forceField.getBoolean("FD_TWO_SIDED", twoSidedFiniteDifference);
  }

  /**
   * Gets the default coprocessor device, ignoring any CUDA_DEVICE over-ride. This is either
   * determined by process rank and the availableDevices/CUDA_DEVICES property, or just 0 if neither
   * property is sets.
   *
   * @param props Properties in use.
   * @return Pre-override device index.
   */
  public static int getDefaultDevice(CompositeConfiguration props) {
    String availDeviceProp = props.getString("availableDevices", props.getString("CUDA_DEVICES"));
    if (availDeviceProp == null) {
      int nDevs = props.getInt("numCudaDevices", 1);
      availDeviceProp = IntStream.range(0, nDevs).mapToObj(Integer::toString)
          .collect(Collectors.joining(" "));
    }
    availDeviceProp = availDeviceProp.trim();

    String[] availDevices = availDeviceProp.split("\\s+");
    int nDevs = availDevices.length;
    int[] devs = new int[nDevs];
    for (int i = 0; i < nDevs; i++) {
      devs[i] = Integer.parseInt(availDevices[i]);
    }

    logger.info(format(" Available devices: %d.", nDevs));

    int index = 0;
    try {
      Comm world = Comm.world();
      if (world != null) {
        int size = world.size();

        // Format the host as a CharacterBuf of length 100.
        int messageLen = 100;
        String host = world.host();
        // Truncate to max 100 characters.
        host = host.substring(0, Math.min(messageLen, host.length()));
        // Pad to 100 characters.
        host = format("%-100s", host);
        char[] messageOut = host.toCharArray();
        CharacterBuf out = CharacterBuf.buffer(messageOut);

        // Now create CharacterBuf array for all incoming messages.
        char[][] incoming = new char[size][messageLen];
        CharacterBuf[] in = new CharacterBuf[size];
        for (int i = 0; i < size; i++) {
          in[i] = CharacterBuf.buffer(incoming[i]);
        }

        try {
          world.allGather(out, in);
        } catch (IOException ex) {
          logger.severe(format(" Failure at the allGather step for determining rank: %s\n%s", ex, Utilities.stackTraceToString(ex)));
        }
        int ownIndex = -1;
        int rank = world.rank();
        boolean selfFound = false;

        for (int i = 0; i < size; i++) {
          String hostI = new String(incoming[i]);
          if (hostI.equalsIgnoreCase(host)) {
            ++ownIndex;
            if (i == rank) {
              selfFound = true;
              break;
            }
          }
        }
        if (!selfFound) {
          logger.severe(format(" Rank %d: Could not find any incoming host messages matching self %s!", rank, host.trim()));
        } else {
          index = ownIndex % nDevs;
        }
      }
    } catch (IllegalStateException ise) {
      // Behavior is just to keep index = 0.
    }
    return devs[index];
  }

  /**
   * Create an OpenMM Context.
   *
   * <p>Context.free() must be called to free OpenMM memory.
   *
   * @param integratorName Integrator to use.
   * @param timeStep       Time step.
   * @param temperature    Temperature (K).
   * @param forceCreation  Force a new Context to be created, even if the existing one matches the
   *                       request.
   */
  public void updateContext(String integratorName, double timeStep, double temperature, boolean forceCreation) {
    openMMContext.update(integratorName, timeStep, temperature, forceCreation);
  }

  /**
   * Create an immutable OpenMM State.
   *
   * <p>State.free() must be called to free OpenMM memory.
   *
   * @param mask The State mask.
   * @return Returns the State.
   */
  public OpenMMState getOpenMMState(int mask) {
    return openMMContext.getOpenMMState(mask);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    boolean ffxFFEDestroy = super.destroy();
    free();
    logger.fine(" Destroyed the Context, Integrator, and OpenMMSystem.");
    return ffxFFEDestroy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x, boolean verbose) {

    if (lambdaBondedTerms) {
      return 0.0;
    }

    // Make sure the context has been created.
    openMMContext.update();

    updateParameters(atoms);

    // Unscale the coordinates.
    unscaleCoordinates(x);

    setCoordinates(x);

    OpenMMState openMMState = openMMContext.getOpenMMState(OpenMM_State_Energy);
    double e = openMMState.potentialEnergy;
    openMMState.destroy();

    if (!isFinite(e)) {
      String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      if (lambdaTerm) {
        openMMSystem.printLambdaValues();
      }
      throw new EnergyException(message);
    }

    if (verbose) {
      logger.log(Level.INFO, String.format("\n OpenMM Energy: %14.10g", e));
    }

    // Rescale the coordinates.
    scaleCoordinates(x);

    return e;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    return energyAndGradient(x, g, false);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g, boolean verbose) {
    if (lambdaBondedTerms) {
      return 0.0;
    }

    // Un-scale the coordinates.
    unscaleCoordinates(x);

    // Make sure a context has been created.
    openMMContext.update();

    // long time = -System.nanoTime();
    setCoordinates(x);
    // time += System.nanoTime();
    // logger.info(format(" Load coordinates time %10.6f (sec)", time * 1.0e-9));

    // time = -System.nanoTime();
    OpenMMState openMMState = openMMContext.getOpenMMState(OpenMM_State_Energy | OpenMM_State_Forces);
    double e = openMMState.potentialEnergy;
    g = openMMState.getGradient(g);
    openMMState.destroy();
    // time += System.nanoTime();
    // logger.info(format(" Calculate energy time %10.6f (sec)", time * 1.0e-9));

    if (!isFinite(e)) {
      String message = format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      if (lambdaTerm) {
        openMMSystem.printLambdaValues();
      }
      throw new EnergyException(message);
    }

    // if (vdwLambdaTerm) {
    //    PointerByReference parameterArray = OpenMM_State_getEnergyParameterDerivatives(state);
    //    int numDerives = OpenMM_ParameterArray_getSize(parameterArray);
    //    if (numDerives > 0) {
    //        double vdwdUdL = OpenMM_ParameterArray_get(parameterArray,
    // pointerForString("vdw_lambda")) / OpenMM_KJPerKcal;
    //    }
    // }

    if (maxDebugGradient < Double.MAX_VALUE) {
      boolean extremeGrad = Arrays.stream(g)
          .anyMatch((double gi) -> (gi > maxDebugGradient || gi < -maxDebugGradient));
      if (extremeGrad) {
        File origFile = molecularAssembly.getFile();
        String timeString = LocalDateTime.now()
            .format(DateTimeFormatter.ofPattern("yyyy_MM_dd-HH_mm_ss"));

        String filename = format("%s-LARGEGRAD-%s.pdb",
            FilenameUtils.removeExtension(molecularAssembly.getFile().getName()), timeString);
        PotentialsFunctions ef = new PotentialsUtils();
        filename = ef.versionFile(filename);

        logger.warning(
            format(" Excessively large gradients detected; printing snapshot to file %s", filename));
        ef.saveAsPDB(molecularAssembly, new File(filename));
        molecularAssembly.setFile(origFile);
      }
    }

    if (verbose) {
      logger.log(Level.INFO, format("\n OpenMM Energy: %14.10g", e));
    }

    // Scale the coordinates and gradients.
    scaleCoordinatesAndGradient(x, g);

    return e;
  }

  /**
   * Compute the energy and gradient using the pure Java code path.
   *
   * @param x Input atomic coordinates
   * @param g Storage for the gradient vector.
   * @return The energy (kcal/mol)
   */
  public double energyAndGradientFFX(double[] x, double[] g) {
    return super.energyAndGradient(x, g, false);
  }

  /**
   * Compute the energy and gradient using the pure Java code path.
   *
   * @param x       Input atomic coordinates
   * @param g       Storage for the gradient vector.
   * @param verbose Use verbose logging.
   * @return The energy (kcal/mol)
   */
  public double energyAndGradientFFX(double[] x, double[] g, boolean verbose) {
    return super.energyAndGradient(x, g, verbose);
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x Atomic coordinates.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x) {
    return super.energy(x, false);
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x       Input atomic coordinates
   * @param verbose Use verbose logging.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x, boolean verbose) {
    return super.energy(x, verbose);
  }

  /**
   * Returns the Context instance.
   *
   * @return context
   */
  public OpenMMContext getContext() {
    return openMMContext;
  }

  /**
   * Returns the MolecularAssembly instance.
   *
   * @return molecularAssembly
   */
  public MolecularAssembly getMolecularAssembly() {
    return molecularAssembly;
  }

  /**
   * Set the lambdaTerm flag.
   *
   * @param lambdaTerm The value to set.
   */
  public void setLambdaTerm(boolean lambdaTerm) {
    this.lambdaTerm = lambdaTerm;
  }

  /**
   * Get the lambdaTerm flag.
   *
   * @return lambdaTerm.
   */
  public boolean getLambdaTerm() {
    return lambdaTerm;
  }

  /**
   * Re-compute the gradient using OpenMM and return it.
   *
   * @param g Gradient array.
   */
  @Override
  public double[] getGradient(double[] g) {
    OpenMMState openMMState = openMMContext.getOpenMMState(OpenMM_State_Forces);
    g = openMMState.getGradient(g);
    openMMState.destroy();
    return g;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Platform getPlatform() {
    return platform;
  }

  /**
   * Get a reference to the System instance.
   *
   * @return Java wrapper to an OpenMM system.
   */
  public OpenMMSystem getSystem() {
    return openMMSystem;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getd2EdL2() {
    return 0.0;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getdEdL() {
    // No lambda dependence.
    if (!lambdaTerm) {
      return 0.0;
    }

    // Small optimization to only create the x array once.
    double[] x = new double[getNumberOfVariables()];
    getCoordinates(x);

    double currentLambda = getLambda();
    double width = finiteDifferenceStepSize;
    double ePlus;
    double eMinus;

    if (twoSidedFiniteDifference) {
      if (currentLambda + finiteDifferenceStepSize > 1.0) {
        setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = energy(x);
        setLambda(currentLambda);
        ePlus = energy(x);
      } else if (currentLambda - finiteDifferenceStepSize < 0.0) {
        setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = energy(x);
        setLambda(currentLambda);
        eMinus = energy(x);
      } else {
        // Two sided finite difference estimate of dE/dL.
        setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = energy(x);
        setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = energy(x);
        width *= 2.0;
        setLambda(currentLambda);
      }
    } else {
      // One-sided finite difference estimates of dE/dL
      if (currentLambda + finiteDifferenceStepSize > 1.0) {
        setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = energy(x);
        setLambda(currentLambda);
        ePlus = energy(x);
      } else {
        setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = energy(x);
        setLambda(currentLambda);
        eMinus = energy(x);
      }
    }

    // Compute the finite difference derivative.
    double dEdL = (ePlus - eMinus) / width;

    // logger.info(format(" getdEdL currentLambda: CL=%8.6f L=%8.6f dEdL=%12.6f", currentLambda,
    // lambda, dEdL));
    return dEdL;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void getdEdXdL(double[] gradients) {
    // Note for ForceFieldEnergyOpenMM this method is not implemented.
  }

  /**
   * Update active atoms.
   */
  public void setActiveAtoms() {
    openMMSystem.updateAtomMass();
    // Tests show reinitialization of the OpenMM Context is not necessary to pick up mass changes.
    // context.reinitContext();
  }

  /**
   * Set FFX and OpenMM coordinates for active atoms.
   *
   * @param x Atomic coordinates.
   */
  @Override
  public void setCoordinates(double[] x) {
    // Set both OpenMM and FFX coordinates to x.
    openMMContext.setPositions(x);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    super.setCrystal(crystal);
    openMMContext.setPeriodicBoxVectors(crystal);
  }

  public void setLambdaStart(double lambdaStart) {
    this.lambdaStart = lambdaStart;
  }

  public double getLambdaStart() {
    return lambdaStart;
  }

  public void setTwoSidedFiniteDifference(boolean twoSidedFiniteDifference) {
    this.twoSidedFiniteDifference = twoSidedFiniteDifference;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setLambda(double lambda) {
    if (!lambdaTerm) {
      logger.fine(" Attempting to set lambda for an OpenMMEnergy with lambdaterm false.");
      return;
    }

    // Check for lambda outside the range [0 .. 1].
    if (lambda < 0.0 || lambda > 1.0) {
      String message = format(" Lambda value %8.3f is not in the range [0..1].", lambda);
      logger.warning(message);
      return;
    }

    super.setLambda(lambda);

    // Remove the beginning of the normal Lambda path.
    double mappedLambda = lambda;
    if (lambdaStart > 0) {
      double windowSize = 1.0 - lambdaStart;
      mappedLambda = lambdaStart + lambda * windowSize;
    }

    if (openMMSystem != null) {
      openMMSystem.setLambda(mappedLambda);
      if (atoms != null) {
        List<Atom> atomList = new ArrayList<>();
        for (Atom atom : atoms) {
          if (atom.applyLambda()) {
            atomList.add(atom);
          }
        }
        // Update force field parameters based on defined lambda values.
        updateParameters(atomList.toArray(new Atom[0]));
      } else {
        updateParameters(null);
      }
    }
  }

  /**
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  public void updateParameters(@Nullable Atom[] atoms) {
    if (atoms == null) {
      atoms = this.atoms;
    }
    openMMSystem.updateParameters(atoms);
  }

  /**
   * Free OpenMM memory for the Context, Integrator and System.
   */
  private void free() {
    if (openMMContext != null) {
      openMMContext.free();
      openMMContext = null;
    }
    if (openMMSystem != null) {
      openMMSystem.free();
      openMMSystem = null;
    }
  }

}
