package ffx.xray.refine;

import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.toEnumForm;
import static java.lang.String.format;

/**
 * Different refinement mode selection types.
 */
public enum RefinementMode {

  /**
   * refine B factors only (if anisotropic, refined as such)
   */
  BFACTORS,
  /**
   * refine B factors and occupancies
   */
  BFACTORS_AND_OCCUPANCIES,
  /**
   * refine coordinates only
   */
  COORDINATES,
  /**
   * refine coordinates and B factors (if anisotropic, refined as such)
   */
  COORDINATES_AND_BFACTORS,
  /**
   * refine all
   */
  COORDINATES_AND_BFACTORS_AND_OCCUPANCIES,
  /**
   * refine coordinates and occupancies
   */
  COORDINATES_AND_OCCUPANCIES,
  /**
   * refine occupancies only (alternate conformers are constrained)
   */
  OCCUPANCIES;

  /**
   * Determines whether the refinement mode includes coordinates.
   *
   * @return true if the refinement mode involves coordinates, false otherwise.
   */
  public boolean includesCoordinates() {
    return this == COORDINATES ||
        this == COORDINATES_AND_BFACTORS ||
        this == COORDINATES_AND_OCCUPANCIES ||
        this == COORDINATES_AND_BFACTORS_AND_OCCUPANCIES;
  }

  /**
   * Determines whether the refinement mode includes B-factors.
   *
   * @return true if the refinement mode involves B-factors, false otherwise.
   */
  public boolean includesBFactors() {
    return this == BFACTORS ||
        this == BFACTORS_AND_OCCUPANCIES ||
        this == COORDINATES_AND_BFACTORS ||
        this == COORDINATES_AND_BFACTORS_AND_OCCUPANCIES;
  }

  /**
   * Determines whether the refinement mode includes occupancies.
   *
   * @return true if the refinement mode involves occupancies, false otherwise.
   */
  public boolean includesOccupancies() {
    return this == OCCUPANCIES ||
        this == BFACTORS_AND_OCCUPANCIES ||
        this == COORDINATES_AND_OCCUPANCIES ||
        this == COORDINATES_AND_BFACTORS_AND_OCCUPANCIES;
  }

  public String toString() {
    return switch (this) {
      case BFACTORS -> " Mode: B-Factors";
      case BFACTORS_AND_OCCUPANCIES -> " Mode: B-Factors and Occupancies";
      case COORDINATES -> " Mode: Coordinates";
      case COORDINATES_AND_BFACTORS -> " Mode: Coordinates and B-Factors";
      case COORDINATES_AND_BFACTORS_AND_OCCUPANCIES -> " Mode: Coordinates, B-Factors and Occupancies";
      case COORDINATES_AND_OCCUPANCIES -> " Mode: Coordinates and Occupancies";
      case OCCUPANCIES -> " Mode: Occupancies";
    };
  }

  /**
   * Returns the default epsilon value for the current refinement mode.
   *
   * @return The default epsilon value for the current refinement mode.
   */
  public double getDefaultEps() {
    return getDefaultEps(false);
  }

  /**
   * Returns the default epsilon value for the current refinement mode.
   * The returned value depends on whether the refinement mode includes anisotropic
   * displacement parameters (ANISOU) and the specific refinement mode.
   *
   * @param hasAnisou Indicates whether the refinement mode includes anisotropic displacement parameters.
   * @return The default epsilon value for the current refinement mode.
   */
  public double getDefaultEps(boolean hasAnisou) {
    return switch (this) {
      case BFACTORS, BFACTORS_AND_OCCUPANCIES -> {
        if (hasAnisou) {
          yield 20.0;
        } else {
          yield 0.01;
        }
      }
      case COORDINATES -> 0.4;
      case COORDINATES_AND_BFACTORS, COORDINATES_AND_BFACTORS_AND_OCCUPANCIES -> {
        if (hasAnisou) {
          yield 20.0;
        } else {
          yield 0.2;
        }
      }
      case COORDINATES_AND_OCCUPANCIES -> 0.2;
      case OCCUPANCIES -> 0.1;
    };
  }

  /**
   * Parse a string into a refinement mode.
   *
   * @param mode Refinement mode string.
   * @return An instance of RefinementMode.
   */
  public static RefinementMode parseMode(String mode) {
    try {
      return RefinementMode.valueOf(toEnumForm(mode));
    } catch (Exception e) {
      Logger logger = Logger.getLogger(RefinementMode.class.getName());
      logger.info(format(" Could not parse %s as a refinement mode; defaulting to coordinates.", mode));
      return RefinementMode.COORDINATES;
    }
  }

}
