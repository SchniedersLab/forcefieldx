package ffx.algorithms.dynamics;

/**
 * Define the verbosity of the MolecularDynamics class.
 */
public enum MDVerbosity {
  VERBOSE(false), QUIET(true), SILENT(true);

  private final boolean isQuiet;

  MDVerbosity(boolean isQuiet) {
    this.isQuiet = isQuiet;
  }

  public boolean isQuiet() {
    return isQuiet;
  }
}
