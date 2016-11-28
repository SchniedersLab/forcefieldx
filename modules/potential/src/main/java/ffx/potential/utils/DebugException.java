package ffx.potential.utils;

/**
 * 
 * @author slucore
 */
public class DebugException extends RuntimeException {
    private final String message;
    public DebugException() {
        this.message = "Error not specified.";
    }
    public DebugException(String message) {
        this.message = message;
    }
    @Override
    public String getMessage() {
        return "Encountered code under construction. " + message;
    }
}
