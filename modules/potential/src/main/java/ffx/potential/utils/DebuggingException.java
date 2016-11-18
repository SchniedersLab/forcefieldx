package ffx.potential.utils;

/**
 * 
 * @author slucore
 */
public class DebuggingException extends RuntimeException {
    private final String location;
    public DebuggingException() {
        this.location = "unspecified";
    }
    public DebuggingException(String location) {
        this.location = location;
    }
    @Override
    public String getMessage() {
        return String.format("Encountered code under construction @ location: %s.", location);
    }
}
