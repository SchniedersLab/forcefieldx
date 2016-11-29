package ffx.potential.utils;

/**
 * 
 * @author slucore
 */
public class SystemTemperatureException extends RuntimeException {
    final double tooHot;
    public SystemTemperatureException(double temperature) {
        tooHot = temperature;
    }
    @Override
    public String getMessage() {
        return String.format("System unstable due to high temperature: %.2f", tooHot);
    }
}
