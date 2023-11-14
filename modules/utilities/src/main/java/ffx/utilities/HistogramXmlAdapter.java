package ffx.utilities;

import javax.xml.bind.annotation.adapters.XmlAdapter;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.logging.Logger;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;

/**
 * Write/read Histogram data to/from a String.
 * The String is an array of lines separated by System dependent line separators.
 * Each line (row) corresponds to a fixed lambda.
 * The columns are a space-separated list of weights for the dUdL bins.
 *
 * @author Michael J. Schnieders
 * @see javax.xml.bind.annotation.adapters.XmlAdapter
 * @since 1.0
 */
public class HistogramXmlAdapter extends XmlAdapter<String, double[][]> {

  /**
   * The logger for this class.
   */
  public static final Logger logger = Logger.getLogger(HistogramXmlAdapter.class.getName());

  /**
   * Convert the histogram data string to a 2D double array.
   *
   * @param value String containing the Histogram data.
   * @return The data.
   */
  @Override
  public double[][] unmarshal(String value) {
    double[][] data;
    try {
      String[] lines = value.split(System.lineSeparator());
      int lambdaBins = lines.length;
      data = new double[lambdaBins][];
      for (int i = 0; i < lambdaBins; i++) {
        String[] countToks = lines[i].split(" +");
        int dUdLBins = countToks.length;
        data[i] = new double[dUdLBins];
        double[] row = data[i];
        for (int j = 0; j < dUdLBins; j++) {
          row[j] = parseDouble(countToks[j]);
        }
      }
    } catch (Exception e) {
      logger.warning(" Returning a null histogram due to an XML parsing exception:\n " + e);
      data = null;
    }
    return data;
  }

  /**
   * Convert the 2D histogram double array into a String.
   *
   * @param data The histogram data.
   * @return A String representation of the data.
   */
  @Override
  public String marshal(double[][] data) {
    StringBuilder sb = new StringBuilder();
    for (double[] row : data) {
      sb.append(format("%.16e", row[0]));
      int dUdLBins = row.length;
      for (int j = 1; j < dUdLBins; j++) {
        sb.append(" ").append(format("%.16e", row[j]));
      }
      sb.append(System.lineSeparator());
    }
    return sb.toString();
  }
}

