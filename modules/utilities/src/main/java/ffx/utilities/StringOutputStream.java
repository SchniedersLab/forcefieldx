package ffx.utilities;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.nio.charset.Charset;

/**
 * <p>StringOutputStream class.</p>
 *
 * @author Michael J. Schnieders
 */
public class StringOutputStream extends PrintStream {

    ByteArrayOutputStream baos = null;
    Charset charset = null;

    /**
     * <p>Constructor for StringOutputStream.</p>
     *
     * @param baos a {@link java.io.ByteArrayOutputStream} object.
     * @throws java.io.UnsupportedEncodingException if any.
     */
    public StringOutputStream(ByteArrayOutputStream baos) throws UnsupportedEncodingException {
        super(baos, true, Charset.defaultCharset().displayName());
        this.baos = baos;
        charset = Charset.defaultCharset();
    }

    /**
     * <p>Constructor for StringOutputStream.</p>
     *
     * @param baos    a {@link java.io.ByteArrayOutputStream} object.
     * @param charset a {@link java.nio.charset.Charset} object.
     * @throws java.io.UnsupportedEncodingException if any.
     */
    public StringOutputStream(ByteArrayOutputStream baos, Charset charset) throws UnsupportedEncodingException {
        super(baos, true, charset.displayName());
        this.baos = baos;
        this.charset = charset;
    }

    /**
     * <p>toString.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toString() {
        return new String(baos.toByteArray(), charset);
    }

    /**
     * <p>close.</p>
     */
    public void close() {
        super.close();
        try {
            baos.close();
        } catch (Exception e) {
            //
        }
    }

}
