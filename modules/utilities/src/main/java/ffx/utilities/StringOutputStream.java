package ffx.utilities;

import java.io.PrintStream;
import java.io.ByteArrayOutputStream;
import java.nio.charset.Charset;
import java.io.UnsupportedEncodingException;

public class StringOutputStream extends PrintStream {

    ByteArrayOutputStream baos = null;
    Charset charset = null;

    public StringOutputStream(ByteArrayOutputStream baos) throws UnsupportedEncodingException {
        super(baos, true, Charset.defaultCharset().displayName());
        this.baos = baos;
        charset = Charset.defaultCharset();
    }

    public StringOutputStream(ByteArrayOutputStream baos, Charset charset) throws UnsupportedEncodingException  {
        super(baos, true, charset.displayName());
        this.baos = baos;
        this.charset = charset;
    }

    public String toString() {
        return new String(baos.toByteArray(), charset);
    }

    public void close() {
        super.close();
        try {
            baos.close();
        } catch (Exception e) {
            //
        }
    }

}
