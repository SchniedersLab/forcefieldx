//******************************************************************************
//
// File:    HttpResponse.java
// Package: edu.rit.http
// Unit:    Class edu.rit.http.HttpResponse
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.http;

import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.Socket;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;

// For unit test main program
// import java.net.InetSocketAddress;
// import java.net.ServerSocket;
/**
 * Class HttpResponse encapsulates an HTTP response returned to a web browser.
 * <P>
 * Only HTTP/1.0 responses are supported. This means that only one HTTP response
 * message can be sent over the connection to the web browser; the connection is
 * closed after sending the HTTP response message.
 * <P>
 * To send an HTTP response message:
 * <OL TYPE=1>
 * <LI>
 * Create an instance of class HttpResponse, giving the socket object
 * representing the open connection to the web browser.
 * <BR>&nbsp;
 * <LI>
 * Call methods to set the status code and headers as necessary.
 * <BR>&nbsp;
 * <LI>
 * Call the <code>getPrintWriter()</code> method, and write the entity body to the
 * print writer that is returned.
 * <BR>&nbsp;
 * <LI>
 * Call the <code>close()</code> method.
 * </OL>
 *
 * @author Alan Kaminsky
 * @version 29-Jul-2010
 */
public class HttpResponse {

// Exported enumerations.
    /**
     * Enumeration HttpResponse.Status enumerates the status codes for an HTTP
     * response message.
     *
     * @author Alan Kaminsky
     * @version 10-Oct-2006
     */
    public static enum Status {

        /**
         * The request has succeeded.
         */
        STATUS_200_OK("200 OK"),
        /**
         * The request has been fulfilled and resulted in a new resource being
         * created.
         */
        STATUS_201_CREATED("201 Created"),
        /**
         * The request has been accepted for processing, but the processing has
         * not been completed.
         */
        STATUS_202_ACCEPTED("202 Accepted"),
        /**
         * The server has fulfilled the request but there is no new information
         * to send back.
         */
        STATUS_204_NO_CONTENT("204 No content"),
        /**
         * The requested resource has been assigned a new permanent URL and any
         * future references to this resource should be done using that URL.
         */
        STATUS_301_MOVED_PERMANENTLY("301 Moved Permanently"),
        /**
         * The requested resource resides temporarily under a different URL.
         */
        STATUS_302_MOVED_TEMPORARILY("302 Moved Temporarily"),
        /**
         * If the client has performed a conditional GET request and access is
         * allowed, but the document has not been modified since the date and
         * time specified in the If-Modified-Since field, the server must
         * respond with this status code and not send an Entity-Body to the
         * client.
         */
        STATUS_304_NOT_MODIFIED("304 Not Modified"),
        /**
         * The request could not be understood by the server due to malformed
         * syntax.
         */
        STATUS_400_BAD_REQUEST("400 Bad Request"),
        /**
         * The request requires user authentication.
         */
        STATUS_401_UNAUTHORIZED("401 Unauthorized"),
        /**
         * The server understood the request, but is refusing to fulfill it.
         */
        STATUS_403_FORBIDDEN("403 Forbidden"),
        /**
         * The server has not found anything matching the Request-URI.
         */
        STATUS_404_NOT_FOUND("404 Not Found"),
        /**
         * The server encountered an unexpected condition which prevented it
         * from fulfilling the request.
         */
        STATUS_500_INTERNAL_SERVER_ERROR("500 Internal Server Error"),
        /**
         * The server does not support the functionality required to fulfill the
         * request.
         */
        STATUS_501_NOT_IMPLEMENTED("501 Not Implemented"),
        /**
         * The server, while acting as a gateway or proxy, received an invalid
         * response from the upstream server it accessed in attempting to
         * fulfill the request.
         */
        STATUS_502_BAD_GATEWAY("502 Bad Gateway"),
        /**
         * The server is currently unable to handle the request due to a
         * temporary overloading or maintenance of the server.
         */
        STATUS_503_SERVICE_UNAVAILABLE("503 Service Unavailable");

        private final String stringForm;

        /**
         * Construct a new Status value.
         *
         * @param stringForm String form.
         */
        Status(String stringForm) {
            this.stringForm = stringForm;
        }

        /**
         * Returns a string version of this Status value.
         *
         * @return String version.
         */
        public String toString() {
            return stringForm;
        }
    }

// Hidden data members.
    private Socket mySocket;

    private Status myStatusCode = Status.STATUS_200_OK;
    private String myContentType = "text/html";
    private Charset myCharset = Charset.defaultCharset();

    private Map<String, String> myHeaderMap = new HashMap<String, String>();

    private PrintWriter myPrintWriter;

// Exported constructors.
    /**
     * Construct a new HTTP response. The response is written to the output
     * stream of the given socket.
     *
     * @param theSocket Socket.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theSocket</code> is null.
     */
    public HttpResponse(Socket theSocket) {
        if (theSocket == null) {
            throw new NullPointerException("HttpResponse(): theSocket is null");
        }
        mySocket = theSocket;
        recordContentType();
    }

// Exported operations.
    /**
     * Set this HTTP response's status code. If not set, the default status code
     * is STATUS_200_OK.
     *
     * @param theStatusCode Status code.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theStatusCode</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if the HTTP
     * response headers have already been written to the socket output stream.
     */
    public void setStatusCode(Status theStatusCode) {
        if (theStatusCode == null) {
            throw new NullPointerException("HttpResponse.setStatusCode(): theStatusCode is null");
        }
        if (myPrintWriter != null) {
            throw new IllegalStateException("HttpResponse.setStatusCode(): Headers already written");
        }
        myStatusCode = theStatusCode;
    }

    /**
     * Set this HTTP response's content type. If not set, the default content
     * type is <code>"text/html"</code>.
     *
     * @param theContentType Content type.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theContentType</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>theContentType</code> is zero length.
     * @exception IllegalStateException (unchecked exception) Thrown if the HTTP
     * response headers have already been written to the socket output stream.
     */
    public void setContentType(String theContentType) {
        if (theContentType == null) {
            throw new NullPointerException("HttpResponse.setContentType(): theContentType is null");
        }
        if (theContentType.length() == 0) {
            throw new IllegalArgumentException("HttpResponse.setContentType(): theContentType is zero length");
        }
        if (myPrintWriter != null) {
            throw new IllegalStateException("HttpResponse.setContentType(): Headers already written");
        }
        myContentType = theContentType;
        recordContentType();
    }

    /**
     * Set this HTTP response's character set. If not set, the default character
     * set is the platform default character set.
     *
     * @param theCharset Character set.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theCharset</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if the HTTP
     * response headers have already been written to the socket output stream.
     */
    public void setCharset(Charset theCharset) {
        if (theCharset == null) {
            throw new NullPointerException("HttpResponse.setCharset(): theCharset is null");
        }
        if (myPrintWriter != null) {
            throw new IllegalStateException("HttpResponse.setCharset(): Headers already written");
        }
        myCharset = theCharset;
        recordContentType();
    }

    /**
     * Set this HTTP response's content length. If not set, the default is to
     * omit the Content-Length header; closing this HTTP response marks the end
     * of the entity body.
     *
     * @param theContentLength Content length.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>theContentLength</code> is less than 0.
     * @exception IllegalStateException (unchecked exception) Thrown if the HTTP
     * response headers have already been written to the socket output stream.
     */
    public void setContentLength(int theContentLength) {
        if (theContentLength < 0) {
            throw new IllegalArgumentException("HttpResponse.setContentLength(): theContentLength < 0");
        }
        if (myPrintWriter != null) {
            throw new IllegalStateException("HttpResponse.setContentLength(): Headers already written");
        }
        myHeaderMap.put("Content-Length", "" + theContentLength);
    }

    /**
     * Set the given header in this HTTP response.
     *
     * @param theHeaderName Header name.
     * @param theHeaderValue Header value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theHeaderName</code> or
     * <code>theHeaderValue</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>theHeaderName</code> or
     * <code>theHeaderValue</code> is zero length.
     * @exception IllegalStateException (unchecked exception) Thrown if the HTTP
     * response headers have already been written to the socket output stream.
     */
    public void setHeader(String theHeaderName,
            String theHeaderValue) {
        if (theHeaderName == null) {
            throw new NullPointerException("HttpResponse.setHeader(): theHeaderName is null");
        }
        if (theHeaderName.length() == 0) {
            throw new IllegalArgumentException("HttpResponse.setHeader(): theHeaderName is zero length");
        }
        if (theHeaderValue == null) {
            throw new NullPointerException("HttpResponse.setHeader(): theHeaderValue is null");
        }
        if (theHeaderValue.length() == 0) {
            throw new IllegalArgumentException("HttpResponse.setHeader(): theHeaderValue is zero length");
        }
        if (myPrintWriter != null) {
            throw new IllegalStateException("HttpResponse.setHeader(): Headers already written");
        }
        myHeaderMap.put(theHeaderName, theHeaderValue);
    }

    /**
     * Obtain the print writer for writing the entity body to this HTTP
     * response. As a side effect, the HTTP response headers are written to the
     * socket output stream.
     *
     * @return Print writer.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public PrintWriter getPrintWriter()
            throws IOException {
        return writeHeaders();
    }

    /**
     * Close this HTTP response. If necessary, the HTTP response headers are
     * written to the socket output stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void close()
            throws IOException {
        writeHeaders();
        myPrintWriter.close();
        mySocket.close();
    }

// Hidden operations.
    /**
     * Record the Content-Type header.
     */
    private void recordContentType() {
        myHeaderMap.put("Content-Type",
                myContentType + "; charset=" + myCharset);
    }

    /**
     * Write the headers to the socket output stream if not already written, and
     * return the print writer for writing to the socket output stream.
     *
     * @return Print writer.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private PrintWriter writeHeaders()
            throws IOException {
        if (myPrintWriter == null) {
            myPrintWriter
                    = new PrintWriter(new OutputStreamWriter(mySocket.getOutputStream(),
                                    myCharset));
            myPrintWriter.write("HTTP/1.0 " + myStatusCode + "\r\n");
            for (Map.Entry<String, String> entry : myHeaderMap.entrySet()) {
                myPrintWriter.write(entry.getKey() + ": " + entry.getValue() + "\r\n");
            }
            myPrintWriter.write("\r\n");
        }
        return myPrintWriter;
    }

// Unit test main program.
//	/**
//	 * Unit test main program. The program listens for connections to
//	 * localhost:8080. The program reads each HTTP request from a web browser
//	 * and merely echoes the request data back to the browser.
//	 * <P>
//	 * Usage: java edu.rit.http.HttpResponse
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		ServerSocket serversocket = new ServerSocket();
//		serversocket.bind (new InetSocketAddress ("localhost", 8080));
//		for (;;)
//			{
//			Socket socket = serversocket.accept();
//			HttpRequest request = new HttpRequest (socket);
//			HttpResponse response = new HttpResponse (socket);
//			if (request.isValid())
//				{
//				PrintWriter out = response.getPrintWriter();
//				out.println ("<HTML>");
//				out.println ("<HEAD>");
//				out.println ("</HEAD>");
//				out.println ("<BODY>");
//				out.println ("<UL>");
//				out.println ("<LI>");
//				out.print   ("Method = <code>\"");
//				out.print   (request.getMethod());
//				out.println ("\"</code>");
//				out.println ("<LI>");
//				out.print   ("URI = <code>\"");
//				out.print   (request.getUri());
//				out.println ("\"</code>");
//				out.println ("<LI>");
//				out.print   ("Version = <code>\"");
//				out.print   (request.getHttpVersion());
//				out.println ("\"</code>");
//				for (Map.Entry<String,String> entry : request.getHeaders())
//					{
//					out.println ("<LI>");
//					out.print   ("Header name = <code>\"");
//					out.print   (entry.getKey());
//					out.print   ("\"</code>, value = <code>\"");
//					out.print   (entry.getValue());
//					out.println ("\"</code>");
//					}
//				out.println ("</UL>");
//				out.println ("</BODY>");
//				out.println ("</HTML>");
//				}
//			else
//				{
//				response.setStatusCode (Status.STATUS_400_BAD_REQUEST);
//				PrintWriter out = response.getPrintWriter();
//				out.println ("<HTML>");
//				out.println ("<HEAD>");
//				out.println ("</HEAD>");
//				out.println ("<BODY>");
//				out.println ("<P>");
//				out.println ("400 Bad Request");
//				out.println ("<P>");
//				out.println ("You idiot.");
//				out.println ("</BODY>");
//				out.println ("</HTML>");
//				}
//			response.close();
//			}
//		}
}
