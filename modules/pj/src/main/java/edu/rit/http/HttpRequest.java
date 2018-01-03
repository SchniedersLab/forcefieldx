//******************************************************************************
//
// File:    HttpRequest.java
// Package: edu.rit.http
// Unit:    Class edu.rit.http.HttpRequest
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.http;

import java.io.IOException;
import java.net.Socket;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

// For unit test main program
// import java.io.PrintWriter;
// import java.net.InetSocketAddress;
// import java.net.ServerSocket;
// import java.nio.charset.Charset;
/**
 * Class HttpRequest encapsulates an HTTP request received from a web browser.
 * <P>
 * HTTP/1.0 and HTTP/1.1 requests are supported. The obsolete HTTP/0.9 requests
 * are <I>not</I> supported.
 * <P>
 * This class provides methods for examining the request line and the headers.
 * This class does <I>not</I> support reading the entity body if any.
 * <P>
 * To receive an HTTP request message:
 * <OL TYPE=1>
 * <LI>
 * Call the <TT>isValid()</TT> method.
 * <BR>&nbsp;
 * <LI>
 * If <TT>isValid()</TT> returns false, send an HTTP response message indicating
 * the error.
 * <BR>&nbsp;
 * <LI>
 * If <TT>isValid()</TT> returns true, call the other methods to examine the
 * contents of the HTTP request message, and send an appropriate HTTP response
 * message.
 * </OL>
 *
 * @author Alan Kaminsky
 * @version 29-Jul-2010
 */
public class HttpRequest {

// Exported constants.
    /**
     * The GET method string, <TT>"GET"</TT>.
     */
    public static final String GET_METHOD = "GET";

    /**
     * The HEAD method string, <TT>"HEAD"</TT>.
     */
    public static final String HEAD_METHOD = "HEAD";

    /**
     * The POST method string, <TT>"POST"</TT>.
     */
    public static final String POST_METHOD = "POST";

    /**
     * The HTTP/1.0 version string <TT>"HTTP/1.0"</TT>.
     */
    public static final String HTTP_1_0_VERSION = "HTTP/1.0";

    /**
     * The HTTP/1.1 version string, <TT>"HTTP/1.1"</TT>.
     */
    public static final String HTTP_1_1_VERSION = "HTTP/1.1";

// Hidden data members.
    private Socket mySocket;

    private String myMethod;
    private String myUri;
    private String myHttpVersion;

    private Map<String, String> myHeaderMap
            = new HashMap<String, String>();
    private Map<String, String> myUnmodifiableHeaderMap
            = Collections.unmodifiableMap(myHeaderMap);

    private boolean iamValid;

// Exported constructors.
    /**
     * Construct a new HTTP request. The request is read from the input stream
     * of the given socket.
     *
     * @param theSocket Socket.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theSocket</TT> is null.
     */
    public HttpRequest(Socket theSocket) {
        if (theSocket == null) {
            throw new NullPointerException("HttpRequest(): theSocket is null");
        }
        mySocket = theSocket;
    }

// Exported operations.
    /**
     * Determine if this HTTP request is valid. If the data read from the input
     * stream of the socket given to the constructor represents a valid HTTP
     * request message, true is returned, otherwise false is returned. If an I/O
     * exception is thrown while reading the input, this HTTP request is marked
     * as invalid, but the I/O exception is not propagated to the caller.
     *
     * @return True if this HTTP request is valid, false otherwise.
     */
    public boolean isValid() {
        parse();
        return iamValid;
    }

    /**
     * Obtain this HTTP request's method.
     *
     * @return Method string, e.g. <TT>"GET"</TT>, <TT>"POST"</TT>.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * HTTP request is invalid.
     */
    public String getMethod() {
        if (!isValid()) {
            throw new IllegalStateException("HTTP request is invalid");
        }
        return myMethod;
    }

    /**
     * Obtain this HTTP request's URI.
     *
     * @return URI string.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * HTTP request is invalid.
     */
    public String getUri() {
        if (!isValid()) {
            throw new IllegalStateException("HTTP request is invalid");
        }
        return myUri;
    }

    /**
     * Obtain this HTTP request's version.
     *
     * @return HTTP version string, e.g. <TT>"HTTP/1.0"</TT>,
     * <TT>"HTTP/1.1"</TT>.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * HTTP request is invalid.
     */
    public String getHttpVersion() {
        if (!isValid()) {
            throw new IllegalStateException("HTTP request is invalid");
        }
        return myHttpVersion;
    }

    /**
     * Obtain the value of the given header in this HTTP request.
     *
     * @param theHeaderName Header name.
     * @return Header value, or null if there is no header for
     * <TT>theHeaderName</TT>.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * HTTP request is invalid.
     */
    public String getHeader(String theHeaderName) {
        if (!isValid()) {
            throw new IllegalStateException("HTTP request is invalid");
        }
        return myHeaderMap.get(theHeaderName);
    }

    /**
     * Obtain a collection of all the headers in this HTTP request. The returned
     * object is an unmodifiable collection of zero or more map entries. Each
     * map entry's key is the header name. Each map entry's value is the
     * corresponding header value.
     *
     * @return Unmodifiable collection of header name-value mappings.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * HTTP request is invalid.
     */
    public Collection<Map.Entry<String, String>> getHeaders() {
        if (!isValid()) {
            throw new IllegalStateException("HTTP request is invalid");
        }
        return myUnmodifiableHeaderMap.entrySet();
    }

// Hidden operations.
    /**
     * Parse the input data read from this HTTP request's socket.
     */
    private void parse() {
        // Early return if already parsed.
        if (myMethod != null) {
            return;
        }

        // Assume the request is invalid.
        iamValid = false;
        myMethod = "";
        myUri = "";
        myHttpVersion = "";

        try {
            // Set up to scan lines from the socket input stream.
            Scanner scanner = new Scanner(mySocket.getInputStream());

            // Read the first line. If none, invalid.
            if (!scanner.hasNextLine()) {
                return;
            }
            String line = scanner.nextLine();

            // Parse the first line.
            Scanner linescanner = new Scanner(line);
            if (!linescanner.hasNext()) {
                return;
            }
            String method = linescanner.next();
            if (!linescanner.hasNext()) {
                return;
            }
            String uri = linescanner.next();
            if (!linescanner.hasNext()) {
                return;
            }
            String httpVersion = linescanner.next();
            if (linescanner.hasNext()) {
                return;
            }

            // Read remaining lines if any until an empty line.
            String headerName = null;
            String headerValue = "";
            for (;;) {
                if (!scanner.hasNextLine()) {
                    return;
                }
                line = scanner.nextLine();
                if (line.length() == 0) {
                    break;
                }

                // Check whether line is starting or continuing a header.
                if (Character.isWhitespace(line.charAt(0))) {
                    // Continuing previous header.
                    if (headerName == null) {
                        return;
                    }
                    headerValue += line;
                } else {
                    // Starting new header. Record previous header if any.
                    if (headerName != null) {
                        myHeaderMap.put(headerName, headerValue);
                        headerName = null;
                        headerValue = "";
                    }

                    // Parse header name and value.
                    int i = line.indexOf(':');
                    if (i <= 0) {
                        return;
                    }
                    if (i >= line.length() - 1) {
                        return;
                    }
                    if (!Character.isWhitespace(line.charAt(i + 1))) {
                        return;
                    }
                    headerName = line.substring(0, i);
                    headerValue += line.substring(i + 2);
                }
            }

            // If we get here, all is well. Record final header if any.
            if (headerName != null) {
                myHeaderMap.put(headerName, headerValue);
            }

            // Record method, URI, and HTTP version.
            myMethod = method;
            myUri = uri;
            myHttpVersion = httpVersion;

            // Mark it valid.
            iamValid = true;
        } catch (IOException exc) {
            // Leave it marked invalid.
        }
    }

// Unit test main program.
//	/**
//	 * Unit test main program. The program listens for connections to
//	 * localhost:8080. The program reads each HTTP request from a web browser
//	 * and merely echoes the request data back to the browser.
//	 * <P>
//	 * Usage: java edu.rit.http.HttpRequest
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
//			PrintWriter out = new PrintWriter (socket.getOutputStream());
//			out.print ("HTTP/1.0 200 OK\r\n");
//			out.print ("Content-Type: text/plain; charset=");
//			out.print (Charset.defaultCharset() + "\r\n");
//			out.print ("\r\n");
//			if (request.isValid())
//				{
//				out.print ("Method = \"" + request.getMethod() + "\"\r\n");
//				out.print ("URI = \"" + request.getUri() + "\"\r\n");
//				out.print ("Version = \"" + request.getHttpVersion() + "\"\r\n");
//				for (Map.Entry<String,String> entry : request.getHeaders())
//					{
//					out.print ("Header name = \"" + entry.getKey());
//					out.print ("\", value = \"" + entry.getValue() + "\"\r\n");
//					}
//				}
//			else
//				{
//				out.print ("Invalid request\r\n");
//				}
//			out.close();
//			socket.close();
//			}
//		}
}
