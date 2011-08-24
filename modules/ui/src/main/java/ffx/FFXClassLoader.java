/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.security.ProtectionDomain;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Class loader able to load classes and DLLs with a higher priority from a given set of JARs.
 * Its bytecode is Java 1.1 compatible to be loadable by old JVMs.
 * @author Emmanuel Puybaret
 */
public class FFXClassLoader extends ClassLoader {

    private final ProtectionDomain protectionDomain;
    private final String[] applicationPackages;
    private final Map extensionDlls = new HashMap();
    private JarFile[] extensionJars = null;

    /**
     * Creates a class loader. It will consider JARs and DLLs of <code>extensionJarsAndDlls</code>
     * as classpath and libclasspath elements with a higher priority than the ones of default classpath,
     * and will load itself all the classes belonging to packages of <code>applicationPackages</code>.
     */
    public FFXClassLoader(ClassLoader parent,
                          ProtectionDomain protectionDomain,
                          String[] extensionJarsAndDlls,
                          String[] applicationPackages) {
        super(parent);
        this.protectionDomain = protectionDomain;
        this.applicationPackages = applicationPackages;
        
        // Compute DLLs prefix and suffix
        String dllSuffix;
        String dllPrefix;

        String osName = System.getProperty("os.name");
        if (osName.startsWith("Windows")) {
            dllSuffix = ".dll";
            dllPrefix = "";
        } else if (osName.startsWith("Mac OS X")) {
            dllSuffix = ".jnilib";
            dllPrefix = "lib";
        } else {
            dllSuffix = ".so";
            dllPrefix = "lib";
        }

        // Find extension Jars and DLLs
        ArrayList extensionJarList = new ArrayList();
        for (int i = 0; i < extensionJarsAndDlls.length; i++) {
            String extensionJarOrDll = extensionJarsAndDlls[i];
            try {
                URL extensionJarOrDllUrl = getResource(extensionJarOrDll);
                if (extensionJarOrDllUrl != null) {
                    if (extensionJarOrDll.endsWith(".jar")) {
                        // Copy jar to a tmp file
                        String extensionJar = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(), ".jar");
                        // Add tmp file to extension jars list
                        extensionJarList.add(new JarFile(extensionJar, false));
                    } else if (extensionJarOrDll.endsWith(dllSuffix)) {
                        int lastSlashIndex = extensionJarOrDll.lastIndexOf('/');
                        // Copy DLL to a tmp file
                        String extensionDll = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(), dllSuffix);
                        // Add tmp file to extension DLLs map
                        this.extensionDlls.put(extensionJarOrDll.substring(lastSlashIndex + 1 + dllPrefix.length(),
                                                                           extensionJarOrDll.indexOf(dllSuffix)), extensionDll);
                    }
                } else {
                    System.out.println(" File not extracted: " + extensionJarOrDll);
                }
            } catch (IOException ex) {
                throw new RuntimeException(" Couldn't extract extension jar.\n", ex);
            }
        }

        // Add extension jars from the ffx/conf directory
        try {
        String confDir = System.getProperty("basedir") + File.separator + "conf";
        File confFile = new File(confDir);
        File files[] = confFile.listFiles();        
        for (File jar : files) {
            try {
                String ext = jar.getCanonicalPath().toUpperCase();
                if (ext.endsWith("JAR")) {
                    System.out.println(" Loading extension jar... " + jar);
                    extensionJarList.add(new JarFile(jar, false));
                }
            } catch (Exception e) {
                throw new RuntimeException(" Couldn't load extension jar.\n", e);
            }
        } } catch (Exception e) {
            
        }

        // Create extensionJars array
        if (extensionJarList.size() > 0) {
            this.extensionJars = (JarFile[]) extensionJarList.toArray(new JarFile[extensionJarList.size()]);
        }
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     */
    public static String copyInputStreamToTmpFile(InputStream input,
                                                  String suffix) throws IOException {
        File tmpFile = null;
        try {
            tmpFile = File.createTempFile("tmp.", suffix);
        } catch (Exception e) {
            System.out.println(" Could not extract a Force Field X library.");
            System.err.println(e.toString());
            System.exit(-1);
        }
        tmpFile.deleteOnExit();
        OutputStream output = null;
        try {
            output = new BufferedOutputStream(new FileOutputStream(tmpFile));
            byte[] buffer = new byte[8192];
            int size;
            while ((size = input.read(buffer)) != -1) {
                output.write(buffer, 0, size);
            }
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
            }
        }
        return tmpFile.toString();
    }

    /**
     * Finds and defines the given class among the extension JARs
     * given in constructor, then among resources.
     */
    @Override
    protected Class findClass(String name) throws ClassNotFoundException {
        // Build class file from its name
        String classFile = name.replace('.', '/') + ".class";
        InputStream classInputStream = null;
        if (this.extensionJars != null) {
            // Check if searched class is an extension class
            for (int i = 0; i < this.extensionJars.length; i++) {
                JarFile extensionJar = this.extensionJars[i];
                JarEntry jarEntry = extensionJar.getJarEntry(classFile);
                if (jarEntry != null) {
                    try {
                        classInputStream = extensionJar.getInputStream(jarEntry);
                    } catch (IOException ex) {
                        throw new ClassNotFoundException("Couldn't read class " + name, ex);
                    }
                }
            }
        }
        // If it's not an extension class, search if its an application
        // class that can be read from resources
        if (classInputStream == null) {
            URL url = getResource(classFile);
            if (url == null) {
                throw new ClassNotFoundException("Class " + name);
            }
            try {
                classInputStream = url.openStream();
            } catch (IOException ex) {
                throw new ClassNotFoundException("Couldn't read class " + name, ex);
            }
        }

        try {
            // Read class input content to a byte array
            ByteArrayOutputStream out = new ByteArrayOutputStream();
            BufferedInputStream in = new BufferedInputStream(classInputStream);
            byte[] buffer = new byte[8192];
            int size;
            while ((size = in.read(buffer)) != -1) {
                out.write(buffer, 0, size);
            }
            in.close();
            // Define class
            return defineClass(name, out.toByteArray(), 0, out.size(),
                               this.protectionDomain);
        } catch (IOException ex) {
            throw new ClassNotFoundException("Class " + name, ex);
        }
    }

    /**
     * Returns the library path of an extension DLL.
     */
    @Override
    protected String findLibrary(String libname) {
        return (String) this.extensionDlls.get(libname);
    }

    /**
     * Returns the URL of the given resource searching first if it exists among
     * the extension JARs given in constructor.
     */
    @Override
    protected URL findResource(String name) {
        if (extensionJars != null) {
            // Try to find if resource belongs to one of the extracted jars
            for (int i = 0; i < extensionJars.length; i++) {
                JarFile extensionJar = extensionJars[i];
                JarEntry jarEntry = extensionJar.getJarEntry(name);
                if (jarEntry != null) {
                    String path = "jar:file:" + extensionJar.getName() + "!/" + jarEntry.getName();
                    try {
                        return new URL(path);
                    } catch (MalformedURLException ex) {
                        System.out.println(path + "\n" + ex.toString());
                    }
                }
            }
        }
        return super.findResource(name);
    }

    /**
     * Loads a class with this class loader if its package belongs to <code>applicationPackages</code>
     * given in constructor.
     */
    @Override
    protected Class loadClass(String name, boolean resolve) throws ClassNotFoundException {
        // If no extension jars couldn't be found
        if (this.extensionJars == null) {
            // Let default class loader do its job
            return super.loadClass(name, resolve);
        }
        // Check if the class has already been loaded
        Class loadedClass = findLoadedClass(name);
        if (loadedClass == null) {
            try {
                // Try to find if class belongs to one of the application packages
                for (int i = 0; i < this.applicationPackages.length; i++) {
                    String applicationPackage = this.applicationPackages[i];
                    int applicationPackageLength = applicationPackage.length();
                    if ((applicationPackageLength == 0
                         && name.indexOf('.') == 0)
                        || (applicationPackageLength > 0
                            && name.startsWith(applicationPackage))) {
                        loadedClass = findClass(name);
                        break;
                    }
                }
            } catch (ClassNotFoundException ex) {
                // Let a chance to class to be loaded by default implementation
            }
            if (loadedClass == null) {
                loadedClass = super.loadClass(name, resolve);
            }
        }
        if (resolve) {
            resolveClass(loadedClass);
        }
        return loadedClass;
    }
}
