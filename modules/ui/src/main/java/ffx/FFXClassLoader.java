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
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

/**
 * Class loader able to load classes and DLLs with a higher priority from a given set of JARs.
 * Its bytecode is Java 1.1 compatible to be loadable by old JVMs.
 *
 * @author Emmanuel Puybaret
 * @version $Id: $
 */
public class FFXClassLoader extends ClassLoader {

    private final ProtectionDomain protectionDomain;
    private final Map extensionDlls = new HashMap();
    private JarFile[] extensionJars = null;
    private final String[] applicationPackages = {"ffx",
                                                  "javax.media.j3d",
                                                  "javax.vecmath",
                                                  "com.sun.j3d",
                                                  "com.sun.opengl",
                                                  "com.sun.gluegen.runtime",
                                                  "javax.media.opengl",
                                                  "groovy",
                                                  "org.codehaus.groovy",
                                                  "org.apache.commons.cli",
                                                  "org.apache.commons.configuration",
                                                  "org.apache.commons.io",
                                                  "org.apache.commons.lang",
                                                  "org.apache.commons.lang3",
                                                  "org.apache.commons.math",
                                                  "edu.rit.pj",
                                                  "jcuda"};
    static final List<String> ffxFiles;

    static {
        ffxFiles = new ArrayList<String>(Arrays.asList(new String[]{
                    "com.kenai.ffx/algorithms.jar",
                    "com.kenai.ffx/autoparm.jar",
                    "com.kenai.ffx/binding.jar",
                    "com.kenai.ffx/crystal.jar",
                    "com.kenai.ffx/numerics.jar",
                    "com.kenai.ffx/potentials.jar",
                    "com.kenai.ffx/ui.jar",
                    "com.kenai.ffx/utilities.jar",
                    "com.kenai.ffx/xray.jar",
                    // Force Field X Extensions
                    "com.kenai.ffx/algorithms-ext.jar",
                    "com.kenai.ffx/xray-ext.jar",
                    // Groovy
                    "org.codehaus.groovy/groovy-all.jar",
                    // CUDA
                    "jcuda/jcuda-all.jar",
                    // Parallel Java
                    "edu.rit.pj/pj.jar",
                    // Java3D 1.5.2 and (so far) 1.6.0 depend on JOGL v. 1.1.1
                    "java3d/j3dcore.jar",
                    "java3d/j3dutils.jar",
                    "java3d/j3dvrml.jar",
                    "java3d/vecmath.jar",
                    // GLUEGEN and JOGL v. 1.1.1
                    "net.java.dev.jogl/jogl.jar",
                    "net.java.dev.jogl/gluegen-rt.jar",
                    // JOGAMP GLUEGEN, JOGL and JOCL v. 2.0
                    // "org.jogamp/gluegen-rt.jar",
                    // "org.jogamp/jogl.jar",
                    // "org.jogamp/nativewindow.jar",
                    // Apache Commons
                    "commons-beanutils/commons-beanutils.jar",
                    "commons-cli/commons-cli.jar",
                    "commons-collections/commons-collections.jar",
                    "commons-configuration/commons-configuration.jar",
                    "commons-digester/commons-digester.jar",
                    "commons-io/commons-io.jar",
                    "commons-lang/commons-lang.jar",
                    "commons-lang/commons-lang3.jar",
                    "commons-logging/commons-logging.jar",
                    "commons-math/commons-math.jar",
                    // Mac OS X Extensions
                    "macosx/AppleJavaExtensions.jar",
                    // Java Help
                    "javax.help/javahelp.jar"
                }));

        String osName = System.getProperty("os.name").toUpperCase();
        String osArch = System.getProperty("sun.arch.data.model");
        boolean x86_64 = "64".equals(osArch);
        if ("MAC OS X".equals(osName)) {
            // JOGL Universal Binaries
            ffxFiles.add("universal/libgluegen-rt.jnilib");
            ffxFiles.add("universal/libjogl.jnilib");
            ffxFiles.add("universal/libjogl_awt.jnilib");
            ffxFiles.add("universal/libjogl_cg.jnilib");
            if (x86_64) {
                // JCUDA
                ffxFiles.add("64-bit/libJCudaDriver-apple-x86_64.jnilib");
                ffxFiles.add("64-bit/libJCudaRuntime-apple-x86_64.jnilib");
                ffxFiles.add("64-bit/libJCufft-apple-x86_64.jnilib");
            }
        } else if ("LINUX".equals(osName)) {
            if (x86_64) {
                // JOGL
                ffxFiles.add("64-bit/libgluegen-rt.so");
                ffxFiles.add("64-bit/libjogl.so");
                ffxFiles.add("64-bit/libjogl_awt.so");
                ffxFiles.add("64-bit/libjogl_cg.so");
                // JCUDA
                ffxFiles.add("64-bit/libJCudaDriver-linux-x86_64.so");
                ffxFiles.add("64-bit/libJCudaRuntime-linux-x86_64.so");
                ffxFiles.add("64-bit/libJCufft-linux-x86_64.so");
            } else {
                ffxFiles.add("32-bit/libgluegen-rt.so");
                ffxFiles.add("32-bit/libjogl.so");
                ffxFiles.add("32-bit/libjogl_awt.so");
                ffxFiles.add("32-bit/libjogl_cg.so");
            }
        } else if (osName.startsWith("WINDOWS")) {
            if (x86_64) {
                // JOGL
                ffxFiles.add("64-bit/gluegen-rt.dll");
                ffxFiles.add("64-bit/jogl.dll");
                ffxFiles.add("64-bit/jogl_cg.dll");
                ffxFiles.add("64-bit/jogl_awt.dll");
                // JCUDA
                ffxFiles.add("64-bit/JCudaDriver-linux-x86_64.dll");
                ffxFiles.add("64-bit/JCudaRuntime-linux-x86_64.dll");
                ffxFiles.add("64-bit/JCufft-linux-x86_64.dll");
            } else {
                ffxFiles.add("32-bit/gluegen-rt.dll");
                ffxFiles.add("32-bit/jogl.dll");
                ffxFiles.add("32-bit/jogl_awt.dll");
                ffxFiles.add("32-bit/jogl_cg.dll");
            }
        }
    }

    /**
     * Force Field X custom class loader considers JARs and DLLs of
     * <code>extensionJarsAndDlls</code> as classpath and libclasspath elements
     * with a higher priority than the ones of default classpath. It will load
     * itself all the classes belonging to packages of <code>applicationPackages</code>.
     *
     * @param parent a {@link java.lang.ClassLoader} object.
     */
    public FFXClassLoader(ClassLoader parent) {
        super(parent);
        protectionDomain = FFXClassLoader.class.getProtectionDomain();
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     *
     * @param input a {@link java.io.InputStream} object.
     * @param suffix a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public static String copyInputStreamToTmpFile(InputStream input, String name,
                                                  String suffix) throws IOException {
        File tmpFile = null;
        name = "ffx." + name + ".";
        try {
            tmpFile = File.createTempFile(name, suffix);
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
     * {@inheritDoc}
     *
     * Finds and defines the given class among the extension JARs
     * given in constructor, then among resources.
     */
    @Override
    protected Class findClass(String name) throws ClassNotFoundException {
        if (!extensionsLoaded) {
            loadExtensions();
        }

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
     * {@inheritDoc}
     *
     * Returns the library path of an extension DLL.
     */
    @Override
    protected String findLibrary(String libname) {
        if (!extensionsLoaded) {
            loadExtensions();
        }

        return (String) this.extensionDlls.get(libname);
    }

    /**
     * {@inheritDoc}
     *
     * Returns the URL of the given resource searching first if it exists among
     * the extension JARs given in constructor.
     */
    @Override
    protected URL findResource(String name) {
        if (!extensionsLoaded) {
            loadExtensions();
        }

        if (name.equals("List all scripts")) {
            listScripts();
        }

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

    protected void listScripts() {
        if (this.extensionJars != null) {

            List<String> scripts = new ArrayList<String>();

            // Check if searched class is an extension class
            for (int i = 0; i < this.extensionJars.length; i++) {
                JarFile extensionJar = this.extensionJars[i];
                Enumeration<JarEntry> entries = extensionJar.entries();
                while (entries.hasMoreElements()) {
                    JarEntry entry = entries.nextElement();
                    String name = entry.getName();
                    if (name.startsWith("ffx") && name.endsWith(".groovy")) {
                        name = name.replace('/', '.');
                        name = name.replace("ffx.scripts.", "");
                        name = name.replace(".groovy", "");
                        scripts.add(name);
                    }
                }
            }

            String[] scriptArray = scripts.toArray(new String[scripts.size()]);
            Arrays.sort(scriptArray);
            for (String script : scriptArray) {
                System.out.println(" " + script);
            }
        }
    }

    /**
     * {@inheritDoc}
     *
     * Loads a class with this class loader if its package belongs to <code>applicationPackages</code>
     * given in constructor.
     */
    @Override
    protected Class loadClass(String name, boolean resolve) throws ClassNotFoundException {
        if (!extensionsLoaded) {
            loadExtensions();
        }

        // If no extension jars were found use the super.loadClass method.
        if (this.extensionJars == null) {
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
                // Do Nothin.
            }
            // Try to load the class via the default implementation.
            if (loadedClass == null) {
                loadedClass = super.loadClass(name, resolve);
            }
        }
        if (resolve) {
            resolveClass(loadedClass);
        }
        return loadedClass;
    }
    private boolean extensionsLoaded = false;

    private void loadExtensions() {
        if (extensionsLoaded) {
            return;
        }
        extensionsLoaded = true;

        String extensionJarsAndDlls[] = ffxFiles.toArray(new String[ffxFiles.size()]);

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
                    int lastSlashIndex = extensionJarOrDll.lastIndexOf('/');
                    if (extensionJarOrDll.endsWith(".jar")) {
                        int start = lastSlashIndex + 1;
                        int end = extensionJarOrDll.indexOf(".jar");
                        String name = extensionJarOrDll.substring(start, end);
                        // Copy jar to a tmp file
                        String extensionJar = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(),
                                                                       name, ".jar");
                        // Add extracted file to the extension jars list
                        extensionJarList.add(new JarFile(extensionJar, false));
                    } else if (extensionJarOrDll.endsWith(dllSuffix)) {
                        int start = lastSlashIndex + 1 + dllPrefix.length();
                        int end = extensionJarOrDll.indexOf(dllSuffix);
                        String name = extensionJarOrDll.substring(start, end);
                        // Copy DLL to a tmp file
                        String extensionDll = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(),
                                                                       name, dllSuffix);
                        // Add extracted file to extension DLLs map
                        extensionDlls.put(name, extensionDll);
                    }
                }
            } catch (IOException ex) {
                throw new RuntimeException(" Couldn't extract extension jar.\n", ex);
            }
        }

        // Create extensionJars array
        if (extensionJarList.size() > 0) {
            this.extensionJars = (JarFile[]) extensionJarList.toArray(new JarFile[extensionJarList.size()]);
        }
    }
}
