/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
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
 * Class loader able to load classes and DLLs with a higher priority from a
 * given set of JARs. Its bytecode is Java 1.1 compatible to be loadable by old
 * JVMs.
 *
 * @author Michael J. Schnieders; derived from work by Emmanuel Puybaret
 *
 */
public class FFXClassLoader extends URLClassLoader {

    private final ProtectionDomain protectionDomain;
    private final Map<String, String> extensionDlls = new HashMap<>();
    private JarFile[] extensionJars = null;
    private final String[] applicationPackages = {"ffx",
        "javax.media.j3d",
        "javax.vecmath",
        "javax.media.opengl",
        "com.sun.j3d",
        "groovy",
        "org.codehaus.groovy",
        "org.apache.commons.cli",
        "org.apache.commons.configuration",
        "org.apache.commons.io",
        "org.apache.commons.lang",
        "org.apache.commons.lang3",
        "org.apache.commons.math3",
        "org.jogamp",
        "edu.rit.pj",
        "jcuda"};
    static final List<String> FFX_FILES;
    private static String gluegen = null;
    private static String jogl = null;
    private static String jocl = null;
    private static String nativeExtension = null;

    static {
        FFX_FILES = new ArrayList<>(Arrays.asList(new String[]{
            "edu.uiowa.eng.ffx/algorithms.jar",
            "edu.uiowa.eng.ffx/autoparm.jar",
            "edu.uiowa.eng.ffx/binding.jar",
            "edu.uiowa.eng.ffx/crystal.jar",
            "edu.uiowa.eng.ffx/numerics.jar",
            "edu.uiowa.eng.ffx/potential.jar",
            "edu.uiowa.eng.ffx/ui.jar",
            "edu.uiowa.eng.ffx/utilities.jar",
            "edu.uiowa.eng.ffx/xray.jar",
            // Force Field X Extensions
            "edu.uiowa.eng.ffx/algorithms-ext.jar",
            "edu.uiowa.eng.ffx/xray-ext.jar",
            // Groovy
            "org.codehaus.groovy/groovy-all.jar",
            // CUDA
            "jcuda/jcuda-all.jar",
            // Parallel Java
            "edu.rit.pj/pj.jar",
            // Java3D 1.6.2 (depends on JOGL v. 2.2.0)
            "java3d/j3dcore.jar",
            "java3d/j3dutils.jar",
            "java3d/j3dvrml.jar",
            "java3d/vecmath.jar",
            // JOGAMP GLUEGEN, JOGL and JOCL v. 2.2.0
            "org.jogamp.gluegen/gluegen-rt.jar",
            "org.jogamp.gluegen/gluegen-rt-main.jar",
            "org.jogamp.jogl/jogl-all.jar",
            "org.jogamp.jogl/jogl-all-main.jar",
            "org.jogamp.jocl/jocl.jar",
            "org.jogamp.jocl/jocl-main.jar",
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
            "commons-math/commons-math3.jar",
            // Mac OS X Extensions
            "macosx/AppleJavaExtensions.jar",
            // Java Help
            "javax.help/javahelp.jar",
            // BioJava
            "org.biojava/biojava3-core.jar",
            "org.biojava/core.jar",
            "org.biojava/bytecode.jar",
            "org.biojava/biojava3-structure.jar",
            "org.biojava/biojava3-alignment.jar",
            "org.biojava/biojava3-phylo.jar",
            // Lars Behnke's hierarchical-clustering-java
            "com.apporiented/hierarchical-clustering.jar"
        }));

        String osName = System.getProperty("os.name").toUpperCase();
        String osArch = System.getProperty("sun.arch.data.model");
        final boolean x8664 = "64".equals(osArch);

        if ("MAC OS X".equals(osName)) {
            nativeExtension = "-natives-macosx-universal.jar";
            // Gluegen Runtime Universal Binaries
            FFX_FILES.add("org.jogamp.gluegen/gluegen-rt-natives-macosx-universal.jar");
            // JOGL Universal Binaries
            FFX_FILES.add("org.jogamp.jogl/jogl-all-natives-macosx-universal.jar");
            // JOCL Universal Binaries
            FFX_FILES.add("org.jogamp.jocl/jocl-natives-macosx-universal.jar");
            if (x8664) {
                // JCUDA
                FFX_FILES.add("64-bit/libJCudaDriver-apple-x86_64.jnilib");
                FFX_FILES.add("64-bit/libJCudaRuntime-apple-x86_64.jnilib");
                FFX_FILES.add("64-bit/libJCufft-apple-x86_64.jnilib");
            }
        } else if ("LINUX".equals(osName)) {
            if (x8664) {
                nativeExtension = "-natives-linux-amd64.jar";
                FFX_FILES.add("64-bit/libJCufft-linux-x86_64.so");
                // Gluegen Runtime
                FFX_FILES.add("org.jogamp.gluegen/gluegen-rt-natives-linux-amd64.jar");
                // JOGL
                FFX_FILES.add("org.jogamp.jogl/jogl-all-natives-linux-amd64.jar");
                // JOCL
                FFX_FILES.add("org.jogamp.jocl/jocl-natives-linux-amd64.jar");
                // JCUDA
                FFX_FILES.add("64-bit/libJCudaDriver-linux-x86_64.so");
                FFX_FILES.add("64-bit/libJCudaRuntime-linux-x86_64.so");
            } else {
                nativeExtension = "-natives-linux-i586.jar";
                // Gluegen Runtime
                FFX_FILES.add("org.jogamp.gluegen/gluegen-rt-natives-linux-i586.jar");
                // JOGL
                FFX_FILES.add("org.jogamp.jogl/jogl-all-natives-linux-i586.jar");
                // JOCL
                FFX_FILES.add("org.jogamp.jocl/jocl-natives-linux-i586.jar");
            }
        } else if (osName.startsWith(
                "WINDOWS")) {
            if (x8664) {
                nativeExtension = "-natives-windows-amd64.jar";
                // Gluegen Runtime
                FFX_FILES.add("org.jogamp.glugen/gluegen-rt-natives-windows-amd64.jar");
                // JOGL
                FFX_FILES.add("org.jogamp.jogl/jogl-all-natives-windows-amd64.jar");
                // JOCL
                FFX_FILES.add("org.jogamp.jocl/jocl-natives-windows-amd64.jar");
                // JCUDA
                FFX_FILES.add("64-bit/JCudaDriver-linux-x86_64.dll");
                FFX_FILES.add("64-bit/JCudaRuntime-linux-x86_64.dll");
                FFX_FILES.add("64-bit/JCufft-linux-x86_64.dll");
            } else {
                nativeExtension = "-natives-windows-i586.jar";
                // Gluegen Runtime
                FFX_FILES.add("org.jogamp.gluegen/gluegen-rt-natives-windows-i586.jar");
                // JOGL
                FFX_FILES.add("org.jogamp.jogl/jogl-all-natives-windows-i586.jar");
                // JOCL
                FFX_FILES.add("org.jogamp.jocl/jocl-natives-windows-i586.jar");
            }
        }
    }

    /**
     * Force Field X custom class loader considers JARs and DLLs of
     * <code>extensionJarsAndDlls</code> as classpath and libclasspath elements
     * with a higher priority than the ones of default classpath. It will load
     * itself all the classes belonging to packages of
     * <code>applicationPackages</code>.
     *
     * @param parent a {@link java.lang.ClassLoader} object.
     */
    public FFXClassLoader(final ClassLoader parent) {
        super(new URL[0], parent);
        protectionDomain
                = FFXClassLoader.class
                .getProtectionDomain();
    }

    /**
     * Implementation of this method is to allow use of the NetBeans JFluid
     * profiler.
     *
     * @param value
     */
    private void appendToClassPathForInstrumentation(String value) {
        String tempDir = System.getProperty("java.io.tmpdir");
        File toDir = new File(tempDir + "deployed");
        toDir.mkdir();
        toDir = new File(tempDir + "deployed/jdk16");
        toDir.mkdir();
        toDir = new File(tempDir + "deployed/jdk16/mac");
        toDir.mkdir();

        String prof = tempDir + "deployed/jdk16/mac/libprofilerinterface.jnilib";
        File toFile = new File(prof);
        prof = "/Applications/NetBeans/NetBeans 8.0.app/Contents/Resources/NetBeans/profiler/lib/deployed/jdk16/mac/libprofilerinterface.jnilib";
        File fromFile = new File(prof);

        InputStream input = null;
        OutputStream output = null;
        try {
            input = new FileInputStream(fromFile);
            output = new BufferedOutputStream(new FileOutputStream(toFile));
            byte[] buffer = new byte[8192];
            int size;
            while ((size = input.read(buffer)) != -1) {
                output.write(buffer, 0, size);
            }
        } catch (Exception ex) {
            System.out.println(ex.toString());
        } finally {
            try {
                if (input != null) {
                    input.close();
                }
                if (output != null) {
                    output.close();
                }
            } catch (Exception e) {
                System.out.println(e.toString());
            }
        }
        toFile.deleteOnExit();
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     *
     * @param input a {@link java.io.InputStream} object.
     * @param name a {@link java.lang.String} object.
     * @param suffix a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public static String copyInputStreamToTmpFile(final InputStream input,
            String name, final String suffix) throws IOException {
        File tmpFile = null;

        if (name.contains("gluegen-rt") && name.contains("natives")) {
            tmpFile = new File(gluegen + nativeExtension);
        } else if (name.contains("jogl-all") && name.contains("natives")) {
            tmpFile = new File(jogl + nativeExtension);
        } else if (name.contains("jocl") && name.contains("natives")) {
            tmpFile = new File(jocl + nativeExtension);
        } else {
            try {
                name = "ffx." + name + ".";
                tmpFile = File.createTempFile(name, suffix);
            } catch (Exception e) {
                System.out.println(" Could not extract a Force Field X library.");
                System.err.println(e.toString());
                System.exit(-1);
            }
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
     * Finds and defines the given class among the extension JARs given in
     * constructor, then among resources.
     */
    @Override
    protected Class findClass(String name) throws ClassNotFoundException {

        /*
         if (name.startsWith("com.jogamp")) {
         System.out.println(" Class requested:" + name);
         } */
        if (!extensionsLoaded) {
            loadExtensions();
        }

        // Build class file from its name
        String classFile = name.replace('.', '/') + ".class";
        InputStream classInputStream = null;
        if (extensionJars != null) {
            for (JarFile extensionJar : extensionJars) {
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

        String path = (String) this.extensionDlls.get(libname);

        if (path == null) {
            path = super.findLibrary(libname);
        }

        return path;
    }

    /**
     * {@inheritDoc}
     *
     * Returns the URL of the given resource searching first if it exists among
     * the extension JARs given in constructor.
     */
    @Override
    public URL findResource(String name) {
        if (!extensionsLoaded) {
            loadExtensions();
        }

        if (name.equals("List all scripts")) {
            listScripts();
        }

        if (extensionJars != null) {
            for (JarFile extensionJar : extensionJars) {
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
        if (extensionJars != null) {

            List<String> scripts = new ArrayList<>();

            for (JarFile extensionJar : extensionJars) {
                Enumeration<JarEntry> entries = extensionJar.entries();
                while (entries.hasMoreElements()) {
                    JarEntry entry = entries.nextElement();
                    String name = entry.getName();
                    // System.out.println(name);
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
     * Loads a class with this class loader if its package belongs to
     * <code>applicationPackages</code> given in constructor.
     */
    @Override
    protected Class loadClass(String name, boolean resolve) throws ClassNotFoundException {
        if (!extensionsLoaded) {
            loadExtensions();
        }

        // If no extension jars were found use the super.loadClass method.
        if (extensionJars == null) {
            return super.loadClass(name, resolve);
        }

        // Check if the class has already been loaded
        Class loadedClass = findLoadedClass(name);
        if (loadedClass == null) {
            try {
                for (String applicationPackage : applicationPackages) {
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

        String extensionJarsAndDlls[] = FFX_FILES.toArray(new String[FFX_FILES.size()]);

        // Compute DLLs prefix and suffix
        String dllSuffix;
        String dllSuffix2 = null;
        String dllPrefix;

        String osName = System.getProperty("os.name");
        if (osName.startsWith("Windows")) {
            dllSuffix = ".dll";
            dllPrefix = "";
        } else if (osName.startsWith("Mac OS X")) {
            dllSuffix = ".jnilib";
            dllSuffix2 = ".dylib";
            dllPrefix = "lib";
        } else {
            dllSuffix = ".so";
            dllPrefix = "lib";
        }

        // Find extension Jars and DLLs
        ArrayList<JarFile> extensionJarList = new ArrayList<>();
        for (String extensionJarOrDll : extensionJarsAndDlls) {
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
                        if (name.equals("gluegen-rt")) {
                            gluegen = extensionJar.substring(0, extensionJar.length() - 4);
                        } else if (name.equals("jogl-all")) {
                            jogl = extensionJar.substring(0, extensionJar.length() - 4);
                        } else if (name.equals("jocl")) {
                            jocl = extensionJar.substring(0, extensionJar.length() - 4);
                        }
                    } else if (extensionJarOrDll.endsWith(dllSuffix)) {
                        int start = lastSlashIndex + 1 + dllPrefix.length();
                        int end = extensionJarOrDll.indexOf(dllSuffix);
                        String name = extensionJarOrDll.substring(start, end);
                        // Copy DLL to a tmp file
                        String extensionDll = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(),
                                name, dllSuffix);
                        // Add extracted file to extension DLLs map
                        extensionDlls.put(name, extensionDll);
                    } else if (dllSuffix2 != null && extensionJarOrDll.endsWith(dllSuffix2)) {
                        int start = lastSlashIndex + 1 + dllPrefix.length();
                        int end = extensionJarOrDll.indexOf(dllSuffix2);
                        String name = extensionJarOrDll.substring(start, end);
                        // Copy DLL to a tmp file
                        String extensionDll = copyInputStreamToTmpFile(extensionJarOrDllUrl.openStream(),
                                name, dllSuffix2);
                        // Add extracted file to extension DLLs map
                        extensionDlls.put(name, extensionDll);
                    }
                }
            } catch (IOException ex) {
                System.out.println(extensionJarOrDll);
                throw new RuntimeException(" Couldn't extract extension jar.\n", ex);
            }
        }

        // Create extensionJars array
        if (extensionJarList.size() > 0) {
            extensionJars = (JarFile[]) extensionJarList.toArray(new JarFile[extensionJarList.size()]);
        }
    }
}
