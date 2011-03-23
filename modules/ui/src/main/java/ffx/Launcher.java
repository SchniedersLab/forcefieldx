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

import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This bootstrap class loads Force Field X application classes from jars
 * in classpath or from extension jars stored as resources.
 *
 * @author Michael J. Schnieders<br>
 *         Derived from a similar class by Emmanuel Puybaret
 *
 * @since 1.0
 */
public class Launcher {

    public static void main(String[] args) throws MalformedURLException, IllegalAccessException,
            InvocationTargetException, NoSuchMethodException, ClassNotFoundException {
        Class ffxBootstrapClass = Launcher.class;

        List<String> ffxFiles = new ArrayList<String>(Arrays.asList(new String[]{
                    "com.kenai.ffx/algorithms.jar",
                    "com.kenai.ffx/potentials.jar",
                    "com.kenai.ffx/crystal.jar",
                    "com.kenai.ffx/numerics.jar",
                    "com.kenai.ffx/ui.jar",
                    "com.kenai.ffx/utilities.jar",
                    "com.kenai.ffx/xray.jar",
                    "com.kenai.ffx/autoparm.jar",
                    "org.codehaus.groovy/groovy-all.jar",
                    "edu.rit.pj/pj.jar",
                    "jcuda/jcuda-all.jar",
                    "java3d/j3dcore.jar",
                    "java3d/j3dutils.jar",
                    "java3d/j3dvrml.jar",
                    "java3d/vecmath.jar",
                    "net.java.dev.jogl/jogl.jar",
                    "net.java.dev.jogl/gluegen-rt.jar",
                    "net.java.dev.jogl/jocl.jar",
                    "commons-beanutils/commons-beanutils.jar",
                    "commons-collections/commons-collections.jar",
                    "commons-configuration/commons-configuration.jar",
                    "commons-digester/commons-digester.jar",
                    "commons-io/commons-io.jar",
                    "commons-lang/commons-lang.jar",
                    "commons-logging/commons-logging.jar",
                    "commons-math/commons-math.jar",
                    "org.apache.ant/ant.jar",
                    "org.apache.ant/ant-launcher.jar",
                    "macosx/AppleJavaExtensions.jar",
                    "javax.help/javahelp.jar",
                    "junit/junit.jar",
                    "jline/jline.jar"
                }));
        String osName = System.getProperty("os.name").toUpperCase();
        String osArch = System.getProperty("os.arch").toUpperCase();
        boolean x86_64 = "64".equals(System.getProperty("sun.arch.data.model"));
        if ("MAC OS X".equals(osName)) {
            // Mac OS X JOGL Universal Binaries
            ffxFiles.add("jogl/libgluegen-rt.jnilib");
            ffxFiles.add("jogl/libjogl.jnilib");
            ffxFiles.add("jogl/libjogl_awt.jnilib");
            ffxFiles.add("jogl/libjogl_cg.jnilib");
            if (x86_64) {
                // Mac OS X Jcuda 64-bit Binaries
                ffxFiles.add("jcuda-64/libJCudaDriver-apple-x86_64.jnilib");
                ffxFiles.add("jcuda-64/libJCudaRuntime-apple-x86_64.jnilib");
                ffxFiles.add("jcuda-64/libJCufft-apple-x86_64.jnilib");
            }
        } else if ("LINUX".equals(osName)) {
            if (x86_64) {
                ffxFiles.add("jogl-64/libgluegen-rt.so");
                ffxFiles.add("jogl-64/libjogl.so");
                ffxFiles.add("jogl-64/libjogl_awt.so");
                ffxFiles.add("jogl-64/libjogl_cg.so");
                if ("X86_64".equals(osArch)) {
                    ffxFiles.add("jcuda-64/libJCudaDriver-apple-x86_64.so");
                    ffxFiles.add("jcuda-64/libJCudaRuntime-apple-x86_64.so");
                    ffxFiles.add("jcuda-64/libJCufft-apple-x86_64.so");
                } else if ("AMD64".equals(osArch)) {
                    ffxFiles.add("jcuda-64/libJCudaDriver-apple-amd64.so");
                    ffxFiles.add("jcuda-64/libJCudaRuntime-apple-amd64.so");
                    ffxFiles.add("jcuda-64/libJCufft-apple-amd64.so");
                }
            } else {
                ffxFiles.add("jogl-32/libgluegen-rt.so");
                ffxFiles.add("jogl-32/libjogl.so");
                ffxFiles.add("jogl-32/libjogl_awt.so");
                ffxFiles.add("jogl-32/libjogl_cg.so");
            }
        } else if (osName.startsWith("WINDOWS")) {
            if (x86_64) {
                ffxFiles.add("jogl-64/gluegen-rt.dll");
                ffxFiles.add("jogl-64/jogl.dll");
                ffxFiles.add("jogl-64/jogl_awt.dll");
                ffxFiles.add("jogl-64/jogl_cg.dll");
                ffxFiles.add("jcuda-64/JCudaDriver-apple-x86_64.dll");
                ffxFiles.add("jcuda-64/JCudaRuntime-apple-x86_64.dll");
                ffxFiles.add("jcuda-64/JCufft-apple-x86_64.dll");
            } else {
                ffxFiles.add("jogl-32/gluegen-rt.dll");
                ffxFiles.add("jogl-32/jogl.dll");
                ffxFiles.add("jogl-32/jogl_awt.dll");
                ffxFiles.add("jogl-32/jogl_cg.dll");
            }
        }

        String[] applicationPackages = {"ffx",
            // Java 3D packages
            "javax.media.j3d",
            "javax.vecmath",
            "com.sun.j3d",
            "com.sun.opengl",
            "com.sun.gluegen.runtime",
            "javax.media.opengl",
            "groovy",
            "org.codehaus.groovy",
            "org.apache.commons.configuration",
            "org.apache.commons.io",
            "org.apache.commons.lang",
            "org.apache.commons.math",
            "edu.rit.pj",
            "jcuda"};
        ClassLoader classLoader = new FFXClassLoader(
                ffxBootstrapClass.getClassLoader(),
                ffxBootstrapClass.getProtectionDomain(),
                ffxFiles.toArray(new String[ffxFiles.size()]), applicationPackages);

        String applicationClassName = "ffx.Main";
        Class applicationClass = classLoader.loadClass(applicationClassName);
        Method applicationClassMain =
                applicationClass.getMethod("main", Array.newInstance(String.class, 0).getClass());
        // Call application class main method with reflection
        applicationClassMain.invoke(null, new Object[]{args});
    }
    private String packages[] = {
        "ffx.algorithms",
        "ffx.crystal",
        "ffx.numerics",
        "ffx.numerics.fft",
        "ffx.potential",
        "ffx.potential.bonded",
        "ffx.potential.nonbonded",
        "ffx.potential.parameters",
        "ffx.potential.parameters.amoeba",
        "ffx.potential.parsers",
        "ffx.potential.structures",
        "ffx.ui",
        "ffx.ui.behaviors",
        "ffx.ui.commands",
        "ffx.ui.icons",
        "ffx.ui.macosx",
        "ffx.ui.properties",
        "ffx.utilities",
        "ffx.xray",
        "ffx.xray.structures"};
}
