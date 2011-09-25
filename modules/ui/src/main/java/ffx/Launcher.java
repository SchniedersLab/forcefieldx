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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.rit.pj.Comm;

/**
 * This bootstrap class loads Force Field X application classes from jars
 * in classpath or from extension jars stored as resources.
 *
 * @author Michael J. Schnieders<br>
 *         Derived from a similar class by Emmanuel Puybaret
 * @since 1.0
 * @version $Id: $
 */
public class Launcher {
    
    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.ClassNotFoundException if any.
     * @throws java.lang.NoSuchMethodException if any.
     * @throws java.lang.IllegalAccessException if any.
     * @throws java.lang.IllegalArgumentException if any.
     * @throws java.lang.reflect.InvocationTargetException if any.
     */
    public static void main(String[] args) throws ClassNotFoundException, NoSuchMethodException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {

        List<String> ffxFiles = new ArrayList<String>(Arrays.asList(new String[]{
                    "com.kenai.ffx/algorithms.jar",
                    "com.kenai.ffx/autoparm.jar",
                    "com.kenai.ffx/binding.jar",
                    "com.kenai.ffx/crystal.jar",
                    "com.kenai.ffx/numerics.jar",
                    "com.kenai.ffx/potentials.jar",
                    "com.kenai.ffx/ui.jar",
                    "com.kenai.ffx/utilities.jar",
                    "com.kenai.ffx/xray.jar",
                    "org.codehaus.groovy/groovy-all.jar",
                    "jcuda/jcuda-all.jar",
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
                    "commons-beanutils/commons-beanutils.jar",
                    "commons-collections/commons-collections.jar",
                    "commons-configuration/commons-configuration.jar",
                    "commons-digester/commons-digester.jar",
                    "commons-io/commons-io.jar",
                    "commons-lang/commons-lang.jar",
                    "commons-logging/commons-logging.jar",
                    "commons-math/commons-math.jar",
                    "macosx/AppleJavaExtensions.jar",
                    "javax.help/javahelp.jar"
                }));

        String osName = System.getProperty("os.name").toUpperCase();
        String osArch = System.getProperty("os.arch").toUpperCase();
        boolean x86_64 = "64".equals(System.getProperty("sun.arch.data.model"));
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
                if (osArch.equalsIgnoreCase("x86_64")) {
                    ffxFiles.add("64-bit/libJCudaDriver-linux-x86_64.so");
                    ffxFiles.add("64-bit/libJCudaRuntime-linux-x86_64.so");
                    ffxFiles.add("64-bit/libJCufft-linux-x86_64.so");
                } else if (osArch.equalsIgnoreCase("amd64")) {
                    ffxFiles.add("64-bit/libJCudaDriver-linux-amd64.so");
                    ffxFiles.add("64-bit/libJCudaRuntime-linux-amd64.so");
                    ffxFiles.add("64-bit/libJCufft-linux-amd64.so");
                }
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

        Class bootStrap = Launcher.class;
        ClassLoader classLoader = new FFXClassLoader(bootStrap.getClassLoader());
        
        Class ffxMain = classLoader.loadClass("ffx.Main");
        Method main = ffxMain.getMethod("main", Array.newInstance(String.class, 0).getClass());
        main.invoke(null, new Object[]{args});
        
    }
}
