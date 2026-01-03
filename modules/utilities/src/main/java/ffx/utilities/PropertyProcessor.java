// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.utilities;

import static java.lang.String.format;
import static java.nio.file.StandardOpenOption.CREATE;
import static java.nio.file.StandardOpenOption.TRUNCATE_EXISTING;
import static java.nio.file.StandardOpenOption.WRITE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import javax.annotation.processing.*;
import javax.lang.model.SourceVersion;
import javax.lang.model.element.Element;
import javax.lang.model.element.TypeElement;
import javax.lang.model.type.MirroredTypeException;
import javax.lang.model.type.TypeMirror;


/**
 * Log out FFXProperty Annotations for documentation purposes.
 *
 * @author Michael J. Schnieders
 */
@SupportedAnnotationTypes({"ffx.utilities.FFXProperty", "ffx.utilities.FFXProperties"})
@SupportedOptions({"propertyDir"})
@SupportedSourceVersion(SourceVersion.RELEASE_25)
public class PropertyProcessor extends AbstractProcessor {

  /**
   * Default constructor.
   */
  public PropertyProcessor() {
    super();
  }

  /**
   * Processes a set of annotation types on type elements originating from the prior round and
   * returns whether these annotations are claimed by this processor. If true is returned, the
   * annotations are claimed and subsequent processors will not be asked to process them; if false is
   * returned, the annotations are unclaimed and subsequent processors may be asked to process them.
   * A processor may always return the same boolean value or may vary the result based on chosen
   * criteria. The input set will be empty if the processor supports "*" and the root elements have
   * no annotations. A Processor must gracefully handle an empty set of annotations.
   */
  @Override
  public boolean process(Set<? extends TypeElement> annotations, RoundEnvironment roundEnv) {
    if (annotations.isEmpty()) {
      return true;
    }

    String propertyDir = processingEnv.getOptions().get("propertyDir");
    // System.out.println(" PropertyDir: " + propertyDir);

    // Create the Property directory
    Path propertyPath = Paths.get(propertyDir);
    if (!Files.exists(propertyPath)) {
      try {
        Files.createDirectory(propertyPath);
      } catch (IOException exception) {
        System.out.println(" Exception creating Property directory:\n " + exception);
      }
    }

    // TypeElement[] typeElements = annotations.toArray(new TypeElement[0]);
    // System.out.println(" Total Annotations Types: " + typeElements.length);
    // System.out.println(" First Annotation Type:   " + typeElements[0]);

    // Collect fields and classes annotated once.
    List<FFXProperty> propertyList = new ArrayList<>();
    Set<? extends Element> annotatedProperties = roundEnv.getElementsAnnotatedWith(FFXProperty.class);
    for (Element element : annotatedProperties) {
      FFXProperty ffxProperty = element.getAnnotation(FFXProperty.class);
      propertyList.add(ffxProperty);
    }
    // Collect classes annotated more than once.
    Set<? extends Element> propertyArrays = roundEnv.getElementsAnnotatedWith(FFXProperties.class);
    for (Element element : propertyArrays) {
      FFXProperties ffxProperties = element.getAnnotation(FFXProperties.class);
      FFXProperty[] properties = ffxProperties.value();
      Collections.addAll(propertyList, properties);
    }

    // System.out.println(" FFXProperty Annotations Processed: " + propertyList.size());
    // Loop over FFXProperty annotations.
    for (FFXProperty ffxProperty : propertyList) {
      StringBuilder sb = new StringBuilder(format("\n=== %s\n", ffxProperty.name()));
      sb.append("[%collapsible]\n====\n");
      // This causes a MirroredTypeException, which has a TypeMirror whose value is the class.
      try {
        Class<?> clazz = ffxProperty.clazz();
      } catch (MirroredTypeException e) {
        TypeMirror typeMirror = e.getTypeMirror();
        String type = typeMirror.toString();
        type = type.replace("java.lang.", "");
        sb.append(format("  Type:         %s\n", type));
      }
      String defaultValue = ffxProperty.defaultValue();
      if (defaultValue != null && !defaultValue.isEmpty()) {
        sb.append(format("  Default:      %s\n", ffxProperty.defaultValue()));
      }
      sb.append(format("  Definition:   %s\n", ffxProperty.description()));
      sb.append("====\n");

      // Create the Property Group subdirectory.
      PropertyGroup group = ffxProperty.propertyGroup();
      Path groupPath = Paths.get(propertyPath.toString(), group.name());

      if (!Files.exists(groupPath)) {
        try {
          Files.createDirectory(groupPath);
        } catch (IOException exception) {
          System.out.println(" Exception creating property group sub-directory:\n " + exception);
        }
      }

      // Create or update the file for this Property.
      Path adocPath = Paths.get(groupPath.toString(), ffxProperty.name() + ".adoc");
      try {
        Files.writeString(adocPath, sb.toString(), CREATE, WRITE, TRUNCATE_EXISTING);
      } catch (Exception e) {
        System.out.println(" Exception writing property:\n " + e);
      }
      // System.out.println(sb);
    }

    return true;
  }
}
