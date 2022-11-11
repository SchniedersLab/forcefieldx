// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import java.util.Set;
import javax.annotation.processing.AbstractProcessor;
import javax.annotation.processing.RoundEnvironment;
import javax.annotation.processing.SupportedAnnotationTypes;
import javax.annotation.processing.SupportedOptions;
import javax.lang.model.element.Element;
import javax.lang.model.element.TypeElement;
import javax.lang.model.type.MirroredTypeException;
import javax.lang.model.type.TypeMirror;


/**
 * Log out FFXKeyword Annotations for documentation purposes.
 */
@SupportedAnnotationTypes("ffx.utilities.FFXKeyword")
@SupportedOptions({"keywordDir"})
public class KeywordProcessor extends AbstractProcessor {

  /**
   * {@inheritdoc}
   */
  @Override
  public boolean process(Set<? extends TypeElement> annotations, RoundEnvironment roundEnv) {
    if (annotations.size() == 0) {
      return true;
    }

    String keywordDir = processingEnv.getOptions().get("keywordDir");
    // System.out.println(" KeywordDir: " + keywordDir);

    // Create the Keyword directory
    Path keyPath = Paths.get(keywordDir);
    if (!Files.exists(keyPath)) {
      try {
        Files.createDirectory(keyPath);
      } catch (IOException exception) {
        System.out.println(" Exception creating keyword directory:\n " + exception);
      }
    }

    // TypeElement[] typeElements = annotations.toArray(new TypeElement[0]);
    // System.out.println(" Total Annotations Types: " + typeElements.length);
    // System.out.println(" First Annotation Type:   " + typeElements[0]);

    Set<? extends Element> annotatedKeywords = roundEnv.getElementsAnnotatedWith(FFXKeyword.class);
    // System.out.println(" FFXKeyword Annotations Processed: " + annotatedKeywords.size());
    // Loop over FFXKeyword annotations.
    for (Element element : annotatedKeywords) {
      FFXKeyword ffxKeyword = element.getAnnotation(FFXKeyword.class);
      StringBuilder sb = new StringBuilder(format("\n=== %s\n", ffxKeyword.name()));
      sb.append(format("  Definition:   %s\n", ffxKeyword.description()));
      // This causes a MirroredTypeException, which has a TypeMirror whose value is the class.
      try {
        Class clazz = ffxKeyword.clazz();
      } catch (MirroredTypeException e) {
        TypeMirror typeMirror = e.getTypeMirror();
        String type = typeMirror.toString();
        type = type.replace("java.lang.", "");
        sb.append(format("  Type:         %s\n", type));
      }
      String defaultValue = ffxKeyword.defaultValue();
      if (defaultValue != null && !defaultValue.equals("")) {
        sb.append(format("  Default:      %s\n", ffxKeyword.defaultValue()));
      }

      // Create the Keyword Group sub-directory
      KeywordGroup group = ffxKeyword.keywordGroup();
      Path groupPath = Paths.get(keyPath.toString(), group.name());

      if (!Files.exists(groupPath)) {
        try {
          Files.createDirectory(groupPath);
        } catch (IOException exception) {
          System.out.println(" Exception creating keyword group sub-directory:\n " + exception);
        }
      }

      // Create or update the file for this Keyword.
      Path keywordPath = Paths.get(groupPath.toString(), ffxKeyword.name() + ".adoc");
      try {
        Files.writeString(keywordPath, sb.toString(), CREATE, WRITE, TRUNCATE_EXISTING);
      } catch (Exception e) {
        System.out.println(" Exception writing keyword:\n " + e);
      }

      // System.out.println(sb);
    }

    return true;
  }
}
