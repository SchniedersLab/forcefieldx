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
package ffx.openmm.ffm.bindings;

import java.lang.foreign.AddressLayout;
import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.GroupLayout;
import java.lang.foreign.Linker;
import java.lang.foreign.MemoryLayout;
import java.lang.foreign.PaddingLayout;
import java.lang.foreign.SequenceLayout;
import java.lang.foreign.StructLayout;
import java.lang.foreign.ValueLayout;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.util.Arrays;
import java.util.stream.Collectors;

public class OpenMMNative$shared {

  OpenMMNative$shared() {
    // Should not be called directly
  }

  public static final ValueLayout.OfBoolean C_BOOL = (ValueLayout.OfBoolean) Linker.nativeLinker().canonicalLayouts().get("bool");
  public static final ValueLayout.OfByte C_CHAR = (ValueLayout.OfByte) Linker.nativeLinker().canonicalLayouts().get("char");
  public static final ValueLayout.OfShort C_SHORT = (ValueLayout.OfShort) Linker.nativeLinker().canonicalLayouts().get("short");
  public static final ValueLayout.OfInt C_INT = (ValueLayout.OfInt) Linker.nativeLinker().canonicalLayouts().get("int");
  public static final ValueLayout.OfLong C_LONG_LONG = (ValueLayout.OfLong) Linker.nativeLinker().canonicalLayouts().get("long long");
  public static final ValueLayout.OfFloat C_FLOAT = (ValueLayout.OfFloat) Linker.nativeLinker().canonicalLayouts().get("float");
  public static final ValueLayout.OfDouble C_DOUBLE = (ValueLayout.OfDouble) Linker.nativeLinker().canonicalLayouts().get("double");
  public static final AddressLayout C_POINTER = ((AddressLayout) Linker.nativeLinker().canonicalLayouts().get("void*"))
      .withTargetLayout(MemoryLayout.sequenceLayout(java.lang.Long.MAX_VALUE, C_CHAR));
  public static final ValueLayout.OfLong C_LONG = (ValueLayout.OfLong) Linker.nativeLinker().canonicalLayouts().get("long");

  static final boolean TRACE_DOWNCALLS = Boolean.getBoolean("jextract.trace.downcalls");

  static void traceDowncall(String name, Object... args) {
    String traceArgs = Arrays.stream(args)
        .map(Object::toString)
        .collect(Collectors.joining(", "));
    System.out.printf("%s(%s)\n", name, traceArgs);
  }

  static MethodHandle upcallHandle(Class<?> fi, String name, FunctionDescriptor fdesc) {
    try {
      return MethodHandles.lookup().findVirtual(fi, name, fdesc.toMethodType());
    } catch (ReflectiveOperationException ex) {
      throw new AssertionError(ex);
    }
  }

  static MemoryLayout align(MemoryLayout layout, long align) {
    return switch (layout) {
      case PaddingLayout p -> p;
      case ValueLayout v -> v.withByteAlignment(align);
      case GroupLayout g -> {
        MemoryLayout[] alignedMembers = g.memberLayouts().stream()
            .map(m -> align(m, align)).toArray(MemoryLayout[]::new);
        yield g instanceof StructLayout ?
            MemoryLayout.structLayout(alignedMembers) : MemoryLayout.unionLayout(alignedMembers);
      }
      case SequenceLayout s -> MemoryLayout.sequenceLayout(s.elementCount(), align(s.elementLayout(), align));
    };
  }
}

