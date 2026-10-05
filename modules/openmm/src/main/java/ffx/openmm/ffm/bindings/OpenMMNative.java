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

import java.lang.foreign.FunctionDescriptor;
import java.lang.foreign.Linker;
import java.lang.foreign.MemoryLayout;
import java.lang.foreign.MemorySegment;
import java.lang.invoke.MethodHandle;

public class OpenMMNative extends OpenMMNative_1 {

  OpenMMNative() {
    // Should not be called directly
  }

  private static class OpenMM_RMSDForce_getReferencePositions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_getReferencePositions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_Vec3Array *OpenMM_RMSDForce_getReferencePositions(const OpenMM_RMSDForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_getReferencePositions$descriptor() {
    return OpenMM_RMSDForce_getReferencePositions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_Vec3Array *OpenMM_RMSDForce_getReferencePositions(const OpenMM_RMSDForce *target)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_getReferencePositions$handle() {
    return OpenMM_RMSDForce_getReferencePositions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_Vec3Array *OpenMM_RMSDForce_getReferencePositions(const OpenMM_RMSDForce *target)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_getReferencePositions$address() {
    return OpenMM_RMSDForce_getReferencePositions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_Vec3Array *OpenMM_RMSDForce_getReferencePositions(const OpenMM_RMSDForce *target)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_getReferencePositions(MemorySegment target) {
    var mh$ = OpenMM_RMSDForce_getReferencePositions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_getReferencePositions", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_RMSDForce_setReferencePositions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_setReferencePositions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setReferencePositions(OpenMM_RMSDForce *target, const OpenMM_Vec3Array *positions)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_setReferencePositions$descriptor() {
    return OpenMM_RMSDForce_setReferencePositions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setReferencePositions(OpenMM_RMSDForce *target, const OpenMM_Vec3Array *positions)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_setReferencePositions$handle() {
    return OpenMM_RMSDForce_setReferencePositions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setReferencePositions(OpenMM_RMSDForce *target, const OpenMM_Vec3Array *positions)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_setReferencePositions$address() {
    return OpenMM_RMSDForce_setReferencePositions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setReferencePositions(OpenMM_RMSDForce *target, const OpenMM_Vec3Array *positions)
   *}
   */
  public static void OpenMM_RMSDForce_setReferencePositions(MemorySegment target, MemorySegment positions) {
    var mh$ = OpenMM_RMSDForce_setReferencePositions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_setReferencePositions", target, positions);
      }
      mh$.invokeExact(target, positions);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_RMSDForce_getParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_getParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_IntArray *OpenMM_RMSDForce_getParticles(const OpenMM_RMSDForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_getParticles$descriptor() {
    return OpenMM_RMSDForce_getParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_IntArray *OpenMM_RMSDForce_getParticles(const OpenMM_RMSDForce *target)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_getParticles$handle() {
    return OpenMM_RMSDForce_getParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_IntArray *OpenMM_RMSDForce_getParticles(const OpenMM_RMSDForce *target)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_getParticles$address() {
    return OpenMM_RMSDForce_getParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_IntArray *OpenMM_RMSDForce_getParticles(const OpenMM_RMSDForce *target)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_getParticles(MemorySegment target) {
    var mh$ = OpenMM_RMSDForce_getParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_getParticles", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_RMSDForce_setParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_setParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setParticles(OpenMM_RMSDForce *target, const OpenMM_IntArray *particles)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_setParticles$descriptor() {
    return OpenMM_RMSDForce_setParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setParticles(OpenMM_RMSDForce *target, const OpenMM_IntArray *particles)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_setParticles$handle() {
    return OpenMM_RMSDForce_setParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setParticles(OpenMM_RMSDForce *target, const OpenMM_IntArray *particles)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_setParticles$address() {
    return OpenMM_RMSDForce_setParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_setParticles(OpenMM_RMSDForce *target, const OpenMM_IntArray *particles)
   *}
   */
  public static void OpenMM_RMSDForce_setParticles(MemorySegment target, MemorySegment particles) {
    var mh$ = OpenMM_RMSDForce_setParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_setParticles", target, particles);
      }
      mh$.invokeExact(target, particles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_RMSDForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_updateParametersInContext(OpenMM_RMSDForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_updateParametersInContext$descriptor() {
    return OpenMM_RMSDForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_updateParametersInContext(OpenMM_RMSDForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_updateParametersInContext$handle() {
    return OpenMM_RMSDForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_updateParametersInContext(OpenMM_RMSDForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_updateParametersInContext$address() {
    return OpenMM_RMSDForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_RMSDForce_updateParametersInContext(OpenMM_RMSDForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_RMSDForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_RMSDForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_RMSDForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_RMSDForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_RMSDForce_usesPeriodicBoundaryConditions(const OpenMM_RMSDForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_RMSDForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_RMSDForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_RMSDForce_usesPeriodicBoundaryConditions(const OpenMM_RMSDForce *target)
   *}
   */
  public static MethodHandle OpenMM_RMSDForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_RMSDForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_RMSDForce_usesPeriodicBoundaryConditions(const OpenMM_RMSDForce *target)
   *}
   */
  public static MemorySegment OpenMM_RMSDForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_RMSDForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_RMSDForce_usesPeriodicBoundaryConditions(const OpenMM_RMSDForce *target)
   *}
   */
  public static int OpenMM_RMSDForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_RMSDForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_RMSDForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_CustomExternalForce *OpenMM_CustomExternalForce_create(const char *energy)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_create$descriptor() {
    return OpenMM_CustomExternalForce_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_CustomExternalForce *OpenMM_CustomExternalForce_create(const char *energy)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_create$handle() {
    return OpenMM_CustomExternalForce_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_CustomExternalForce *OpenMM_CustomExternalForce_create(const char *energy)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_create$address() {
    return OpenMM_CustomExternalForce_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_CustomExternalForce *OpenMM_CustomExternalForce_create(const char *energy)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_create(MemorySegment energy) {
    var mh$ = OpenMM_CustomExternalForce_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_create", energy);
      }
      return (MemorySegment) mh$.invokeExact(energy);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_destroy(OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_destroy$descriptor() {
    return OpenMM_CustomExternalForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_destroy(OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_destroy$handle() {
    return OpenMM_CustomExternalForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_destroy(OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_destroy$address() {
    return OpenMM_CustomExternalForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_destroy(OpenMM_CustomExternalForce *target)
   *}
   */
  public static void OpenMM_CustomExternalForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumParticles(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getNumParticles$descriptor() {
    return OpenMM_CustomExternalForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumParticles(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getNumParticles$handle() {
    return OpenMM_CustomExternalForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumParticles(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getNumParticles$address() {
    return OpenMM_CustomExternalForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumParticles(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static int OpenMM_CustomExternalForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getNumPerParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getNumPerParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumPerParticleParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getNumPerParticleParameters$descriptor() {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumPerParticleParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getNumPerParticleParameters$handle() {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumPerParticleParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getNumPerParticleParameters$address() {
    return OpenMM_CustomExternalForce_getNumPerParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumPerParticleParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static int OpenMM_CustomExternalForce_getNumPerParticleParameters(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_getNumPerParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getNumPerParticleParameters", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getNumGlobalParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getNumGlobalParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumGlobalParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getNumGlobalParameters$descriptor() {
    return OpenMM_CustomExternalForce_getNumGlobalParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumGlobalParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getNumGlobalParameters$handle() {
    return OpenMM_CustomExternalForce_getNumGlobalParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumGlobalParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getNumGlobalParameters$address() {
    return OpenMM_CustomExternalForce_getNumGlobalParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_getNumGlobalParameters(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static int OpenMM_CustomExternalForce_getNumGlobalParameters(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_getNumGlobalParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getNumGlobalParameters", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getEnergyFunction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getEnergyFunction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getEnergyFunction(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getEnergyFunction$descriptor() {
    return OpenMM_CustomExternalForce_getEnergyFunction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getEnergyFunction(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getEnergyFunction$handle() {
    return OpenMM_CustomExternalForce_getEnergyFunction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getEnergyFunction(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getEnergyFunction$address() {
    return OpenMM_CustomExternalForce_getEnergyFunction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getEnergyFunction(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getEnergyFunction(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_getEnergyFunction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getEnergyFunction", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_setEnergyFunction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_setEnergyFunction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setEnergyFunction(OpenMM_CustomExternalForce *target, const char *energy)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_setEnergyFunction$descriptor() {
    return OpenMM_CustomExternalForce_setEnergyFunction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setEnergyFunction(OpenMM_CustomExternalForce *target, const char *energy)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_setEnergyFunction$handle() {
    return OpenMM_CustomExternalForce_setEnergyFunction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setEnergyFunction(OpenMM_CustomExternalForce *target, const char *energy)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_setEnergyFunction$address() {
    return OpenMM_CustomExternalForce_setEnergyFunction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setEnergyFunction(OpenMM_CustomExternalForce *target, const char *energy)
   *}
   */
  public static void OpenMM_CustomExternalForce_setEnergyFunction(MemorySegment target, MemorySegment energy) {
    var mh$ = OpenMM_CustomExternalForce_setEnergyFunction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_setEnergyFunction", target, energy);
      }
      mh$.invokeExact(target, energy);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_addPerParticleParameter {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_addPerParticleParameter");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addPerParticleParameter(OpenMM_CustomExternalForce *target, const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_addPerParticleParameter$descriptor() {
    return OpenMM_CustomExternalForce_addPerParticleParameter.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addPerParticleParameter(OpenMM_CustomExternalForce *target, const char *name)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_addPerParticleParameter$handle() {
    return OpenMM_CustomExternalForce_addPerParticleParameter.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addPerParticleParameter(OpenMM_CustomExternalForce *target, const char *name)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_addPerParticleParameter$address() {
    return OpenMM_CustomExternalForce_addPerParticleParameter.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addPerParticleParameter(OpenMM_CustomExternalForce *target, const char *name)
   *}
   */
  public static int OpenMM_CustomExternalForce_addPerParticleParameter(MemorySegment target, MemorySegment name) {
    var mh$ = OpenMM_CustomExternalForce_addPerParticleParameter.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_addPerParticleParameter", target, name);
      }
      return (int) mh$.invokeExact(target, name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getPerParticleParameterName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getPerParticleParameterName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getPerParticleParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getPerParticleParameterName$descriptor() {
    return OpenMM_CustomExternalForce_getPerParticleParameterName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getPerParticleParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getPerParticleParameterName$handle() {
    return OpenMM_CustomExternalForce_getPerParticleParameterName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getPerParticleParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getPerParticleParameterName$address() {
    return OpenMM_CustomExternalForce_getPerParticleParameterName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getPerParticleParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getPerParticleParameterName(MemorySegment target, int index) {
    var mh$ = OpenMM_CustomExternalForce_getPerParticleParameterName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getPerParticleParameterName", target, index);
      }
      return (MemorySegment) mh$.invokeExact(target, index);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_setPerParticleParameterName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_setPerParticleParameterName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setPerParticleParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_setPerParticleParameterName$descriptor() {
    return OpenMM_CustomExternalForce_setPerParticleParameterName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setPerParticleParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_setPerParticleParameterName$handle() {
    return OpenMM_CustomExternalForce_setPerParticleParameterName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setPerParticleParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_setPerParticleParameterName$address() {
    return OpenMM_CustomExternalForce_setPerParticleParameterName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setPerParticleParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static void OpenMM_CustomExternalForce_setPerParticleParameterName(MemorySegment target, int index, MemorySegment name) {
    var mh$ = OpenMM_CustomExternalForce_setPerParticleParameterName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_setPerParticleParameterName", target, index, name);
      }
      mh$.invokeExact(target, index, name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_addGlobalParameter {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_addGlobalParameter");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addGlobalParameter(OpenMM_CustomExternalForce *target, const char *name, double defaultValue)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_addGlobalParameter$descriptor() {
    return OpenMM_CustomExternalForce_addGlobalParameter.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addGlobalParameter(OpenMM_CustomExternalForce *target, const char *name, double defaultValue)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_addGlobalParameter$handle() {
    return OpenMM_CustomExternalForce_addGlobalParameter.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addGlobalParameter(OpenMM_CustomExternalForce *target, const char *name, double defaultValue)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_addGlobalParameter$address() {
    return OpenMM_CustomExternalForce_addGlobalParameter.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addGlobalParameter(OpenMM_CustomExternalForce *target, const char *name, double defaultValue)
   *}
   */
  public static int OpenMM_CustomExternalForce_addGlobalParameter(MemorySegment target, MemorySegment name, double defaultValue) {
    var mh$ = OpenMM_CustomExternalForce_addGlobalParameter.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_addGlobalParameter", target, name, defaultValue);
      }
      return (int) mh$.invokeExact(target, name, defaultValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getGlobalParameterName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getGlobalParameterName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getGlobalParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getGlobalParameterName$descriptor() {
    return OpenMM_CustomExternalForce_getGlobalParameterName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getGlobalParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getGlobalParameterName$handle() {
    return OpenMM_CustomExternalForce_getGlobalParameterName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getGlobalParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getGlobalParameterName$address() {
    return OpenMM_CustomExternalForce_getGlobalParameterName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_CustomExternalForce_getGlobalParameterName(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getGlobalParameterName(MemorySegment target, int index) {
    var mh$ = OpenMM_CustomExternalForce_getGlobalParameterName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getGlobalParameterName", target, index);
      }
      return (MemorySegment) mh$.invokeExact(target, index);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_setGlobalParameterName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_setGlobalParameterName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_setGlobalParameterName$descriptor() {
    return OpenMM_CustomExternalForce_setGlobalParameterName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_setGlobalParameterName$handle() {
    return OpenMM_CustomExternalForce_setGlobalParameterName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_setGlobalParameterName$address() {
    return OpenMM_CustomExternalForce_setGlobalParameterName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterName(OpenMM_CustomExternalForce *target, int index, const char *name)
   *}
   */
  public static void OpenMM_CustomExternalForce_setGlobalParameterName(MemorySegment target, int index, MemorySegment name) {
    var mh$ = OpenMM_CustomExternalForce_setGlobalParameterName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_setGlobalParameterName", target, index, name);
      }
      mh$.invokeExact(target, index, name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getGlobalParameterDefaultValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getGlobalParameterDefaultValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getGlobalParameterDefaultValue$descriptor() {
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getGlobalParameterDefaultValue$handle() {
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getGlobalParameterDefaultValue$address() {
    return OpenMM_CustomExternalForce_getGlobalParameterDefaultValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(const OpenMM_CustomExternalForce *target, int index)
   *}
   */
  public static double OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(MemorySegment target, int index) {
    var mh$ = OpenMM_CustomExternalForce_getGlobalParameterDefaultValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getGlobalParameterDefaultValue", target, index);
      }
      return (double) mh$.invokeExact(target, index);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_setGlobalParameterDefaultValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_setGlobalParameterDefaultValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(OpenMM_CustomExternalForce *target, int index, double defaultValue)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_setGlobalParameterDefaultValue$descriptor() {
    return OpenMM_CustomExternalForce_setGlobalParameterDefaultValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(OpenMM_CustomExternalForce *target, int index, double defaultValue)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_setGlobalParameterDefaultValue$handle() {
    return OpenMM_CustomExternalForce_setGlobalParameterDefaultValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(OpenMM_CustomExternalForce *target, int index, double defaultValue)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_setGlobalParameterDefaultValue$address() {
    return OpenMM_CustomExternalForce_setGlobalParameterDefaultValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(OpenMM_CustomExternalForce *target, int index, double defaultValue)
   *}
   */
  public static void OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(MemorySegment target, int index, double defaultValue) {
    var mh$ = OpenMM_CustomExternalForce_setGlobalParameterDefaultValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_setGlobalParameterDefaultValue", target, index, defaultValue);
      }
      mh$.invokeExact(target, index, defaultValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addParticle(OpenMM_CustomExternalForce *target, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_addParticle$descriptor() {
    return OpenMM_CustomExternalForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addParticle(OpenMM_CustomExternalForce *target, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_addParticle$handle() {
    return OpenMM_CustomExternalForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addParticle(OpenMM_CustomExternalForce *target, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_addParticle$address() {
    return OpenMM_CustomExternalForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CustomExternalForce_addParticle(OpenMM_CustomExternalForce *target, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static int OpenMM_CustomExternalForce_addParticle(MemorySegment target, int particle, MemorySegment parameters) {
    var mh$ = OpenMM_CustomExternalForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_addParticle", target, particle, parameters);
      }
      return (int) mh$.invokeExact(target, particle, parameters);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_getParticleParameters(const OpenMM_CustomExternalForce *target, int index, int *particle, OpenMM_DoubleArray *parameters)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_getParticleParameters$descriptor() {
    return OpenMM_CustomExternalForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_getParticleParameters(const OpenMM_CustomExternalForce *target, int index, int *particle, OpenMM_DoubleArray *parameters)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_getParticleParameters$handle() {
    return OpenMM_CustomExternalForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_getParticleParameters(const OpenMM_CustomExternalForce *target, int index, int *particle, OpenMM_DoubleArray *parameters)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_getParticleParameters$address() {
    return OpenMM_CustomExternalForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_getParticleParameters(const OpenMM_CustomExternalForce *target, int index, int *particle, OpenMM_DoubleArray *parameters)
   *}
   */
  public static void OpenMM_CustomExternalForce_getParticleParameters(MemorySegment target, int index, MemorySegment particle, MemorySegment parameters) {
    var mh$ = OpenMM_CustomExternalForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_getParticleParameters", target, index, particle, parameters);
      }
      mh$.invokeExact(target, index, particle, parameters);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setParticleParameters(OpenMM_CustomExternalForce *target, int index, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_setParticleParameters$descriptor() {
    return OpenMM_CustomExternalForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setParticleParameters(OpenMM_CustomExternalForce *target, int index, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_setParticleParameters$handle() {
    return OpenMM_CustomExternalForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setParticleParameters(OpenMM_CustomExternalForce *target, int index, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_setParticleParameters$address() {
    return OpenMM_CustomExternalForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_setParticleParameters(OpenMM_CustomExternalForce *target, int index, int particle, const OpenMM_DoubleArray *parameters)
   *}
   */
  public static void OpenMM_CustomExternalForce_setParticleParameters(MemorySegment target, int index, int particle, MemorySegment parameters) {
    var mh$ = OpenMM_CustomExternalForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_setParticleParameters", target, index, particle, parameters);
      }
      mh$.invokeExact(target, index, particle, parameters);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_updateParametersInContext(OpenMM_CustomExternalForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_updateParametersInContext$descriptor() {
    return OpenMM_CustomExternalForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_updateParametersInContext(OpenMM_CustomExternalForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_updateParametersInContext$handle() {
    return OpenMM_CustomExternalForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_updateParametersInContext(OpenMM_CustomExternalForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_updateParametersInContext$address() {
    return OpenMM_CustomExternalForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CustomExternalForce_updateParametersInContext(OpenMM_CustomExternalForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_CustomExternalForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_CustomExternalForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MethodHandle OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static MemorySegment OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(const OpenMM_CustomExternalForce *target)
   *}
   */
  public static int OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Continuous2DFunction_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Continuous2DFunction_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax, OpenMM_Boolean periodic)
   *}
   */
  public static FunctionDescriptor OpenMM_Continuous2DFunction_create$descriptor() {
    return OpenMM_Continuous2DFunction_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax, OpenMM_Boolean periodic)
   *}
   */
  public static MethodHandle OpenMM_Continuous2DFunction_create$handle() {
    return OpenMM_Continuous2DFunction_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax, OpenMM_Boolean periodic)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_create$address() {
    return OpenMM_Continuous2DFunction_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_create(int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax, OpenMM_Boolean periodic)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_create(int xsize, int ysize, MemorySegment values, double xmin, double xmax, double ymin, double ymax, int periodic) {
    var mh$ = OpenMM_Continuous2DFunction_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Continuous2DFunction_create", xsize, ysize, values, xmin, xmax, ymin, ymax, periodic);
      }
      return (MemorySegment) mh$.invokeExact(xsize, ysize, values, xmin, xmax, ymin, ymax, periodic);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Continuous2DFunction_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Continuous2DFunction_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_destroy(OpenMM_Continuous2DFunction *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Continuous2DFunction_destroy$descriptor() {
    return OpenMM_Continuous2DFunction_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_destroy(OpenMM_Continuous2DFunction *target)
   *}
   */
  public static MethodHandle OpenMM_Continuous2DFunction_destroy$handle() {
    return OpenMM_Continuous2DFunction_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_destroy(OpenMM_Continuous2DFunction *target)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_destroy$address() {
    return OpenMM_Continuous2DFunction_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_destroy(OpenMM_Continuous2DFunction *target)
   *}
   */
  public static void OpenMM_Continuous2DFunction_destroy(MemorySegment target) {
    var mh$ = OpenMM_Continuous2DFunction_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Continuous2DFunction_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Continuous2DFunction_getFunctionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Continuous2DFunction_getFunctionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_getFunctionParameters(const OpenMM_Continuous2DFunction *target, int *xsize, int *ysize, OpenMM_DoubleArray *values, double *xmin, double *xmax, double *ymin, double *ymax)
   *}
   */
  public static FunctionDescriptor OpenMM_Continuous2DFunction_getFunctionParameters$descriptor() {
    return OpenMM_Continuous2DFunction_getFunctionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_getFunctionParameters(const OpenMM_Continuous2DFunction *target, int *xsize, int *ysize, OpenMM_DoubleArray *values, double *xmin, double *xmax, double *ymin, double *ymax)
   *}
   */
  public static MethodHandle OpenMM_Continuous2DFunction_getFunctionParameters$handle() {
    return OpenMM_Continuous2DFunction_getFunctionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_getFunctionParameters(const OpenMM_Continuous2DFunction *target, int *xsize, int *ysize, OpenMM_DoubleArray *values, double *xmin, double *xmax, double *ymin, double *ymax)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_getFunctionParameters$address() {
    return OpenMM_Continuous2DFunction_getFunctionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_getFunctionParameters(const OpenMM_Continuous2DFunction *target, int *xsize, int *ysize, OpenMM_DoubleArray *values, double *xmin, double *xmax, double *ymin, double *ymax)
   *}
   */
  public static void OpenMM_Continuous2DFunction_getFunctionParameters(MemorySegment target, MemorySegment xsize, MemorySegment ysize, MemorySegment values, MemorySegment xmin, MemorySegment xmax, MemorySegment ymin, MemorySegment ymax) {
    var mh$ = OpenMM_Continuous2DFunction_getFunctionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Continuous2DFunction_getFunctionParameters", target, xsize, ysize, values, xmin, xmax, ymin, ymax);
      }
      mh$.invokeExact(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Continuous2DFunction_setFunctionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Continuous2DFunction_setFunctionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_setFunctionParameters(OpenMM_Continuous2DFunction *target, int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax)
   *}
   */
  public static FunctionDescriptor OpenMM_Continuous2DFunction_setFunctionParameters$descriptor() {
    return OpenMM_Continuous2DFunction_setFunctionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_setFunctionParameters(OpenMM_Continuous2DFunction *target, int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax)
   *}
   */
  public static MethodHandle OpenMM_Continuous2DFunction_setFunctionParameters$handle() {
    return OpenMM_Continuous2DFunction_setFunctionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_setFunctionParameters(OpenMM_Continuous2DFunction *target, int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_setFunctionParameters$address() {
    return OpenMM_Continuous2DFunction_setFunctionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Continuous2DFunction_setFunctionParameters(OpenMM_Continuous2DFunction *target, int xsize, int ysize, const OpenMM_DoubleArray *values, double xmin, double xmax, double ymin, double ymax)
   *}
   */
  public static void OpenMM_Continuous2DFunction_setFunctionParameters(MemorySegment target, int xsize, int ysize, MemorySegment values, double xmin, double xmax, double ymin, double ymax) {
    var mh$ = OpenMM_Continuous2DFunction_setFunctionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Continuous2DFunction_setFunctionParameters", target, xsize, ysize, values, xmin, xmax, ymin, ymax);
      }
      mh$.invokeExact(target, xsize, ysize, values, xmin, xmax, ymin, ymax);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Continuous2DFunction_Copy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Continuous2DFunction_Copy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_Copy(const OpenMM_Continuous2DFunction *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Continuous2DFunction_Copy$descriptor() {
    return OpenMM_Continuous2DFunction_Copy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_Copy(const OpenMM_Continuous2DFunction *target)
   *}
   */
  public static MethodHandle OpenMM_Continuous2DFunction_Copy$handle() {
    return OpenMM_Continuous2DFunction_Copy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_Copy(const OpenMM_Continuous2DFunction *target)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_Copy$address() {
    return OpenMM_Continuous2DFunction_Copy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Continuous2DFunction *OpenMM_Continuous2DFunction_Copy(const OpenMM_Continuous2DFunction *target)
   *}
   */
  public static MemorySegment OpenMM_Continuous2DFunction_Copy(MemorySegment target) {
    var mh$ = OpenMM_Continuous2DFunction_Copy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Continuous2DFunction_Copy", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CMMotionRemover_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CMMotionRemover_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_CMMotionRemover *OpenMM_CMMotionRemover_create(int frequency)
   *}
   */
  public static FunctionDescriptor OpenMM_CMMotionRemover_create$descriptor() {
    return OpenMM_CMMotionRemover_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_CMMotionRemover *OpenMM_CMMotionRemover_create(int frequency)
   *}
   */
  public static MethodHandle OpenMM_CMMotionRemover_create$handle() {
    return OpenMM_CMMotionRemover_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_CMMotionRemover *OpenMM_CMMotionRemover_create(int frequency)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_create$address() {
    return OpenMM_CMMotionRemover_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_CMMotionRemover *OpenMM_CMMotionRemover_create(int frequency)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_create(int frequency) {
    var mh$ = OpenMM_CMMotionRemover_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CMMotionRemover_create", frequency);
      }
      return (MemorySegment) mh$.invokeExact(frequency);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CMMotionRemover_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CMMotionRemover_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CMMotionRemover_destroy$descriptor() {
    return OpenMM_CMMotionRemover_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover *target)
   *}
   */
  public static MethodHandle OpenMM_CMMotionRemover_destroy$handle() {
    return OpenMM_CMMotionRemover_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover *target)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_destroy$address() {
    return OpenMM_CMMotionRemover_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_destroy(OpenMM_CMMotionRemover *target)
   *}
   */
  public static void OpenMM_CMMotionRemover_destroy(MemorySegment target) {
    var mh$ = OpenMM_CMMotionRemover_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CMMotionRemover_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CMMotionRemover_getFrequency {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CMMotionRemover_getFrequency");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CMMotionRemover_getFrequency$descriptor() {
    return OpenMM_CMMotionRemover_getFrequency.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static MethodHandle OpenMM_CMMotionRemover_getFrequency$handle() {
    return OpenMM_CMMotionRemover_getFrequency.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_getFrequency$address() {
    return OpenMM_CMMotionRemover_getFrequency.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_CMMotionRemover_getFrequency(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static int OpenMM_CMMotionRemover_getFrequency(MemorySegment target) {
    var mh$ = OpenMM_CMMotionRemover_getFrequency.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CMMotionRemover_getFrequency", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CMMotionRemover_setFrequency {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CMMotionRemover_setFrequency");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover *target, int freq)
   *}
   */
  public static FunctionDescriptor OpenMM_CMMotionRemover_setFrequency$descriptor() {
    return OpenMM_CMMotionRemover_setFrequency.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover *target, int freq)
   *}
   */
  public static MethodHandle OpenMM_CMMotionRemover_setFrequency$handle() {
    return OpenMM_CMMotionRemover_setFrequency.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover *target, int freq)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_setFrequency$address() {
    return OpenMM_CMMotionRemover_setFrequency.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_CMMotionRemover_setFrequency(OpenMM_CMMotionRemover *target, int freq)
   *}
   */
  public static void OpenMM_CMMotionRemover_setFrequency(MemorySegment target, int freq) {
    var mh$ = OpenMM_CMMotionRemover_setFrequency.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CMMotionRemover_setFrequency", target, freq);
      }
      mh$.invokeExact(target, freq);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static FunctionDescriptor OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static MethodHandle OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static MemorySegment OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions$address() {
    return OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(const OpenMM_CMMotionRemover *target)
   *}
   */
  public static int OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_destroy(OpenMM_Platform *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_destroy$descriptor() {
    return OpenMM_Platform_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_destroy(OpenMM_Platform *target)
   *}
   */
  public static MethodHandle OpenMM_Platform_destroy$handle() {
    return OpenMM_Platform_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_destroy(OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_destroy$address() {
    return OpenMM_Platform_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Platform_destroy(OpenMM_Platform *target)
   *}
   */
  public static void OpenMM_Platform_destroy(MemorySegment target) {
    var mh$ = OpenMM_Platform_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_registerPlatform {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_registerPlatform");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_registerPlatform(OpenMM_Platform *platform)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_registerPlatform$descriptor() {
    return OpenMM_Platform_registerPlatform.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_registerPlatform(OpenMM_Platform *platform)
   *}
   */
  public static MethodHandle OpenMM_Platform_registerPlatform$handle() {
    return OpenMM_Platform_registerPlatform.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_registerPlatform(OpenMM_Platform *platform)
   *}
   */
  public static MemorySegment OpenMM_Platform_registerPlatform$address() {
    return OpenMM_Platform_registerPlatform.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Platform_registerPlatform(OpenMM_Platform *platform)
   *}
   */
  public static void OpenMM_Platform_registerPlatform(MemorySegment platform) {
    var mh$ = OpenMM_Platform_registerPlatform.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_registerPlatform", platform);
      }
      mh$.invokeExact(platform);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern int OpenMM_Platform_getNumPlatforms()
   *}
   */
  public static class OpenMM_Platform_getNumPlatforms {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getNumPlatforms");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_Platform_getNumPlatforms(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern int OpenMM_Platform_getNumPlatforms()
     *}
     */
    public static OpenMM_Platform_getNumPlatforms makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_Platform_getNumPlatforms(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public int apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_Platform_getNumPlatforms", x0);
        }
        return (int) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_Platform_getPlatform {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPlatform");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform(int index)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPlatform$descriptor() {
    return OpenMM_Platform_getPlatform.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform(int index)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPlatform$handle() {
    return OpenMM_Platform_getPlatform.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform(int index)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatform$address() {
    return OpenMM_Platform_getPlatform.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform(int index)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatform(int index) {
    var mh$ = OpenMM_Platform_getPlatform.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPlatform", index);
      }
      return (MemorySegment) mh$.invokeExact(index);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getPlatform_1 {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPlatform_1");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform_1(const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPlatform_1$descriptor() {
    return OpenMM_Platform_getPlatform_1.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform_1(const char *name)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPlatform_1$handle() {
    return OpenMM_Platform_getPlatform_1.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform_1(const char *name)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatform_1$address() {
    return OpenMM_Platform_getPlatform_1.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatform_1(const char *name)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatform_1(MemorySegment name) {
    var mh$ = OpenMM_Platform_getPlatform_1.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPlatform_1", name);
      }
      return (MemorySegment) mh$.invokeExact(name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getPlatformByName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPlatformByName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatformByName(const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPlatformByName$descriptor() {
    return OpenMM_Platform_getPlatformByName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatformByName(const char *name)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPlatformByName$handle() {
    return OpenMM_Platform_getPlatformByName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatformByName(const char *name)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatformByName$address() {
    return OpenMM_Platform_getPlatformByName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_getPlatformByName(const char *name)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPlatformByName(MemorySegment name) {
    var mh$ = OpenMM_Platform_getPlatformByName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPlatformByName", name);
      }
      return (MemorySegment) mh$.invokeExact(name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_findPlatform {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_findPlatform");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_findPlatform(const OpenMM_StringArray *kernelNames)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_findPlatform$descriptor() {
    return OpenMM_Platform_findPlatform.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_findPlatform(const OpenMM_StringArray *kernelNames)
   *}
   */
  public static MethodHandle OpenMM_Platform_findPlatform$handle() {
    return OpenMM_Platform_findPlatform.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_findPlatform(const OpenMM_StringArray *kernelNames)
   *}
   */
  public static MemorySegment OpenMM_Platform_findPlatform$address() {
    return OpenMM_Platform_findPlatform.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Platform *OpenMM_Platform_findPlatform(const OpenMM_StringArray *kernelNames)
   *}
   */
  public static MemorySegment OpenMM_Platform_findPlatform(MemorySegment kernelNames) {
    var mh$ = OpenMM_Platform_findPlatform.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_findPlatform", kernelNames);
      }
      return (MemorySegment) mh$.invokeExact(kernelNames);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_loadPluginLibrary {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_loadPluginLibrary");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_loadPluginLibrary(const char *file)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_loadPluginLibrary$descriptor() {
    return OpenMM_Platform_loadPluginLibrary.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_loadPluginLibrary(const char *file)
   *}
   */
  public static MethodHandle OpenMM_Platform_loadPluginLibrary$handle() {
    return OpenMM_Platform_loadPluginLibrary.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_loadPluginLibrary(const char *file)
   *}
   */
  public static MemorySegment OpenMM_Platform_loadPluginLibrary$address() {
    return OpenMM_Platform_loadPluginLibrary.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Platform_loadPluginLibrary(const char *file)
   *}
   */
  public static void OpenMM_Platform_loadPluginLibrary(MemorySegment file) {
    var mh$ = OpenMM_Platform_loadPluginLibrary.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_loadPluginLibrary", file);
      }
      mh$.invokeExact(file);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getDefaultPluginsDirectory()
   *}
   */
  public static class OpenMM_Platform_getDefaultPluginsDirectory {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getDefaultPluginsDirectory");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_Platform_getDefaultPluginsDirectory(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern const char *OpenMM_Platform_getDefaultPluginsDirectory()
     *}
     */
    public static OpenMM_Platform_getDefaultPluginsDirectory makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_Platform_getDefaultPluginsDirectory(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_Platform_getDefaultPluginsDirectory", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getOpenMMVersion()
   *}
   */
  public static class OpenMM_Platform_getOpenMMVersion {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getOpenMMVersion");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_Platform_getOpenMMVersion(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern const char *OpenMM_Platform_getOpenMMVersion()
     *}
     */
    public static OpenMM_Platform_getOpenMMVersion makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_Platform_getOpenMMVersion(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_Platform_getOpenMMVersion", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_Platform_getName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getName(const OpenMM_Platform *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getName$descriptor() {
    return OpenMM_Platform_getName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getName(const OpenMM_Platform *target)
   *}
   */
  public static MethodHandle OpenMM_Platform_getName$handle() {
    return OpenMM_Platform_getName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getName(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_getName$address() {
    return OpenMM_Platform_getName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getName(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_getName(MemorySegment target) {
    var mh$ = OpenMM_Platform_getName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getName", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getSpeed {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getSpeed");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_Platform_getSpeed(const OpenMM_Platform *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getSpeed$descriptor() {
    return OpenMM_Platform_getSpeed.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_Platform_getSpeed(const OpenMM_Platform *target)
   *}
   */
  public static MethodHandle OpenMM_Platform_getSpeed$handle() {
    return OpenMM_Platform_getSpeed.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_Platform_getSpeed(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_getSpeed$address() {
    return OpenMM_Platform_getSpeed.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_Platform_getSpeed(const OpenMM_Platform *target)
   *}
   */
  public static double OpenMM_Platform_getSpeed(MemorySegment target) {
    var mh$ = OpenMM_Platform_getSpeed.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getSpeed", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_supportsDoublePrecision {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_supportsDoublePrecision");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_supportsDoublePrecision$descriptor() {
    return OpenMM_Platform_supportsDoublePrecision.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform *target)
   *}
   */
  public static MethodHandle OpenMM_Platform_supportsDoublePrecision$handle() {
    return OpenMM_Platform_supportsDoublePrecision.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_supportsDoublePrecision$address() {
    return OpenMM_Platform_supportsDoublePrecision.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsDoublePrecision(const OpenMM_Platform *target)
   *}
   */
  public static int OpenMM_Platform_supportsDoublePrecision(MemorySegment target) {
    var mh$ = OpenMM_Platform_supportsDoublePrecision.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_supportsDoublePrecision", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getPropertyNames {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPropertyNames");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_StringArray *OpenMM_Platform_getPropertyNames(const OpenMM_Platform *target)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPropertyNames$descriptor() {
    return OpenMM_Platform_getPropertyNames.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_StringArray *OpenMM_Platform_getPropertyNames(const OpenMM_Platform *target)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPropertyNames$handle() {
    return OpenMM_Platform_getPropertyNames.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_StringArray *OpenMM_Platform_getPropertyNames(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyNames$address() {
    return OpenMM_Platform_getPropertyNames.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_StringArray *OpenMM_Platform_getPropertyNames(const OpenMM_Platform *target)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyNames(MemorySegment target) {
    var mh$ = OpenMM_Platform_getPropertyNames.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPropertyNames", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getPropertyValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPropertyValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyValue(const OpenMM_Platform *target, const OpenMM_Context *context, const char *property)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPropertyValue$descriptor() {
    return OpenMM_Platform_getPropertyValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyValue(const OpenMM_Platform *target, const OpenMM_Context *context, const char *property)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPropertyValue$handle() {
    return OpenMM_Platform_getPropertyValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyValue(const OpenMM_Platform *target, const OpenMM_Context *context, const char *property)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyValue$address() {
    return OpenMM_Platform_getPropertyValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyValue(const OpenMM_Platform *target, const OpenMM_Context *context, const char *property)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyValue(MemorySegment target, MemorySegment context, MemorySegment property) {
    var mh$ = OpenMM_Platform_getPropertyValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPropertyValue", target, context, property);
      }
      return (MemorySegment) mh$.invokeExact(target, context, property);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_setPropertyValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_setPropertyValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyValue(const OpenMM_Platform *target, OpenMM_Context *context, const char *property, const char *value)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_setPropertyValue$descriptor() {
    return OpenMM_Platform_setPropertyValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyValue(const OpenMM_Platform *target, OpenMM_Context *context, const char *property, const char *value)
   *}
   */
  public static MethodHandle OpenMM_Platform_setPropertyValue$handle() {
    return OpenMM_Platform_setPropertyValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyValue(const OpenMM_Platform *target, OpenMM_Context *context, const char *property, const char *value)
   *}
   */
  public static MemorySegment OpenMM_Platform_setPropertyValue$address() {
    return OpenMM_Platform_setPropertyValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyValue(const OpenMM_Platform *target, OpenMM_Context *context, const char *property, const char *value)
   *}
   */
  public static void OpenMM_Platform_setPropertyValue(MemorySegment target, MemorySegment context, MemorySegment property, MemorySegment value) {
    var mh$ = OpenMM_Platform_setPropertyValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_setPropertyValue", target, context, property, value);
      }
      mh$.invokeExact(target, context, property, value);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_getPropertyDefaultValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_getPropertyDefaultValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform *target, const char *property)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_getPropertyDefaultValue$descriptor() {
    return OpenMM_Platform_getPropertyDefaultValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform *target, const char *property)
   *}
   */
  public static MethodHandle OpenMM_Platform_getPropertyDefaultValue$handle() {
    return OpenMM_Platform_getPropertyDefaultValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform *target, const char *property)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyDefaultValue$address() {
    return OpenMM_Platform_getPropertyDefaultValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_Platform_getPropertyDefaultValue(const OpenMM_Platform *target, const char *property)
   *}
   */
  public static MemorySegment OpenMM_Platform_getPropertyDefaultValue(MemorySegment target, MemorySegment property) {
    var mh$ = OpenMM_Platform_getPropertyDefaultValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_getPropertyDefaultValue", target, property);
      }
      return (MemorySegment) mh$.invokeExact(target, property);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_setPropertyDefaultValue {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_setPropertyDefaultValue");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform *target, const char *property, const char *value)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_setPropertyDefaultValue$descriptor() {
    return OpenMM_Platform_setPropertyDefaultValue.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform *target, const char *property, const char *value)
   *}
   */
  public static MethodHandle OpenMM_Platform_setPropertyDefaultValue$handle() {
    return OpenMM_Platform_setPropertyDefaultValue.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform *target, const char *property, const char *value)
   *}
   */
  public static MemorySegment OpenMM_Platform_setPropertyDefaultValue$address() {
    return OpenMM_Platform_setPropertyDefaultValue.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_Platform_setPropertyDefaultValue(OpenMM_Platform *target, const char *property, const char *value)
   *}
   */
  public static void OpenMM_Platform_setPropertyDefaultValue(MemorySegment target, MemorySegment property, MemorySegment value) {
    var mh$ = OpenMM_Platform_setPropertyDefaultValue.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_setPropertyDefaultValue", target, property, value);
      }
      mh$.invokeExact(target, property, value);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_Platform_supportsKernels {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_Platform_supportsKernels");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform *target, const OpenMM_StringArray *kernelNames)
   *}
   */
  public static FunctionDescriptor OpenMM_Platform_supportsKernels$descriptor() {
    return OpenMM_Platform_supportsKernels.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform *target, const OpenMM_StringArray *kernelNames)
   *}
   */
  public static MethodHandle OpenMM_Platform_supportsKernels$handle() {
    return OpenMM_Platform_supportsKernels.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform *target, const OpenMM_StringArray *kernelNames)
   *}
   */
  public static MemorySegment OpenMM_Platform_supportsKernels$address() {
    return OpenMM_Platform_supportsKernels.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_Platform_supportsKernels(const OpenMM_Platform *target, const OpenMM_StringArray *kernelNames)
   *}
   */
  public static int OpenMM_Platform_supportsKernels(MemorySegment target, MemorySegment kernelNames) {
    var mh$ = OpenMM_Platform_supportsKernels.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_Platform_supportsKernels", target, kernelNames);
      }
      return (int) mh$.invokeExact(target, kernelNames);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_3D_DoubleArray_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_3D_DoubleArray_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * OpenMM_3D_DoubleArray *OpenMM_3D_DoubleArray_create(int size1, int size2, int size3)
   *}
   */
  public static FunctionDescriptor OpenMM_3D_DoubleArray_create$descriptor() {
    return OpenMM_3D_DoubleArray_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * OpenMM_3D_DoubleArray *OpenMM_3D_DoubleArray_create(int size1, int size2, int size3)
   *}
   */
  public static MethodHandle OpenMM_3D_DoubleArray_create$handle() {
    return OpenMM_3D_DoubleArray_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * OpenMM_3D_DoubleArray *OpenMM_3D_DoubleArray_create(int size1, int size2, int size3)
   *}
   */
  public static MemorySegment OpenMM_3D_DoubleArray_create$address() {
    return OpenMM_3D_DoubleArray_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * OpenMM_3D_DoubleArray *OpenMM_3D_DoubleArray_create(int size1, int size2, int size3)
   *}
   */
  public static MemorySegment OpenMM_3D_DoubleArray_create(int size1, int size2, int size3) {
    var mh$ = OpenMM_3D_DoubleArray_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_3D_DoubleArray_create", size1, size2, size3);
      }
      return (MemorySegment) mh$.invokeExact(size1, size2, size3);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_3D_DoubleArray_set {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_3D_DoubleArray_set");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray *array, int index1, int index2, OpenMM_DoubleArray *values)
   *}
   */
  public static FunctionDescriptor OpenMM_3D_DoubleArray_set$descriptor() {
    return OpenMM_3D_DoubleArray_set.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray *array, int index1, int index2, OpenMM_DoubleArray *values)
   *}
   */
  public static MethodHandle OpenMM_3D_DoubleArray_set$handle() {
    return OpenMM_3D_DoubleArray_set.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray *array, int index1, int index2, OpenMM_DoubleArray *values)
   *}
   */
  public static MemorySegment OpenMM_3D_DoubleArray_set$address() {
    return OpenMM_3D_DoubleArray_set.ADDR;
  }

  /**
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray *array, int index1, int index2, OpenMM_DoubleArray *values)
   *}
   */
  public static void OpenMM_3D_DoubleArray_set(MemorySegment array, int index1, int index2, MemorySegment values) {
    var mh$ = OpenMM_3D_DoubleArray_set.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_3D_DoubleArray_set", array, index1, index2, values);
      }
      mh$.invokeExact(array, index1, index2, values);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_3D_DoubleArray_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_3D_DoubleArray_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray *array)
   *}
   */
  public static FunctionDescriptor OpenMM_3D_DoubleArray_destroy$descriptor() {
    return OpenMM_3D_DoubleArray_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray *array)
   *}
   */
  public static MethodHandle OpenMM_3D_DoubleArray_destroy$handle() {
    return OpenMM_3D_DoubleArray_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray *array)
   *}
   */
  public static MemorySegment OpenMM_3D_DoubleArray_destroy$address() {
    return OpenMM_3D_DoubleArray_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray *array)
   *}
   */
  public static void OpenMM_3D_DoubleArray_destroy(MemorySegment array) {
    var mh$ = OpenMM_3D_DoubleArray_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_3D_DoubleArray_destroy", array);
      }
      mh$.invokeExact(array);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static final int OpenMM_AmoebaVdwForce_NoCutoff = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_NoCutoff = 0
   *}
   */
  public static int OpenMM_AmoebaVdwForce_NoCutoff() {
    return OpenMM_AmoebaVdwForce_NoCutoff;
  }

  private static final int OpenMM_AmoebaVdwForce_CutoffPeriodic = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_CutoffPeriodic = 1
   *}
   */
  public static int OpenMM_AmoebaVdwForce_CutoffPeriodic() {
    return OpenMM_AmoebaVdwForce_CutoffPeriodic;
  }

  private static final int OpenMM_AmoebaVdwForce_Buffered147 = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_Buffered147 = 0
   *}
   */
  public static int OpenMM_AmoebaVdwForce_Buffered147() {
    return OpenMM_AmoebaVdwForce_Buffered147;
  }

  private static final int OpenMM_AmoebaVdwForce_LennardJones = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_LennardJones = 1
   *}
   */
  public static int OpenMM_AmoebaVdwForce_LennardJones() {
    return OpenMM_AmoebaVdwForce_LennardJones;
  }

  private static final int OpenMM_AmoebaVdwForce_None = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_None = 0
   *}
   */
  public static int OpenMM_AmoebaVdwForce_None() {
    return OpenMM_AmoebaVdwForce_None;
  }

  private static final int OpenMM_AmoebaVdwForce_Decouple = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_Decouple = 1
   *}
   */
  public static int OpenMM_AmoebaVdwForce_Decouple() {
    return OpenMM_AmoebaVdwForce_Decouple;
  }

  private static final int OpenMM_AmoebaVdwForce_Annihilate = (int) 2L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaVdwForce_Annihilate = 2
   *}
   */
  public static int OpenMM_AmoebaVdwForce_Annihilate() {
    return OpenMM_AmoebaVdwForce_Annihilate;
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce *OpenMM_AmoebaVdwForce_create()
   *}
   */
  public static class OpenMM_AmoebaVdwForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_AmoebaVdwForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_AmoebaVdwForce *OpenMM_AmoebaVdwForce_create()
     *}
     */
    public static OpenMM_AmoebaVdwForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_AmoebaVdwForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_AmoebaVdwForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_AmoebaVdwForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_destroy$descriptor() {
    return OpenMM_AmoebaVdwForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_destroy$handle() {
    return OpenMM_AmoebaVdwForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_destroy$address() {
    return OpenMM_AmoebaVdwForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_destroy(OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_Lambda {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_Lambda");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_Lambda(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_Lambda$descriptor() {
    return OpenMM_AmoebaVdwForce_Lambda.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_Lambda(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_Lambda$handle() {
    return OpenMM_AmoebaVdwForce_Lambda.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_Lambda(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_Lambda$address() {
    return OpenMM_AmoebaVdwForce_Lambda.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_Lambda(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_Lambda(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_Lambda.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_Lambda", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setLambdaName {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setLambdaName");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setLambdaName(OpenMM_AmoebaVdwForce *target, const char *name)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setLambdaName$descriptor() {
    return OpenMM_AmoebaVdwForce_setLambdaName.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setLambdaName(OpenMM_AmoebaVdwForce *target, const char *name)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setLambdaName$handle() {
    return OpenMM_AmoebaVdwForce_setLambdaName.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setLambdaName(OpenMM_AmoebaVdwForce *target, const char *name)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setLambdaName$address() {
    return OpenMM_AmoebaVdwForce_setLambdaName.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setLambdaName(OpenMM_AmoebaVdwForce *target, const char *name)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setLambdaName(MemorySegment target, MemorySegment name) {
    var mh$ = OpenMM_AmoebaVdwForce_setLambdaName.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setLambdaName", target, name);
      }
      mh$.invokeExact(target, name);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getNumParticles$descriptor() {
    return OpenMM_AmoebaVdwForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getNumParticles$handle() {
    return OpenMM_AmoebaVdwForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getNumParticles$address() {
    return OpenMM_AmoebaVdwForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticles(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getNumParticleTypes {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getNumParticleTypes");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getNumParticleTypes$descriptor() {
    return OpenMM_AmoebaVdwForce_getNumParticleTypes.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getNumParticleTypes$handle() {
    return OpenMM_AmoebaVdwForce_getNumParticleTypes.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getNumParticleTypes$address() {
    return OpenMM_AmoebaVdwForce_getNumParticleTypes.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getNumParticleTypes(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getNumParticleTypes.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getNumParticleTypes", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getNumTypePairs {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getNumTypePairs");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumTypePairs(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getNumTypePairs$descriptor() {
    return OpenMM_AmoebaVdwForce_getNumTypePairs.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumTypePairs(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getNumTypePairs$handle() {
    return OpenMM_AmoebaVdwForce_getNumTypePairs.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumTypePairs(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getNumTypePairs$address() {
    return OpenMM_AmoebaVdwForce_getNumTypePairs.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getNumTypePairs(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getNumTypePairs(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getNumTypePairs.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getNumTypePairs", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce *target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, int typeIndex, double scaleFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setParticleParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce *target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, int typeIndex, double scaleFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setParticleParameters$handle() {
    return OpenMM_AmoebaVdwForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce *target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, int typeIndex, double scaleFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setParticleParameters$address() {
    return OpenMM_AmoebaVdwForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleParameters(OpenMM_AmoebaVdwForce *target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, int typeIndex, double scaleFactor)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setParticleParameters(MemorySegment target, int particleIndex, int parentIndex, double sigma, double epsilon, double reductionFactor, int isAlchemical, int typeIndex, double scaleFactor) {
    var mh$ = OpenMM_AmoebaVdwForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setParticleParameters", target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, typeIndex, scaleFactor);
      }
      mh$.invokeExact(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, typeIndex, scaleFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce *target, int particleIndex, int *parentIndex, double *sigma, double *epsilon, double *reductionFactor, OpenMM_Boolean *isAlchemical, int *typeIndex, double *scaleFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getParticleParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce *target, int particleIndex, int *parentIndex, double *sigma, double *epsilon, double *reductionFactor, OpenMM_Boolean *isAlchemical, int *typeIndex, double *scaleFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getParticleParameters$handle() {
    return OpenMM_AmoebaVdwForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce *target, int particleIndex, int *parentIndex, double *sigma, double *epsilon, double *reductionFactor, OpenMM_Boolean *isAlchemical, int *typeIndex, double *scaleFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getParticleParameters$address() {
    return OpenMM_AmoebaVdwForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleParameters(const OpenMM_AmoebaVdwForce *target, int particleIndex, int *parentIndex, double *sigma, double *epsilon, double *reductionFactor, OpenMM_Boolean *isAlchemical, int *typeIndex, double *scaleFactor)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_getParticleParameters(MemorySegment target, int particleIndex, MemorySegment parentIndex, MemorySegment sigma, MemorySegment epsilon, MemorySegment reductionFactor, MemorySegment isAlchemical, MemorySegment typeIndex, MemorySegment scaleFactor) {
    var mh$ = OpenMM_AmoebaVdwForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getParticleParameters", target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, typeIndex, scaleFactor);
      }
      mh$.invokeExact(target, particleIndex, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, typeIndex, scaleFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce *target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_addParticle$descriptor() {
    return OpenMM_AmoebaVdwForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce *target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_addParticle$handle() {
    return OpenMM_AmoebaVdwForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce *target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_addParticle$address() {
    return OpenMM_AmoebaVdwForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle(OpenMM_AmoebaVdwForce *target, int parentIndex, double sigma, double epsilon, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_addParticle(MemorySegment target, int parentIndex, double sigma, double epsilon, double reductionFactor, int isAlchemical, double scaleFactor) {
    var mh$ = OpenMM_AmoebaVdwForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_addParticle", target, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, scaleFactor);
      }
      return (int) mh$.invokeExact(target, parentIndex, sigma, epsilon, reductionFactor, isAlchemical, scaleFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_addParticle_1 {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_addParticle_1");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle_1(OpenMM_AmoebaVdwForce *target, int parentIndex, int typeIndex, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_addParticle_1$descriptor() {
    return OpenMM_AmoebaVdwForce_addParticle_1.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle_1(OpenMM_AmoebaVdwForce *target, int parentIndex, int typeIndex, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_addParticle_1$handle() {
    return OpenMM_AmoebaVdwForce_addParticle_1.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle_1(OpenMM_AmoebaVdwForce *target, int parentIndex, int typeIndex, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_addParticle_1$address() {
    return OpenMM_AmoebaVdwForce_addParticle_1.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticle_1(OpenMM_AmoebaVdwForce *target, int parentIndex, int typeIndex, double reductionFactor, OpenMM_Boolean isAlchemical, double scaleFactor)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_addParticle_1(MemorySegment target, int parentIndex, int typeIndex, double reductionFactor, int isAlchemical, double scaleFactor) {
    var mh$ = OpenMM_AmoebaVdwForce_addParticle_1.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_addParticle_1", target, parentIndex, typeIndex, reductionFactor, isAlchemical, scaleFactor);
      }
      return (int) mh$.invokeExact(target, parentIndex, typeIndex, reductionFactor, isAlchemical, scaleFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_addParticleType {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_addParticleType");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticleType(OpenMM_AmoebaVdwForce *target, double sigma, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_addParticleType$descriptor() {
    return OpenMM_AmoebaVdwForce_addParticleType.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticleType(OpenMM_AmoebaVdwForce *target, double sigma, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_addParticleType$handle() {
    return OpenMM_AmoebaVdwForce_addParticleType.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticleType(OpenMM_AmoebaVdwForce *target, double sigma, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_addParticleType$address() {
    return OpenMM_AmoebaVdwForce_addParticleType.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addParticleType(OpenMM_AmoebaVdwForce *target, double sigma, double epsilon)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_addParticleType(MemorySegment target, double sigma, double epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_addParticleType.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_addParticleType", target, sigma, epsilon);
      }
      return (int) mh$.invokeExact(target, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getParticleTypeParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getParticleTypeParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleTypeParameters(const OpenMM_AmoebaVdwForce *target, int typeIndex, double *sigma, double *epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getParticleTypeParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_getParticleTypeParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleTypeParameters(const OpenMM_AmoebaVdwForce *target, int typeIndex, double *sigma, double *epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getParticleTypeParameters$handle() {
    return OpenMM_AmoebaVdwForce_getParticleTypeParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleTypeParameters(const OpenMM_AmoebaVdwForce *target, int typeIndex, double *sigma, double *epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getParticleTypeParameters$address() {
    return OpenMM_AmoebaVdwForce_getParticleTypeParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleTypeParameters(const OpenMM_AmoebaVdwForce *target, int typeIndex, double *sigma, double *epsilon)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_getParticleTypeParameters(MemorySegment target, int typeIndex, MemorySegment sigma, MemorySegment epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_getParticleTypeParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getParticleTypeParameters", target, typeIndex, sigma, epsilon);
      }
      mh$.invokeExact(target, typeIndex, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setParticleTypeParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setParticleTypeParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleTypeParameters(OpenMM_AmoebaVdwForce *target, int typeIndex, double sigma, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setParticleTypeParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_setParticleTypeParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleTypeParameters(OpenMM_AmoebaVdwForce *target, int typeIndex, double sigma, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setParticleTypeParameters$handle() {
    return OpenMM_AmoebaVdwForce_setParticleTypeParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleTypeParameters(OpenMM_AmoebaVdwForce *target, int typeIndex, double sigma, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setParticleTypeParameters$address() {
    return OpenMM_AmoebaVdwForce_setParticleTypeParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleTypeParameters(OpenMM_AmoebaVdwForce *target, int typeIndex, double sigma, double epsilon)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setParticleTypeParameters(MemorySegment target, int typeIndex, double sigma, double epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_setParticleTypeParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setParticleTypeParameters", target, typeIndex, sigma, epsilon);
      }
      mh$.invokeExact(target, typeIndex, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_addTypePair {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_addTypePair");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addTypePair(OpenMM_AmoebaVdwForce *target, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_addTypePair$descriptor() {
    return OpenMM_AmoebaVdwForce_addTypePair.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addTypePair(OpenMM_AmoebaVdwForce *target, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_addTypePair$handle() {
    return OpenMM_AmoebaVdwForce_addTypePair.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addTypePair(OpenMM_AmoebaVdwForce *target, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_addTypePair$address() {
    return OpenMM_AmoebaVdwForce_addTypePair.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_addTypePair(OpenMM_AmoebaVdwForce *target, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_addTypePair(MemorySegment target, int type1, int type2, double sigma, double epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_addTypePair.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_addTypePair", target, type1, type2, sigma, epsilon);
      }
      return (int) mh$.invokeExact(target, type1, type2, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getTypePairParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getTypePairParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getTypePairParameters(const OpenMM_AmoebaVdwForce *target, int pairIndex, int *type1, int *type2, double *sigma, double *epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getTypePairParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_getTypePairParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getTypePairParameters(const OpenMM_AmoebaVdwForce *target, int pairIndex, int *type1, int *type2, double *sigma, double *epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getTypePairParameters$handle() {
    return OpenMM_AmoebaVdwForce_getTypePairParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getTypePairParameters(const OpenMM_AmoebaVdwForce *target, int pairIndex, int *type1, int *type2, double *sigma, double *epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getTypePairParameters$address() {
    return OpenMM_AmoebaVdwForce_getTypePairParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getTypePairParameters(const OpenMM_AmoebaVdwForce *target, int pairIndex, int *type1, int *type2, double *sigma, double *epsilon)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_getTypePairParameters(MemorySegment target, int pairIndex, MemorySegment type1, MemorySegment type2, MemorySegment sigma, MemorySegment epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_getTypePairParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getTypePairParameters", target, pairIndex, type1, type2, sigma, epsilon);
      }
      mh$.invokeExact(target, pairIndex, type1, type2, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setTypePairParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setTypePairParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setTypePairParameters(OpenMM_AmoebaVdwForce *target, int pairIndex, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setTypePairParameters$descriptor() {
    return OpenMM_AmoebaVdwForce_setTypePairParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setTypePairParameters(OpenMM_AmoebaVdwForce *target, int pairIndex, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setTypePairParameters$handle() {
    return OpenMM_AmoebaVdwForce_setTypePairParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setTypePairParameters(OpenMM_AmoebaVdwForce *target, int pairIndex, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setTypePairParameters$address() {
    return OpenMM_AmoebaVdwForce_setTypePairParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setTypePairParameters(OpenMM_AmoebaVdwForce *target, int pairIndex, int type1, int type2, double sigma, double epsilon)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setTypePairParameters(MemorySegment target, int pairIndex, int type1, int type2, double sigma, double epsilon) {
    var mh$ = OpenMM_AmoebaVdwForce_setTypePairParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setTypePairParameters", target, pairIndex, type1, type2, sigma, epsilon);
      }
      mh$.invokeExact(target, pairIndex, type1, type2, sigma, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setSigmaCombiningRule {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setSigmaCombiningRule");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce *target, const char *sigmaCombiningRule)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setSigmaCombiningRule$descriptor() {
    return OpenMM_AmoebaVdwForce_setSigmaCombiningRule.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce *target, const char *sigmaCombiningRule)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setSigmaCombiningRule$handle() {
    return OpenMM_AmoebaVdwForce_setSigmaCombiningRule.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce *target, const char *sigmaCombiningRule)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setSigmaCombiningRule$address() {
    return OpenMM_AmoebaVdwForce_setSigmaCombiningRule.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(OpenMM_AmoebaVdwForce *target, const char *sigmaCombiningRule)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setSigmaCombiningRule(MemorySegment target, MemorySegment sigmaCombiningRule) {
    var mh$ = OpenMM_AmoebaVdwForce_setSigmaCombiningRule.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setSigmaCombiningRule", target, sigmaCombiningRule);
      }
      mh$.invokeExact(target, sigmaCombiningRule);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getSigmaCombiningRule {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getSigmaCombiningRule");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getSigmaCombiningRule$descriptor() {
    return OpenMM_AmoebaVdwForce_getSigmaCombiningRule.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getSigmaCombiningRule$handle() {
    return OpenMM_AmoebaVdwForce_getSigmaCombiningRule.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getSigmaCombiningRule$address() {
    return OpenMM_AmoebaVdwForce_getSigmaCombiningRule.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getSigmaCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getSigmaCombiningRule(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getSigmaCombiningRule.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getSigmaCombiningRule", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setEpsilonCombiningRule {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setEpsilonCombiningRule");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce *target, const char *epsilonCombiningRule)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setEpsilonCombiningRule$descriptor() {
    return OpenMM_AmoebaVdwForce_setEpsilonCombiningRule.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce *target, const char *epsilonCombiningRule)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setEpsilonCombiningRule$handle() {
    return OpenMM_AmoebaVdwForce_setEpsilonCombiningRule.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce *target, const char *epsilonCombiningRule)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setEpsilonCombiningRule$address() {
    return OpenMM_AmoebaVdwForce_setEpsilonCombiningRule.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(OpenMM_AmoebaVdwForce *target, const char *epsilonCombiningRule)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(MemorySegment target, MemorySegment epsilonCombiningRule) {
    var mh$ = OpenMM_AmoebaVdwForce_setEpsilonCombiningRule.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setEpsilonCombiningRule", target, epsilonCombiningRule);
      }
      mh$.invokeExact(target, epsilonCombiningRule);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getEpsilonCombiningRule {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getEpsilonCombiningRule");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getEpsilonCombiningRule$descriptor() {
    return OpenMM_AmoebaVdwForce_getEpsilonCombiningRule.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getEpsilonCombiningRule$handle() {
    return OpenMM_AmoebaVdwForce_getEpsilonCombiningRule.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getEpsilonCombiningRule$address() {
    return OpenMM_AmoebaVdwForce_getEpsilonCombiningRule.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const char *OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getEpsilonCombiningRule.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getEpsilonCombiningRule", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getUseDispersionCorrection {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getUseDispersionCorrection");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getUseDispersionCorrection$descriptor() {
    return OpenMM_AmoebaVdwForce_getUseDispersionCorrection.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getUseDispersionCorrection$handle() {
    return OpenMM_AmoebaVdwForce_getUseDispersionCorrection.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getUseDispersionCorrection$address() {
    return OpenMM_AmoebaVdwForce_getUseDispersionCorrection.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseDispersionCorrection(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getUseDispersionCorrection(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getUseDispersionCorrection.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getUseDispersionCorrection", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setUseDispersionCorrection {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setUseDispersionCorrection");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce *target, OpenMM_Boolean useCorrection)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setUseDispersionCorrection$descriptor() {
    return OpenMM_AmoebaVdwForce_setUseDispersionCorrection.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce *target, OpenMM_Boolean useCorrection)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setUseDispersionCorrection$handle() {
    return OpenMM_AmoebaVdwForce_setUseDispersionCorrection.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce *target, OpenMM_Boolean useCorrection)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setUseDispersionCorrection$address() {
    return OpenMM_AmoebaVdwForce_setUseDispersionCorrection.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(OpenMM_AmoebaVdwForce *target, OpenMM_Boolean useCorrection)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setUseDispersionCorrection(MemorySegment target, int useCorrection) {
    var mh$ = OpenMM_AmoebaVdwForce_setUseDispersionCorrection.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setUseDispersionCorrection", target, useCorrection);
      }
      mh$.invokeExact(target, useCorrection);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getUseParticleTypes {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getUseParticleTypes");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getUseParticleTypes$descriptor() {
    return OpenMM_AmoebaVdwForce_getUseParticleTypes.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getUseParticleTypes$handle() {
    return OpenMM_AmoebaVdwForce_getUseParticleTypes.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getUseParticleTypes$address() {
    return OpenMM_AmoebaVdwForce_getUseParticleTypes.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_getUseParticleTypes(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getUseParticleTypes(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getUseParticleTypes.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getUseParticleTypes", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setParticleExclusions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setParticleExclusions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce *target, int particleIndex, const OpenMM_IntArray *exclusions)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setParticleExclusions$descriptor() {
    return OpenMM_AmoebaVdwForce_setParticleExclusions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce *target, int particleIndex, const OpenMM_IntArray *exclusions)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setParticleExclusions$handle() {
    return OpenMM_AmoebaVdwForce_setParticleExclusions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce *target, int particleIndex, const OpenMM_IntArray *exclusions)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setParticleExclusions$address() {
    return OpenMM_AmoebaVdwForce_setParticleExclusions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setParticleExclusions(OpenMM_AmoebaVdwForce *target, int particleIndex, const OpenMM_IntArray *exclusions)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setParticleExclusions(MemorySegment target, int particleIndex, MemorySegment exclusions) {
    var mh$ = OpenMM_AmoebaVdwForce_setParticleExclusions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setParticleExclusions", target, particleIndex, exclusions);
      }
      mh$.invokeExact(target, particleIndex, exclusions);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getParticleExclusions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getParticleExclusions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce *target, int particleIndex, OpenMM_IntArray *exclusions)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getParticleExclusions$descriptor() {
    return OpenMM_AmoebaVdwForce_getParticleExclusions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce *target, int particleIndex, OpenMM_IntArray *exclusions)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getParticleExclusions$handle() {
    return OpenMM_AmoebaVdwForce_getParticleExclusions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce *target, int particleIndex, OpenMM_IntArray *exclusions)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getParticleExclusions$address() {
    return OpenMM_AmoebaVdwForce_getParticleExclusions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_getParticleExclusions(const OpenMM_AmoebaVdwForce *target, int particleIndex, OpenMM_IntArray *exclusions)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_getParticleExclusions(MemorySegment target, int particleIndex, MemorySegment exclusions) {
    var mh$ = OpenMM_AmoebaVdwForce_getParticleExclusions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getParticleExclusions", target, particleIndex, exclusions);
      }
      mh$.invokeExact(target, particleIndex, exclusions);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getCutoffDistance$descriptor() {
    return OpenMM_AmoebaVdwForce_getCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getCutoffDistance$handle() {
    return OpenMM_AmoebaVdwForce_getCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getCutoffDistance$address() {
    return OpenMM_AmoebaVdwForce_getCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoffDistance(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static double OpenMM_AmoebaVdwForce_getCutoffDistance(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getCutoffDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setCutoffDistance$descriptor() {
    return OpenMM_AmoebaVdwForce_setCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setCutoffDistance$handle() {
    return OpenMM_AmoebaVdwForce_setCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setCutoffDistance$address() {
    return OpenMM_AmoebaVdwForce_setCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoffDistance(OpenMM_AmoebaVdwForce *target, double distance)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setCutoffDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_AmoebaVdwForce_setCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setCutoffDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setCutoff {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setCutoff");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce *target, double cutoff)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setCutoff$descriptor() {
    return OpenMM_AmoebaVdwForce_setCutoff.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce *target, double cutoff)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setCutoff$handle() {
    return OpenMM_AmoebaVdwForce_setCutoff.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce *target, double cutoff)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setCutoff$address() {
    return OpenMM_AmoebaVdwForce_setCutoff.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setCutoff(OpenMM_AmoebaVdwForce *target, double cutoff)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setCutoff(MemorySegment target, double cutoff) {
    var mh$ = OpenMM_AmoebaVdwForce_setCutoff.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setCutoff", target, cutoff);
      }
      mh$.invokeExact(target, cutoff);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getCutoff {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getCutoff");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getCutoff$descriptor() {
    return OpenMM_AmoebaVdwForce_getCutoff.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getCutoff$handle() {
    return OpenMM_AmoebaVdwForce_getCutoff.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getCutoff$address() {
    return OpenMM_AmoebaVdwForce_getCutoff.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getCutoff(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static double OpenMM_AmoebaVdwForce_getCutoff(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getCutoff.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getCutoff", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getNonbondedMethod$descriptor() {
    return OpenMM_AmoebaVdwForce_getNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getNonbondedMethod$handle() {
    return OpenMM_AmoebaVdwForce_getNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getNonbondedMethod$address() {
    return OpenMM_AmoebaVdwForce_getNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_NonbondedMethod OpenMM_AmoebaVdwForce_getNonbondedMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getNonbondedMethod(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getNonbondedMethod", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_NonbondedMethod method)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setNonbondedMethod$descriptor() {
    return OpenMM_AmoebaVdwForce_setNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_NonbondedMethod method)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setNonbondedMethod$handle() {
    return OpenMM_AmoebaVdwForce_setNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_NonbondedMethod method)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setNonbondedMethod$address() {
    return OpenMM_AmoebaVdwForce_setNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setNonbondedMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_NonbondedMethod method)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setNonbondedMethod(MemorySegment target, int method) {
    var mh$ = OpenMM_AmoebaVdwForce_setNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setNonbondedMethod", target, method);
      }
      mh$.invokeExact(target, method);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getPotentialFunction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getPotentialFunction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_PotentialFunction OpenMM_AmoebaVdwForce_getPotentialFunction(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getPotentialFunction$descriptor() {
    return OpenMM_AmoebaVdwForce_getPotentialFunction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_PotentialFunction OpenMM_AmoebaVdwForce_getPotentialFunction(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getPotentialFunction$handle() {
    return OpenMM_AmoebaVdwForce_getPotentialFunction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_PotentialFunction OpenMM_AmoebaVdwForce_getPotentialFunction(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getPotentialFunction$address() {
    return OpenMM_AmoebaVdwForce_getPotentialFunction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_PotentialFunction OpenMM_AmoebaVdwForce_getPotentialFunction(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getPotentialFunction(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getPotentialFunction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getPotentialFunction", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setPotentialFunction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setPotentialFunction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setPotentialFunction(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_PotentialFunction potential)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setPotentialFunction$descriptor() {
    return OpenMM_AmoebaVdwForce_setPotentialFunction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setPotentialFunction(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_PotentialFunction potential)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setPotentialFunction$handle() {
    return OpenMM_AmoebaVdwForce_setPotentialFunction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setPotentialFunction(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_PotentialFunction potential)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setPotentialFunction$address() {
    return OpenMM_AmoebaVdwForce_setPotentialFunction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setPotentialFunction(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_PotentialFunction potential)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setPotentialFunction(MemorySegment target, int potential) {
    var mh$ = OpenMM_AmoebaVdwForce_setPotentialFunction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setPotentialFunction", target, potential);
      }
      mh$.invokeExact(target, potential);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setSoftcorePower {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setSoftcorePower");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce *target, int n)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setSoftcorePower$descriptor() {
    return OpenMM_AmoebaVdwForce_setSoftcorePower.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce *target, int n)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setSoftcorePower$handle() {
    return OpenMM_AmoebaVdwForce_setSoftcorePower.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce *target, int n)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setSoftcorePower$address() {
    return OpenMM_AmoebaVdwForce_setSoftcorePower.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcorePower(OpenMM_AmoebaVdwForce *target, int n)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setSoftcorePower(MemorySegment target, int n) {
    var mh$ = OpenMM_AmoebaVdwForce_setSoftcorePower.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setSoftcorePower", target, n);
      }
      mh$.invokeExact(target, n);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getSoftcorePower {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getSoftcorePower");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getSoftcorePower$descriptor() {
    return OpenMM_AmoebaVdwForce_getSoftcorePower.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getSoftcorePower$handle() {
    return OpenMM_AmoebaVdwForce_getSoftcorePower.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getSoftcorePower$address() {
    return OpenMM_AmoebaVdwForce_getSoftcorePower.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaVdwForce_getSoftcorePower(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getSoftcorePower(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getSoftcorePower.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getSoftcorePower", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setSoftcoreAlpha {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setSoftcoreAlpha");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce *target, double alpha)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setSoftcoreAlpha$descriptor() {
    return OpenMM_AmoebaVdwForce_setSoftcoreAlpha.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce *target, double alpha)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setSoftcoreAlpha$handle() {
    return OpenMM_AmoebaVdwForce_setSoftcoreAlpha.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce *target, double alpha)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setSoftcoreAlpha$address() {
    return OpenMM_AmoebaVdwForce_setSoftcoreAlpha.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(OpenMM_AmoebaVdwForce *target, double alpha)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setSoftcoreAlpha(MemorySegment target, double alpha) {
    var mh$ = OpenMM_AmoebaVdwForce_setSoftcoreAlpha.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setSoftcoreAlpha", target, alpha);
      }
      mh$.invokeExact(target, alpha);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getSoftcoreAlpha {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getSoftcoreAlpha");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getSoftcoreAlpha$descriptor() {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getSoftcoreAlpha$handle() {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getSoftcoreAlpha$address() {
    return OpenMM_AmoebaVdwForce_getSoftcoreAlpha.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static double OpenMM_AmoebaVdwForce_getSoftcoreAlpha(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getSoftcoreAlpha.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getSoftcoreAlpha", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_getAlchemicalMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_getAlchemicalMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_getAlchemicalMethod$descriptor() {
    return OpenMM_AmoebaVdwForce_getAlchemicalMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_getAlchemicalMethod$handle() {
    return OpenMM_AmoebaVdwForce_getAlchemicalMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_getAlchemicalMethod$address() {
    return OpenMM_AmoebaVdwForce_getAlchemicalMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_AmoebaVdwForce_AlchemicalMethod OpenMM_AmoebaVdwForce_getAlchemicalMethod(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_getAlchemicalMethod(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_getAlchemicalMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_getAlchemicalMethod", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_setAlchemicalMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_setAlchemicalMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_AlchemicalMethod method)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_setAlchemicalMethod$descriptor() {
    return OpenMM_AmoebaVdwForce_setAlchemicalMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_AlchemicalMethod method)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_setAlchemicalMethod$handle() {
    return OpenMM_AmoebaVdwForce_setAlchemicalMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_AlchemicalMethod method)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_setAlchemicalMethod$address() {
    return OpenMM_AmoebaVdwForce_setAlchemicalMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_setAlchemicalMethod(OpenMM_AmoebaVdwForce *target, OpenMM_AmoebaVdwForce_AlchemicalMethod method)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_setAlchemicalMethod(MemorySegment target, int method) {
    var mh$ = OpenMM_AmoebaVdwForce_setAlchemicalMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_setAlchemicalMethod", target, method);
      }
      mh$.invokeExact(target, method);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_updateParametersInContext$descriptor() {
    return OpenMM_AmoebaVdwForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_updateParametersInContext$handle() {
    return OpenMM_AmoebaVdwForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_updateParametersInContext$address() {
    return OpenMM_AmoebaVdwForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaVdwForce_updateParametersInContext(OpenMM_AmoebaVdwForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_AmoebaVdwForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_AmoebaVdwForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaVdwForce *target)
   *}
   */
  public static int OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaWcaDispersionForce *OpenMM_AmoebaWcaDispersionForce_create()
   *}
   */
  public static class OpenMM_AmoebaWcaDispersionForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_AmoebaWcaDispersionForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_AmoebaWcaDispersionForce *OpenMM_AmoebaWcaDispersionForce_create()
     *}
     */
    public static OpenMM_AmoebaWcaDispersionForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_AmoebaWcaDispersionForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_AmoebaWcaDispersionForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_destroy$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_destroy$handle() {
    return OpenMM_AmoebaWcaDispersionForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_destroy$address() {
    return OpenMM_AmoebaWcaDispersionForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_destroy(OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getNumParticles$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getNumParticles$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getNumParticles$address() {
    return OpenMM_AmoebaWcaDispersionForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_getNumParticles(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static int OpenMM_AmoebaWcaDispersionForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double radius, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setParticleParameters$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double radius, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setParticleParameters$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double radius, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setParticleParameters$address() {
    return OpenMM_AmoebaWcaDispersionForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double radius, double epsilon)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setParticleParameters(MemorySegment target, int particleIndex, double radius, double epsilon) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setParticleParameters", target, particleIndex, radius, epsilon);
      }
      mh$.invokeExact(target, particleIndex, radius, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double *radius, double *epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getParticleParameters$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double *radius, double *epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getParticleParameters$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double *radius, double *epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getParticleParameters$address() {
    return OpenMM_AmoebaWcaDispersionForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(const OpenMM_AmoebaWcaDispersionForce *target, int particleIndex, double *radius, double *epsilon)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_getParticleParameters(MemorySegment target, int particleIndex, MemorySegment radius, MemorySegment epsilon) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getParticleParameters", target, particleIndex, radius, epsilon);
      }
      mh$.invokeExact(target, particleIndex, radius, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce *target, double radius, double epsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_addParticle$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce *target, double radius, double epsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_addParticle$handle() {
    return OpenMM_AmoebaWcaDispersionForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce *target, double radius, double epsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_addParticle$address() {
    return OpenMM_AmoebaWcaDispersionForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaWcaDispersionForce_addParticle(OpenMM_AmoebaWcaDispersionForce *target, double radius, double epsilon)
   *}
   */
  public static int OpenMM_AmoebaWcaDispersionForce_addParticle(MemorySegment target, double radius, double epsilon) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_addParticle", target, radius, epsilon);
      }
      return (int) mh$.invokeExact(target, radius, epsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_updateParametersInContext$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_updateParametersInContext$handle() {
    return OpenMM_AmoebaWcaDispersionForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_updateParametersInContext$address() {
    return OpenMM_AmoebaWcaDispersionForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(OpenMM_AmoebaWcaDispersionForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getEpso {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getEpso");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getEpso$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getEpso.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getEpso$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getEpso.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getEpso$address() {
    return OpenMM_AmoebaWcaDispersionForce_getEpso.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpso(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getEpso(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getEpso.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getEpso", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getEpsh {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getEpsh");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getEpsh$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getEpsh$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getEpsh$address() {
    return OpenMM_AmoebaWcaDispersionForce_getEpsh.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getEpsh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getEpsh(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getEpsh.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getEpsh", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getRmino {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getRmino");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getRmino$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getRmino.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getRmino$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getRmino.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getRmino$address() {
    return OpenMM_AmoebaWcaDispersionForce_getRmino.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRmino(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getRmino(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getRmino.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getRmino", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getRminh {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getRminh");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getRminh$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getRminh.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getRminh$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getRminh.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getRminh$address() {
    return OpenMM_AmoebaWcaDispersionForce_getRminh.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getRminh(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getRminh(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getRminh.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getRminh", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getAwater {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getAwater");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getAwater$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getAwater.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getAwater$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getAwater.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getAwater$address() {
    return OpenMM_AmoebaWcaDispersionForce_getAwater.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getAwater(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getAwater(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getAwater.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getAwater", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getShctd {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getShctd");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getShctd$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getShctd.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getShctd$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getShctd.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getShctd$address() {
    return OpenMM_AmoebaWcaDispersionForce_getShctd.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getShctd(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getShctd(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getShctd.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getShctd", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getDispoff {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getDispoff");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getDispoff$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getDispoff$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getDispoff$address() {
    return OpenMM_AmoebaWcaDispersionForce_getDispoff.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getDispoff(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getDispoff(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getDispoff.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getDispoff", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_getSlevy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_getSlevy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_getSlevy$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_getSlevy$handle() {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_getSlevy$address() {
    return OpenMM_AmoebaWcaDispersionForce_getSlevy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaWcaDispersionForce_getSlevy(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static double OpenMM_AmoebaWcaDispersionForce_getSlevy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_getSlevy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_getSlevy", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setEpso {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setEpso");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setEpso$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setEpso.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setEpso$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setEpso.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setEpso$address() {
    return OpenMM_AmoebaWcaDispersionForce_setEpso.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpso(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setEpso(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setEpso.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setEpso", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setEpsh {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setEpsh");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setEpsh$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setEpsh.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setEpsh$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setEpsh.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setEpsh$address() {
    return OpenMM_AmoebaWcaDispersionForce_setEpsh.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setEpsh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setEpsh(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setEpsh.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setEpsh", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setRmino {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setRmino");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setRmino$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setRmino.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setRmino$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setRmino.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setRmino$address() {
    return OpenMM_AmoebaWcaDispersionForce_setRmino.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRmino(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setRmino(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setRmino.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setRmino", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setRminh {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setRminh");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setRminh$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setRminh.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setRminh$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setRminh.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setRminh$address() {
    return OpenMM_AmoebaWcaDispersionForce_setRminh.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setRminh(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setRminh(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setRminh.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setRminh", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setAwater {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setAwater");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setAwater$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setAwater.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setAwater$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setAwater.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setAwater$address() {
    return OpenMM_AmoebaWcaDispersionForce_setAwater.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setAwater(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setAwater(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setAwater.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setAwater", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setShctd {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setShctd");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setShctd$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setShctd.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setShctd$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setShctd.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setShctd$address() {
    return OpenMM_AmoebaWcaDispersionForce_setShctd.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setShctd(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setShctd(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setShctd.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setShctd", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setDispoff {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setDispoff");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setDispoff$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setDispoff.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setDispoff$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setDispoff.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setDispoff$address() {
    return OpenMM_AmoebaWcaDispersionForce_setDispoff.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setDispoff(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setDispoff(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setDispoff.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setDispoff", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_setSlevy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_setSlevy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_setSlevy$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_setSlevy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_setSlevy$handle() {
    return OpenMM_AmoebaWcaDispersionForce_setSlevy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_setSlevy$address() {
    return OpenMM_AmoebaWcaDispersionForce_setSlevy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaWcaDispersionForce_setSlevy(OpenMM_AmoebaWcaDispersionForce *target, double inputValue)
   *}
   */
  public static void OpenMM_AmoebaWcaDispersionForce_setSlevy(MemorySegment target, double inputValue) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_setSlevy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_setSlevy", target, inputValue);
      }
      mh$.invokeExact(target, inputValue);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaWcaDispersionForce *target)
   *}
   */
  public static int OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaWcaDispersionForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static final int OpenMM_HippoNonbondedForce_NoCutoff = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_NoCutoff = 0
   *}
   */
  public static int OpenMM_HippoNonbondedForce_NoCutoff() {
    return OpenMM_HippoNonbondedForce_NoCutoff;
  }

  private static final int OpenMM_HippoNonbondedForce_PME = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_PME = 1
   *}
   */
  public static int OpenMM_HippoNonbondedForce_PME() {
    return OpenMM_HippoNonbondedForce_PME;
  }

  private static final int OpenMM_HippoNonbondedForce_ZThenX = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_ZThenX = 0
   *}
   */
  public static int OpenMM_HippoNonbondedForce_ZThenX() {
    return OpenMM_HippoNonbondedForce_ZThenX;
  }

  private static final int OpenMM_HippoNonbondedForce_Bisector = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_Bisector = 1
   *}
   */
  public static int OpenMM_HippoNonbondedForce_Bisector() {
    return OpenMM_HippoNonbondedForce_Bisector;
  }

  private static final int OpenMM_HippoNonbondedForce_ZBisect = (int) 2L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_ZBisect = 2
   *}
   */
  public static int OpenMM_HippoNonbondedForce_ZBisect() {
    return OpenMM_HippoNonbondedForce_ZBisect;
  }

  private static final int OpenMM_HippoNonbondedForce_ThreeFold = (int) 3L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_ThreeFold = 3
   *}
   */
  public static int OpenMM_HippoNonbondedForce_ThreeFold() {
    return OpenMM_HippoNonbondedForce_ThreeFold;
  }

  private static final int OpenMM_HippoNonbondedForce_ZOnly = (int) 4L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_ZOnly = 4
   *}
   */
  public static int OpenMM_HippoNonbondedForce_ZOnly() {
    return OpenMM_HippoNonbondedForce_ZOnly;
  }

  private static final int OpenMM_HippoNonbondedForce_NoAxisType = (int) 5L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_HippoNonbondedForce_NoAxisType = 5
   *}
   */
  public static int OpenMM_HippoNonbondedForce_NoAxisType() {
    return OpenMM_HippoNonbondedForce_NoAxisType;
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_HippoNonbondedForce *OpenMM_HippoNonbondedForce_create()
   *}
   */
  public static class OpenMM_HippoNonbondedForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_HippoNonbondedForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_HippoNonbondedForce *OpenMM_HippoNonbondedForce_create()
     *}
     */
    public static OpenMM_HippoNonbondedForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_HippoNonbondedForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_HippoNonbondedForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_HippoNonbondedForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_destroy$descriptor() {
    return OpenMM_HippoNonbondedForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_destroy$handle() {
    return OpenMM_HippoNonbondedForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_destroy$address() {
    return OpenMM_HippoNonbondedForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_destroy(OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getNumParticles$descriptor() {
    return OpenMM_HippoNonbondedForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getNumParticles$handle() {
    return OpenMM_HippoNonbondedForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getNumParticles$address() {
    return OpenMM_HippoNonbondedForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumParticles(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getNumExceptions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getNumExceptions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getNumExceptions$descriptor() {
    return OpenMM_HippoNonbondedForce_getNumExceptions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getNumExceptions$handle() {
    return OpenMM_HippoNonbondedForce_getNumExceptions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getNumExceptions$address() {
    return OpenMM_HippoNonbondedForce_getNumExceptions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_getNumExceptions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_getNumExceptions(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getNumExceptions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getNumExceptions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getNonbondedMethod$descriptor() {
    return OpenMM_HippoNonbondedForce_getNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getNonbondedMethod$handle() {
    return OpenMM_HippoNonbondedForce_getNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getNonbondedMethod$address() {
    return OpenMM_HippoNonbondedForce_getNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_HippoNonbondedForce_NonbondedMethod OpenMM_HippoNonbondedForce_getNonbondedMethod(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_getNonbondedMethod(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getNonbondedMethod", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce *target, OpenMM_HippoNonbondedForce_NonbondedMethod method)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setNonbondedMethod$descriptor() {
    return OpenMM_HippoNonbondedForce_setNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce *target, OpenMM_HippoNonbondedForce_NonbondedMethod method)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setNonbondedMethod$handle() {
    return OpenMM_HippoNonbondedForce_setNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce *target, OpenMM_HippoNonbondedForce_NonbondedMethod method)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setNonbondedMethod$address() {
    return OpenMM_HippoNonbondedForce_setNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setNonbondedMethod(OpenMM_HippoNonbondedForce *target, OpenMM_HippoNonbondedForce_NonbondedMethod method)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setNonbondedMethod(MemorySegment target, int method) {
    var mh$ = OpenMM_HippoNonbondedForce_setNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setNonbondedMethod", target, method);
      }
      mh$.invokeExact(target, method);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getCutoffDistance$descriptor() {
    return OpenMM_HippoNonbondedForce_getCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getCutoffDistance$handle() {
    return OpenMM_HippoNonbondedForce_getCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getCutoffDistance$address() {
    return OpenMM_HippoNonbondedForce_getCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getCutoffDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static double OpenMM_HippoNonbondedForce_getCutoffDistance(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getCutoffDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setCutoffDistance$descriptor() {
    return OpenMM_HippoNonbondedForce_setCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setCutoffDistance$handle() {
    return OpenMM_HippoNonbondedForce_setCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setCutoffDistance$address() {
    return OpenMM_HippoNonbondedForce_setCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setCutoffDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setCutoffDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_HippoNonbondedForce_setCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setCutoffDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getSwitchingDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getSwitchingDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getSwitchingDistance$descriptor() {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getSwitchingDistance$handle() {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getSwitchingDistance$address() {
    return OpenMM_HippoNonbondedForce_getSwitchingDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getSwitchingDistance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static double OpenMM_HippoNonbondedForce_getSwitchingDistance(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getSwitchingDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getSwitchingDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setSwitchingDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setSwitchingDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setSwitchingDistance$descriptor() {
    return OpenMM_HippoNonbondedForce_setSwitchingDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setSwitchingDistance$handle() {
    return OpenMM_HippoNonbondedForce_setSwitchingDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setSwitchingDistance$address() {
    return OpenMM_HippoNonbondedForce_setSwitchingDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setSwitchingDistance(OpenMM_HippoNonbondedForce *target, double distance)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setSwitchingDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_HippoNonbondedForce_setSwitchingDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setSwitchingDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getExtrapolationCoefficients {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getExtrapolationCoefficients");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getExtrapolationCoefficients$descriptor() {
    return OpenMM_HippoNonbondedForce_getExtrapolationCoefficients.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getExtrapolationCoefficients$handle() {
    return OpenMM_HippoNonbondedForce_getExtrapolationCoefficients.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getExtrapolationCoefficients$address() {
    return OpenMM_HippoNonbondedForce_getExtrapolationCoefficients.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getExtrapolationCoefficients(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getExtrapolationCoefficients.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getExtrapolationCoefficients", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setExtrapolationCoefficients {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setExtrapolationCoefficients");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setExtrapolationCoefficients$descriptor() {
    return OpenMM_HippoNonbondedForce_setExtrapolationCoefficients.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setExtrapolationCoefficients$handle() {
    return OpenMM_HippoNonbondedForce_setExtrapolationCoefficients.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setExtrapolationCoefficients$address() {
    return OpenMM_HippoNonbondedForce_setExtrapolationCoefficients.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(OpenMM_HippoNonbondedForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setExtrapolationCoefficients(MemorySegment target, MemorySegment coefficients) {
    var mh$ = OpenMM_HippoNonbondedForce_setExtrapolationCoefficients.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setExtrapolationCoefficients", target, coefficients);
      }
      mh$.invokeExact(target, coefficients);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getPMEParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_getPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getPMEParameters$handle() {
    return OpenMM_HippoNonbondedForce_getPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getPMEParameters$address() {
    return OpenMM_HippoNonbondedForce_getPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getPMEParameters(MemorySegment target, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_HippoNonbondedForce_getPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getDPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getDPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getDPMEParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_getDPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getDPMEParameters$handle() {
    return OpenMM_HippoNonbondedForce_getDPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getDPMEParameters$address() {
    return OpenMM_HippoNonbondedForce_getDPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParameters(const OpenMM_HippoNonbondedForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getDPMEParameters(MemorySegment target, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_HippoNonbondedForce_getDPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getDPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setPMEParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_setPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setPMEParameters$handle() {
    return OpenMM_HippoNonbondedForce_setPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setPMEParameters$address() {
    return OpenMM_HippoNonbondedForce_setPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setPMEParameters(MemorySegment target, double alpha, int nx, int ny, int nz) {
    var mh$ = OpenMM_HippoNonbondedForce_setPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setDPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setDPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setDPMEParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_setDPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setDPMEParameters$handle() {
    return OpenMM_HippoNonbondedForce_setDPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setDPMEParameters$address() {
    return OpenMM_HippoNonbondedForce_setDPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setDPMEParameters(OpenMM_HippoNonbondedForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setDPMEParameters(MemorySegment target, double alpha, int nx, int ny, int nz) {
    var mh$ = OpenMM_HippoNonbondedForce_setDPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setDPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getPMEParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getPMEParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getPMEParametersInContext$descriptor() {
    return OpenMM_HippoNonbondedForce_getPMEParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getPMEParametersInContext$handle() {
    return OpenMM_HippoNonbondedForce_getPMEParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getPMEParametersInContext$address() {
    return OpenMM_HippoNonbondedForce_getPMEParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getPMEParametersInContext(MemorySegment target, MemorySegment context, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_HippoNonbondedForce_getPMEParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getPMEParametersInContext", target, context, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, context, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getDPMEParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getDPMEParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getDPMEParametersInContext$descriptor() {
    return OpenMM_HippoNonbondedForce_getDPMEParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getDPMEParametersInContext$handle() {
    return OpenMM_HippoNonbondedForce_getDPMEParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getDPMEParametersInContext$address() {
    return OpenMM_HippoNonbondedForce_getDPMEParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(const OpenMM_HippoNonbondedForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getDPMEParametersInContext(MemorySegment target, MemorySegment context, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_HippoNonbondedForce_getDPMEParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getDPMEParametersInContext", target, context, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, context, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce *target, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_addParticle$descriptor() {
    return OpenMM_HippoNonbondedForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce *target, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_addParticle$handle() {
    return OpenMM_HippoNonbondedForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce *target, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_addParticle$address() {
    return OpenMM_HippoNonbondedForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addParticle(OpenMM_HippoNonbondedForce *target, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_addParticle(MemorySegment target, double charge, MemorySegment dipole, MemorySegment quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    var mh$ = OpenMM_HippoNonbondedForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_addParticle", target, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
      }
      return (int) mh$.invokeExact(target, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce *target, int index, double *charge, OpenMM_DoubleArray *dipole, OpenMM_DoubleArray *quadrupole, double *coreCharge, double *alpha, double *epsilon, double *damping, double *c6, double *pauliK, double *pauliQ, double *pauliAlpha, double *polarizability, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getParticleParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce *target, int index, double *charge, OpenMM_DoubleArray *dipole, OpenMM_DoubleArray *quadrupole, double *coreCharge, double *alpha, double *epsilon, double *damping, double *c6, double *pauliK, double *pauliQ, double *pauliAlpha, double *polarizability, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getParticleParameters$handle() {
    return OpenMM_HippoNonbondedForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce *target, int index, double *charge, OpenMM_DoubleArray *dipole, OpenMM_DoubleArray *quadrupole, double *coreCharge, double *alpha, double *epsilon, double *damping, double *c6, double *pauliK, double *pauliQ, double *pauliAlpha, double *polarizability, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getParticleParameters$address() {
    return OpenMM_HippoNonbondedForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getParticleParameters(const OpenMM_HippoNonbondedForce *target, int index, double *charge, OpenMM_DoubleArray *dipole, OpenMM_DoubleArray *quadrupole, double *coreCharge, double *alpha, double *epsilon, double *damping, double *c6, double *pauliK, double *pauliQ, double *pauliAlpha, double *polarizability, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getParticleParameters(MemorySegment target, int index, MemorySegment charge, MemorySegment dipole, MemorySegment quadrupole, MemorySegment coreCharge, MemorySegment alpha, MemorySegment epsilon, MemorySegment damping, MemorySegment c6, MemorySegment pauliK, MemorySegment pauliQ, MemorySegment pauliAlpha, MemorySegment polarizability, MemorySegment axisType, MemorySegment multipoleAtomZ, MemorySegment multipoleAtomX, MemorySegment multipoleAtomY) {
    var mh$ = OpenMM_HippoNonbondedForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getParticleParameters", target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
      }
      mh$.invokeExact(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce *target, int index, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setParticleParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce *target, int index, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setParticleParameters$handle() {
    return OpenMM_HippoNonbondedForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce *target, int index, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setParticleParameters$address() {
    return OpenMM_HippoNonbondedForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setParticleParameters(OpenMM_HippoNonbondedForce *target, int index, double charge, const OpenMM_DoubleArray *dipole, const OpenMM_DoubleArray *quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setParticleParameters(MemorySegment target, int index, double charge, MemorySegment dipole, MemorySegment quadrupole, double coreCharge, double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha, double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) {
    var mh$ = OpenMM_HippoNonbondedForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setParticleParameters", target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
      }
      mh$.invokeExact(target, index, charge, dipole, quadrupole, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_addException {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_addException");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce *target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_addException$descriptor() {
    return OpenMM_HippoNonbondedForce_addException.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce *target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_addException$handle() {
    return OpenMM_HippoNonbondedForce_addException.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce *target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_addException$address() {
    return OpenMM_HippoNonbondedForce_addException.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_HippoNonbondedForce_addException(OpenMM_HippoNonbondedForce *target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, OpenMM_Boolean replace)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_addException(MemorySegment target, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale, int replace) {
    var mh$ = OpenMM_HippoNonbondedForce_addException.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_addException", target, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
      }
      return (int) mh$.invokeExact(target, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale, replace);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getExceptionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getExceptionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce *target, int index, int *particle1, int *particle2, double *multipoleMultipoleScale, double *dipoleMultipoleScale, double *dipoleDipoleScale, double *dispersionScale, double *repulsionScale, double *chargeTransferScale)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getExceptionParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_getExceptionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce *target, int index, int *particle1, int *particle2, double *multipoleMultipoleScale, double *dipoleMultipoleScale, double *dipoleDipoleScale, double *dispersionScale, double *repulsionScale, double *chargeTransferScale)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getExceptionParameters$handle() {
    return OpenMM_HippoNonbondedForce_getExceptionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce *target, int index, int *particle1, int *particle2, double *multipoleMultipoleScale, double *dipoleMultipoleScale, double *dipoleDipoleScale, double *dispersionScale, double *repulsionScale, double *chargeTransferScale)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getExceptionParameters$address() {
    return OpenMM_HippoNonbondedForce_getExceptionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getExceptionParameters(const OpenMM_HippoNonbondedForce *target, int index, int *particle1, int *particle2, double *multipoleMultipoleScale, double *dipoleMultipoleScale, double *dipoleDipoleScale, double *dispersionScale, double *repulsionScale, double *chargeTransferScale)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getExceptionParameters(MemorySegment target, int index, MemorySegment particle1, MemorySegment particle2, MemorySegment multipoleMultipoleScale, MemorySegment dipoleMultipoleScale, MemorySegment dipoleDipoleScale, MemorySegment dispersionScale, MemorySegment repulsionScale, MemorySegment chargeTransferScale) {
    var mh$ = OpenMM_HippoNonbondedForce_getExceptionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getExceptionParameters", target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
      }
      mh$.invokeExact(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setExceptionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setExceptionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce *target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setExceptionParameters$descriptor() {
    return OpenMM_HippoNonbondedForce_setExceptionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce *target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setExceptionParameters$handle() {
    return OpenMM_HippoNonbondedForce_setExceptionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce *target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setExceptionParameters$address() {
    return OpenMM_HippoNonbondedForce_setExceptionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setExceptionParameters(OpenMM_HippoNonbondedForce *target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setExceptionParameters(MemorySegment target, int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale) {
    var mh$ = OpenMM_HippoNonbondedForce_setExceptionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setExceptionParameters", target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
      }
      mh$.invokeExact(target, index, particle1, particle2, multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getEwaldErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getEwaldErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getEwaldErrorTolerance$descriptor() {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getEwaldErrorTolerance$handle() {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getEwaldErrorTolerance$address() {
    return OpenMM_HippoNonbondedForce_getEwaldErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static double OpenMM_HippoNonbondedForce_getEwaldErrorTolerance(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_getEwaldErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getEwaldErrorTolerance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_setEwaldErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_setEwaldErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce *target, double tol)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_setEwaldErrorTolerance$descriptor() {
    return OpenMM_HippoNonbondedForce_setEwaldErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce *target, double tol)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_setEwaldErrorTolerance$handle() {
    return OpenMM_HippoNonbondedForce_setEwaldErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce *target, double tol)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_setEwaldErrorTolerance$address() {
    return OpenMM_HippoNonbondedForce_setEwaldErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(OpenMM_HippoNonbondedForce *target, double tol)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_setEwaldErrorTolerance(MemorySegment target, double tol) {
    var mh$ = OpenMM_HippoNonbondedForce_setEwaldErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_setEwaldErrorTolerance", target, tol);
      }
      mh$.invokeExact(target, tol);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles$descriptor() {
    return OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles$handle() {
    return OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles$address() {
    return OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles(MemorySegment target, MemorySegment context, MemorySegment dipoles) {
    var mh$ = OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getLabFramePermanentDipoles", target, context, dipoles);
      }
      mh$.invokeExact(target, context, dipoles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_getInducedDipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_getInducedDipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_getInducedDipoles$descriptor() {
    return OpenMM_HippoNonbondedForce_getInducedDipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_getInducedDipoles$handle() {
    return OpenMM_HippoNonbondedForce_getInducedDipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_getInducedDipoles$address() {
    return OpenMM_HippoNonbondedForce_getInducedDipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_getInducedDipoles(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_getInducedDipoles(MemorySegment target, MemorySegment context, MemorySegment dipoles) {
    var mh$ = OpenMM_HippoNonbondedForce_getInducedDipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_getInducedDipoles", target, context, dipoles);
      }
      mh$.invokeExact(target, context, dipoles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_updateParametersInContext$descriptor() {
    return OpenMM_HippoNonbondedForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_updateParametersInContext$handle() {
    return OpenMM_HippoNonbondedForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_updateParametersInContext$address() {
    return OpenMM_HippoNonbondedForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_HippoNonbondedForce_updateParametersInContext(OpenMM_HippoNonbondedForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_HippoNonbondedForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_HippoNonbondedForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MethodHandle OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static MemorySegment OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(const OpenMM_HippoNonbondedForce *target)
   *}
   */
  public static int OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_HippoNonbondedForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaGeneralizedKirkwoodForce *OpenMM_AmoebaGeneralizedKirkwoodForce_create()
   *}
   */
  public static class OpenMM_AmoebaGeneralizedKirkwoodForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_AmoebaGeneralizedKirkwoodForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_AmoebaGeneralizedKirkwoodForce *OpenMM_AmoebaGeneralizedKirkwoodForce_create()
     *}
     */
    public static OpenMM_AmoebaGeneralizedKirkwoodForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_AmoebaGeneralizedKirkwoodForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_destroy$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_destroy$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_destroy$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(MemorySegment target, double charge, double radius, double scalingFactor) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle", target, charge, radius, scalingFactor);
      }
      return (int) mh$.invokeExact(target, charge, radius, scalingFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1 {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(MemorySegment target, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1", target, charge, radius, scalingFactor, descreenRadius, neckFactor);
      }
      return (int) mh$.invokeExact(target, charge, radius, scalingFactor, descreenRadius, neckFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double *charge, double *radius, double *scalingFactor, double *descreenRadius, double *neckFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double *charge, double *radius, double *scalingFactor, double *descreenRadius, double *neckFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double *charge, double *radius, double *scalingFactor, double *descreenRadius, double *neckFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double *charge, double *radius, double *scalingFactor, double *descreenRadius, double *neckFactor)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters(MemorySegment target, int index, MemorySegment charge, MemorySegment radius, MemorySegment scalingFactor, MemorySegment descreenRadius, MemorySegment neckFactor) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getParticleParameters", target, index, charge, radius, scalingFactor, descreenRadius, neckFactor);
      }
      mh$.invokeExact(target, index, charge, radius, scalingFactor, descreenRadius, neckFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int index, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(MemorySegment target, int index, double charge, double radius, double scalingFactor, double descreenRadius, double neckFactor) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters", target, index, charge, radius, scalingFactor, descreenRadius, neckFactor);
      }
      mh$.invokeExact(target, index, charge, radius, scalingFactor, descreenRadius, neckFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getSolventDielectric", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(MemorySegment target, double dielectric) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric", target, dielectric);
      }
      mh$.invokeExact(target, dielectric);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getSoluteDielectric", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectric)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(MemorySegment target, double dielectric) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric", target, dielectric);
      }
      mh$.invokeExact(target, dielectric);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhRescaling", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Boolean tanhRescale)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Boolean tanhRescale)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Boolean tanhRescale)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Boolean tanhRescale)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(MemorySegment target, int tanhRescale) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling", target, tanhRescale);
      }
      mh$.invokeExact(target, tanhRescale);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, double *b0, double *b1, double *b2)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, double *b0, double *b1, double *b2)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, double *b0, double *b1, double *b2)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(const OpenMM_AmoebaGeneralizedKirkwoodForce *target, double *b0, double *b1, double *b2)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters(MemorySegment target, MemorySegment b0, MemorySegment b1, MemorySegment b2) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getTanhParameters", target, b0, b1, b2);
      }
      mh$.invokeExact(target, b0, b1, b2);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double b0, double b1, double b2)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double b0, double b1, double b2)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double b0, double b1, double b2)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double b0, double b1, double b2)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(MemorySegment target, double b0, double b1, double b2) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters", target, b0, b1, b2);
      }
      mh$.invokeExact(target, b0, b1, b2);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getDescreenOffset", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double descreenOffet)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double descreenOffet)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double descreenOffet)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double descreenOffet)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset(MemorySegment target, double descreenOffet) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setDescreenOffset", target, descreenOffet);
      }
      mh$.invokeExact(target, descreenOffet);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getIncludeCavityTerm", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int includeCavityTerm)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int includeCavityTerm)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int includeCavityTerm)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(OpenMM_AmoebaGeneralizedKirkwoodForce *target, int includeCavityTerm)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(MemorySegment target, int includeCavityTerm) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm", target, includeCavityTerm);
      }
      mh$.invokeExact(target, includeCavityTerm);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getProbeRadius", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double probeRadius)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double probeRadius)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double probeRadius)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double probeRadius)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(MemorySegment target, double probeRadius) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius", target, probeRadius);
      }
      mh$.invokeExact(target, probeRadius);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getDielectricOffset", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectricOffset)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectricOffset)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectricOffset)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double dielectricOffset)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(MemorySegment target, double dielectricOffset) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset", target, dielectricOffset);
      }
      mh$.invokeExact(target, dielectricOffset);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static double OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_getSurfaceAreaFactor", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double surfaceAreaFactor)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double surfaceAreaFactor)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double surfaceAreaFactor)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(OpenMM_AmoebaGeneralizedKirkwoodForce *target, double surfaceAreaFactor)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(MemorySegment target, double surfaceAreaFactor) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor", target, surfaceAreaFactor);
      }
      mh$.invokeExact(target, surfaceAreaFactor);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(OpenMM_AmoebaGeneralizedKirkwoodForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaGeneralizedKirkwoodForce *target)
   *}
   */
  public static int OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaGeneralizedKirkwoodForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static final int OpenMM_AmoebaMultipoleForce_NoCutoff = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_NoCutoff = 0
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_NoCutoff() {
    return OpenMM_AmoebaMultipoleForce_NoCutoff;
  }

  private static final int OpenMM_AmoebaMultipoleForce_PME = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_PME = 1
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_PME() {
    return OpenMM_AmoebaMultipoleForce_PME;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Mutual = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Mutual = 0
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Mutual() {
    return OpenMM_AmoebaMultipoleForce_Mutual;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Direct = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Direct = 1
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Direct() {
    return OpenMM_AmoebaMultipoleForce_Direct;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Extrapolated = (int) 2L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Extrapolated = 2
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Extrapolated() {
    return OpenMM_AmoebaMultipoleForce_Extrapolated;
  }

  private static final int OpenMM_AmoebaMultipoleForce_ZThenX = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_ZThenX = 0
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_ZThenX() {
    return OpenMM_AmoebaMultipoleForce_ZThenX;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Bisector = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Bisector = 1
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Bisector() {
    return OpenMM_AmoebaMultipoleForce_Bisector;
  }

  private static final int OpenMM_AmoebaMultipoleForce_ZBisect = (int) 2L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_ZBisect = 2
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_ZBisect() {
    return OpenMM_AmoebaMultipoleForce_ZBisect;
  }

  private static final int OpenMM_AmoebaMultipoleForce_ThreeFold = (int) 3L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_ThreeFold = 3
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_ThreeFold() {
    return OpenMM_AmoebaMultipoleForce_ThreeFold;
  }

  private static final int OpenMM_AmoebaMultipoleForce_ZOnly = (int) 4L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_ZOnly = 4
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_ZOnly() {
    return OpenMM_AmoebaMultipoleForce_ZOnly;
  }

  private static final int OpenMM_AmoebaMultipoleForce_NoAxisType = (int) 5L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_NoAxisType = 5
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_NoAxisType() {
    return OpenMM_AmoebaMultipoleForce_NoAxisType;
  }

  private static final int OpenMM_AmoebaMultipoleForce_LastAxisTypeIndex = (int) 6L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_LastAxisTypeIndex = 6
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_LastAxisTypeIndex() {
    return OpenMM_AmoebaMultipoleForce_LastAxisTypeIndex;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Covalent12 = (int) 0L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Covalent12 = 0
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Covalent12() {
    return OpenMM_AmoebaMultipoleForce_Covalent12;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Covalent13 = (int) 1L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Covalent13 = 1
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Covalent13() {
    return OpenMM_AmoebaMultipoleForce_Covalent13;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Covalent14 = (int) 2L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Covalent14 = 2
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Covalent14() {
    return OpenMM_AmoebaMultipoleForce_Covalent14;
  }

  private static final int OpenMM_AmoebaMultipoleForce_Covalent15 = (int) 3L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_Covalent15 = 3
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_Covalent15() {
    return OpenMM_AmoebaMultipoleForce_Covalent15;
  }

  private static final int OpenMM_AmoebaMultipoleForce_PolarizationCovalent11 = (int) 4L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11 = 4
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_PolarizationCovalent11() {
    return OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
  }

  private static final int OpenMM_AmoebaMultipoleForce_PolarizationCovalent12 = (int) 5L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_PolarizationCovalent12 = 5
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_PolarizationCovalent12() {
    return OpenMM_AmoebaMultipoleForce_PolarizationCovalent12;
  }

  private static final int OpenMM_AmoebaMultipoleForce_PolarizationCovalent13 = (int) 6L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_PolarizationCovalent13 = 6
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_PolarizationCovalent13() {
    return OpenMM_AmoebaMultipoleForce_PolarizationCovalent13;
  }

  private static final int OpenMM_AmoebaMultipoleForce_PolarizationCovalent14 = (int) 7L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_PolarizationCovalent14 = 7
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_PolarizationCovalent14() {
    return OpenMM_AmoebaMultipoleForce_PolarizationCovalent14;
  }

  private static final int OpenMM_AmoebaMultipoleForce_CovalentEnd = (int) 8L;

  /**
   * {@snippet lang = c:
   * enum <anonymous>.OpenMM_AmoebaMultipoleForce_CovalentEnd = 8
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_CovalentEnd() {
    return OpenMM_AmoebaMultipoleForce_CovalentEnd;
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce *OpenMM_AmoebaMultipoleForce_create()
   *}
   */
  public static class OpenMM_AmoebaMultipoleForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_AmoebaMultipoleForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_AmoebaMultipoleForce *OpenMM_AmoebaMultipoleForce_create()
     *}
     */
    public static OpenMM_AmoebaMultipoleForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_AmoebaMultipoleForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_AmoebaMultipoleForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_destroy$descriptor() {
    return OpenMM_AmoebaMultipoleForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_destroy$handle() {
    return OpenMM_AmoebaMultipoleForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_destroy$address() {
    return OpenMM_AmoebaMultipoleForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_destroy(OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getNumMultipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getNumMultipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getNumMultipoles$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getNumMultipoles$handle() {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getNumMultipoles$address() {
    return OpenMM_AmoebaMultipoleForce_getNumMultipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getNumMultipoles(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_getNumMultipoles(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getNumMultipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getNumMultipoles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getNonbondedMethod$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getNonbondedMethod$handle() {
    return OpenMM_AmoebaMultipoleForce_getNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getNonbondedMethod$address() {
    return OpenMM_AmoebaMultipoleForce_getNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_NonbondedMethod OpenMM_AmoebaMultipoleForce_getNonbondedMethod(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_getNonbondedMethod(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getNonbondedMethod", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setNonbondedMethod {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setNonbondedMethod");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setNonbondedMethod$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setNonbondedMethod.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setNonbondedMethod$handle() {
    return OpenMM_AmoebaMultipoleForce_setNonbondedMethod.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setNonbondedMethod$address() {
    return OpenMM_AmoebaMultipoleForce_setNonbondedMethod.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_NonbondedMethod method)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setNonbondedMethod(MemorySegment target, int method) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setNonbondedMethod.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setNonbondedMethod", target, method);
      }
      mh$.invokeExact(target, method);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getPolarizationType {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getPolarizationType");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getPolarizationType$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getPolarizationType.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getPolarizationType$handle() {
    return OpenMM_AmoebaMultipoleForce_getPolarizationType.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getPolarizationType$address() {
    return OpenMM_AmoebaMultipoleForce_getPolarizationType.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_AmoebaMultipoleForce_PolarizationType OpenMM_AmoebaMultipoleForce_getPolarizationType(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_getPolarizationType(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getPolarizationType.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getPolarizationType", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setPolarizationType {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setPolarizationType");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_PolarizationType type)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setPolarizationType$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setPolarizationType.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_PolarizationType type)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setPolarizationType$handle() {
    return OpenMM_AmoebaMultipoleForce_setPolarizationType.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_PolarizationType type)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setPolarizationType$address() {
    return OpenMM_AmoebaMultipoleForce_setPolarizationType.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPolarizationType(OpenMM_AmoebaMultipoleForce *target, OpenMM_AmoebaMultipoleForce_PolarizationType type)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setPolarizationType(MemorySegment target, int type) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setPolarizationType.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setPolarizationType", target, type);
      }
      mh$.invokeExact(target, type);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getCutoffDistance$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getCutoffDistance$handle() {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getCutoffDistance$address() {
    return OpenMM_AmoebaMultipoleForce_getCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getCutoffDistance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static double OpenMM_AmoebaMultipoleForce_getCutoffDistance(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getCutoffDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setCutoffDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setCutoffDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setCutoffDistance$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setCutoffDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setCutoffDistance$handle() {
    return OpenMM_AmoebaMultipoleForce_setCutoffDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setCutoffDistance$address() {
    return OpenMM_AmoebaMultipoleForce_setCutoffDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCutoffDistance(OpenMM_AmoebaMultipoleForce *target, double distance)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setCutoffDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setCutoffDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setCutoffDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getPMEParameters$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getPMEParameters$handle() {
    return OpenMM_AmoebaMultipoleForce_getPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getPMEParameters$address() {
    return OpenMM_AmoebaMultipoleForce_getPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParameters(const OpenMM_AmoebaMultipoleForce *target, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getPMEParameters(MemorySegment target, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setPMEParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setPMEParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setPMEParameters$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setPMEParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setPMEParameters$handle() {
    return OpenMM_AmoebaMultipoleForce_setPMEParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setPMEParameters$address() {
    return OpenMM_AmoebaMultipoleForce_setPMEParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPMEParameters(OpenMM_AmoebaMultipoleForce *target, double alpha, int nx, int ny, int nz)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setPMEParameters(MemorySegment target, double alpha, int nx, int ny, int nz) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setPMEParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setPMEParameters", target, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getAEwald {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getAEwald");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getAEwald$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getAEwald.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getAEwald$handle() {
    return OpenMM_AmoebaMultipoleForce_getAEwald.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getAEwald$address() {
    return OpenMM_AmoebaMultipoleForce_getAEwald.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getAEwald(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static double OpenMM_AmoebaMultipoleForce_getAEwald(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getAEwald.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getAEwald", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setAEwald {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setAEwald");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce *target, double aewald)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setAEwald$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setAEwald.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce *target, double aewald)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setAEwald$handle() {
    return OpenMM_AmoebaMultipoleForce_setAEwald.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce *target, double aewald)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setAEwald$address() {
    return OpenMM_AmoebaMultipoleForce_setAEwald.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setAEwald(OpenMM_AmoebaMultipoleForce *target, double aewald)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setAEwald(MemorySegment target, double aewald) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setAEwald.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setAEwald", target, aewald);
      }
      mh$.invokeExact(target, aewald);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder$handle() {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder$address() {
    return OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getPmeBSplineOrder", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getPmeGridDimensions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getPmeGridDimensions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce *target, OpenMM_IntArray *gridDimension)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getPmeGridDimensions$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getPmeGridDimensions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce *target, OpenMM_IntArray *gridDimension)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getPmeGridDimensions$handle() {
    return OpenMM_AmoebaMultipoleForce_getPmeGridDimensions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce *target, OpenMM_IntArray *gridDimension)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getPmeGridDimensions$address() {
    return OpenMM_AmoebaMultipoleForce_getPmeGridDimensions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(const OpenMM_AmoebaMultipoleForce *target, OpenMM_IntArray *gridDimension)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getPmeGridDimensions(MemorySegment target, MemorySegment gridDimension) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getPmeGridDimensions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getPmeGridDimensions", target, gridDimension);
      }
      mh$.invokeExact(target, gridDimension);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setPmeGridDimensions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setPmeGridDimensions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce *target, const OpenMM_IntArray *gridDimension)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setPmeGridDimensions$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setPmeGridDimensions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce *target, const OpenMM_IntArray *gridDimension)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setPmeGridDimensions$handle() {
    return OpenMM_AmoebaMultipoleForce_setPmeGridDimensions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce *target, const OpenMM_IntArray *gridDimension)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setPmeGridDimensions$address() {
    return OpenMM_AmoebaMultipoleForce_setPmeGridDimensions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(OpenMM_AmoebaMultipoleForce *target, const OpenMM_IntArray *gridDimension)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(MemorySegment target, MemorySegment gridDimension) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setPmeGridDimensions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setPmeGridDimensions", target, gridDimension);
      }
      mh$.invokeExact(target, gridDimension);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getPMEParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getPMEParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getPMEParametersInContext$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getPMEParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getPMEParametersInContext$handle() {
    return OpenMM_AmoebaMultipoleForce_getPMEParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getPMEParametersInContext$address() {
    return OpenMM_AmoebaMultipoleForce_getPMEParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(const OpenMM_AmoebaMultipoleForce *target, const OpenMM_Context *context, double *alpha, int *nx, int *ny, int *nz)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getPMEParametersInContext(MemorySegment target, MemorySegment context, MemorySegment alpha, MemorySegment nx, MemorySegment ny, MemorySegment nz) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getPMEParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getPMEParametersInContext", target, context, alpha, nx, ny, nz);
      }
      mh$.invokeExact(target, context, alpha, nx, ny, nz);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_addMultipole {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_addMultipole");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce *target, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_addMultipole$descriptor() {
    return OpenMM_AmoebaMultipoleForce_addMultipole.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce *target, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_addMultipole$handle() {
    return OpenMM_AmoebaMultipoleForce_addMultipole.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce *target, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_addMultipole$address() {
    return OpenMM_AmoebaMultipoleForce_addMultipole.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_addMultipole(OpenMM_AmoebaMultipoleForce *target, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_addMultipole(MemorySegment target, double charge, MemorySegment molecularDipole, MemorySegment molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) {
    var mh$ = OpenMM_AmoebaMultipoleForce_addMultipole.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_addMultipole", target, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
      }
      return (int) mh$.invokeExact(target, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getMultipoleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getMultipoleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce *target, int index, double *charge, OpenMM_DoubleArray *molecularDipole, OpenMM_DoubleArray *molecularQuadrupole, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY, double *thole, double *dampingFactor, double *polarity)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getMultipoleParameters$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getMultipoleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce *target, int index, double *charge, OpenMM_DoubleArray *molecularDipole, OpenMM_DoubleArray *molecularQuadrupole, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY, double *thole, double *dampingFactor, double *polarity)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getMultipoleParameters$handle() {
    return OpenMM_AmoebaMultipoleForce_getMultipoleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce *target, int index, double *charge, OpenMM_DoubleArray *molecularDipole, OpenMM_DoubleArray *molecularQuadrupole, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY, double *thole, double *dampingFactor, double *polarity)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getMultipoleParameters$address() {
    return OpenMM_AmoebaMultipoleForce_getMultipoleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(const OpenMM_AmoebaMultipoleForce *target, int index, double *charge, OpenMM_DoubleArray *molecularDipole, OpenMM_DoubleArray *molecularQuadrupole, int *axisType, int *multipoleAtomZ, int *multipoleAtomX, int *multipoleAtomY, double *thole, double *dampingFactor, double *polarity)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getMultipoleParameters(MemorySegment target, int index, MemorySegment charge, MemorySegment molecularDipole, MemorySegment molecularQuadrupole, MemorySegment axisType, MemorySegment multipoleAtomZ, MemorySegment multipoleAtomX, MemorySegment multipoleAtomY, MemorySegment thole, MemorySegment dampingFactor, MemorySegment polarity) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getMultipoleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getMultipoleParameters", target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
      }
      mh$.invokeExact(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setMultipoleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setMultipoleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce *target, int index, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setMultipoleParameters$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setMultipoleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce *target, int index, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setMultipoleParameters$handle() {
    return OpenMM_AmoebaMultipoleForce_setMultipoleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce *target, int index, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setMultipoleParameters$address() {
    return OpenMM_AmoebaMultipoleForce_setMultipoleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(OpenMM_AmoebaMultipoleForce *target, int index, double charge, const OpenMM_DoubleArray *molecularDipole, const OpenMM_DoubleArray *molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setMultipoleParameters(MemorySegment target, int index, double charge, MemorySegment molecularDipole, MemorySegment molecularQuadrupole, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY, double thole, double dampingFactor, double polarity) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setMultipoleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setMultipoleParameters", target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
      }
      mh$.invokeExact(target, index, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY, thole, dampingFactor, polarity);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setCovalentMap {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setCovalentMap");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setCovalentMap$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setCovalentMap.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setCovalentMap$handle() {
    return OpenMM_AmoebaMultipoleForce_setCovalentMap.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setCovalentMap$address() {
    return OpenMM_AmoebaMultipoleForce_setCovalentMap.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setCovalentMap(OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, const OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setCovalentMap(MemorySegment target, int index, int typeId, MemorySegment covalentAtoms) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setCovalentMap.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setCovalentMap", target, index, typeId, covalentAtoms);
      }
      mh$.invokeExact(target, index, typeId, covalentAtoms);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getCovalentMap {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getCovalentMap");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getCovalentMap$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMap.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getCovalentMap$handle() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMap.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getCovalentMap$address() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMap.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMap(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_AmoebaMultipoleForce_CovalentType typeId, OpenMM_IntArray *covalentAtoms)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getCovalentMap(MemorySegment target, int index, int typeId, MemorySegment covalentAtoms) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getCovalentMap.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getCovalentMap", target, index, typeId, covalentAtoms);
      }
      mh$.invokeExact(target, index, typeId, covalentAtoms);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getCovalentMaps {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getCovalentMaps");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_2D_IntArray *covalentLists)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getCovalentMaps$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMaps.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_2D_IntArray *covalentLists)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getCovalentMaps$handle() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMaps.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_2D_IntArray *covalentLists)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getCovalentMaps$address() {
    return OpenMM_AmoebaMultipoleForce_getCovalentMaps.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getCovalentMaps(const OpenMM_AmoebaMultipoleForce *target, int index, OpenMM_2D_IntArray *covalentLists)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getCovalentMaps(MemorySegment target, int index, MemorySegment covalentLists) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getCovalentMaps.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getCovalentMaps", target, index, covalentLists);
      }
      mh$.invokeExact(target, index, covalentLists);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations$handle() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations$address() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getMutualInducedMaxIterations", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce *target, int inputMutualInducedMaxIterations)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce *target, int inputMutualInducedMaxIterations)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations$handle() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce *target, int inputMutualInducedMaxIterations)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations$address() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(OpenMM_AmoebaMultipoleForce *target, int inputMutualInducedMaxIterations)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(MemorySegment target, int inputMutualInducedMaxIterations) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations", target, inputMutualInducedMaxIterations);
      }
      mh$.invokeExact(target, inputMutualInducedMaxIterations);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon$handle() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon$address() {
    return OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static double OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getMutualInducedTargetEpsilon", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce *target, double inputMutualInducedTargetEpsilon)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce *target, double inputMutualInducedTargetEpsilon)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon$handle() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce *target, double inputMutualInducedTargetEpsilon)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon$address() {
    return OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(OpenMM_AmoebaMultipoleForce *target, double inputMutualInducedTargetEpsilon)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(MemorySegment target, double inputMutualInducedTargetEpsilon) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon", target, inputMutualInducedTargetEpsilon);
      }
      mh$.invokeExact(target, inputMutualInducedTargetEpsilon);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients$handle() {
    return OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients$address() {
    return OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(OpenMM_AmoebaMultipoleForce *target, const OpenMM_DoubleArray *coefficients)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(MemorySegment target, MemorySegment coefficients) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients", target, coefficients);
      }
      mh$.invokeExact(target, coefficients);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients$handle() {
    return OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients$address() {
    return OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_DoubleArray *OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getExtrapolationCoefficients", target);
      }
      return (MemorySegment) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance$handle() {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance$address() {
    return OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static double OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getEwaldErrorTolerance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce *target, double tol)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance$descriptor() {
    return OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce *target, double tol)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance$handle() {
    return OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce *target, double tol)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance$address() {
    return OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(OpenMM_AmoebaMultipoleForce *target, double tol)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(MemorySegment target, double tol) {
    var mh$ = OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance", target, tol);
      }
      mh$.invokeExact(target, tol);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles$handle() {
    return OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles$address() {
    return OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles(MemorySegment target, MemorySegment context, MemorySegment dipoles) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getLabFramePermanentDipoles", target, context, dipoles);
      }
      mh$.invokeExact(target, context, dipoles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getInducedDipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getInducedDipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getInducedDipoles$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getInducedDipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getInducedDipoles$handle() {
    return OpenMM_AmoebaMultipoleForce_getInducedDipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getInducedDipoles$address() {
    return OpenMM_AmoebaMultipoleForce_getInducedDipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getInducedDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getInducedDipoles(MemorySegment target, MemorySegment context, MemorySegment dipoles) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getInducedDipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getInducedDipoles", target, context, dipoles);
      }
      mh$.invokeExact(target, context, dipoles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getTotalDipoles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getTotalDipoles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getTotalDipoles$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getTotalDipoles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getTotalDipoles$handle() {
    return OpenMM_AmoebaMultipoleForce_getTotalDipoles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getTotalDipoles$address() {
    return OpenMM_AmoebaMultipoleForce_getTotalDipoles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getTotalDipoles(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_Vec3Array *dipoles)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getTotalDipoles(MemorySegment target, MemorySegment context, MemorySegment dipoles) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getTotalDipoles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getTotalDipoles", target, context, dipoles);
      }
      mh$.invokeExact(target, context, dipoles);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getElectrostaticPotential {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getElectrostaticPotential");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce *target, const OpenMM_Vec3Array *inputGrid, OpenMM_Context *context, OpenMM_DoubleArray *outputElectrostaticPotential)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getElectrostaticPotential$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getElectrostaticPotential.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce *target, const OpenMM_Vec3Array *inputGrid, OpenMM_Context *context, OpenMM_DoubleArray *outputElectrostaticPotential)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getElectrostaticPotential$handle() {
    return OpenMM_AmoebaMultipoleForce_getElectrostaticPotential.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce *target, const OpenMM_Vec3Array *inputGrid, OpenMM_Context *context, OpenMM_DoubleArray *outputElectrostaticPotential)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getElectrostaticPotential$address() {
    return OpenMM_AmoebaMultipoleForce_getElectrostaticPotential.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(OpenMM_AmoebaMultipoleForce *target, const OpenMM_Vec3Array *inputGrid, OpenMM_Context *context, OpenMM_DoubleArray *outputElectrostaticPotential)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getElectrostaticPotential(MemorySegment target, MemorySegment inputGrid, MemorySegment context, MemorySegment outputElectrostaticPotential) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getElectrostaticPotential.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getElectrostaticPotential", target, inputGrid, context, outputElectrostaticPotential);
      }
      mh$.invokeExact(target, inputGrid, context, outputElectrostaticPotential);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_DoubleArray *outputMultipoleMoments)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments$descriptor() {
    return OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_DoubleArray *outputMultipoleMoments)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments$handle() {
    return OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_DoubleArray *outputMultipoleMoments)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments$address() {
    return OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context, OpenMM_DoubleArray *outputMultipoleMoments)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments(MemorySegment target, MemorySegment context, MemorySegment outputMultipoleMoments) {
    var mh$ = OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_getSystemMultipoleMoments", target, context, outputMultipoleMoments);
      }
      mh$.invokeExact(target, context, outputMultipoleMoments);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_updateParametersInContext$descriptor() {
    return OpenMM_AmoebaMultipoleForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_updateParametersInContext$handle() {
    return OpenMM_AmoebaMultipoleForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_updateParametersInContext$address() {
    return OpenMM_AmoebaMultipoleForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaMultipoleForce_updateParametersInContext(OpenMM_AmoebaMultipoleForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_AmoebaMultipoleForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_AmoebaMultipoleForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaMultipoleForce *target)
   *}
   */
  public static int OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaMultipoleForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_AmoebaTorsionTorsionForce *OpenMM_AmoebaTorsionTorsionForce_create()
   *}
   */
  public static class OpenMM_AmoebaTorsionTorsionForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_AmoebaTorsionTorsionForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_AmoebaTorsionTorsionForce *OpenMM_AmoebaTorsionTorsionForce_create()
     *}
     */
    public static OpenMM_AmoebaTorsionTorsionForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_AmoebaTorsionTorsionForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_AmoebaTorsionTorsionForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_destroy$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_destroy$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_destroy$address() {
    return OpenMM_AmoebaTorsionTorsionForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_destroy(OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static void OpenMM_AmoebaTorsionTorsionForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions$address() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids$address() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static int OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(MemorySegment target) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce *target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce *target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce *target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion$address() {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(OpenMM_AmoebaTorsionTorsionForce *target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static int OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(MemorySegment target, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion", target, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
      }
      return (int) mh$.invokeExact(target, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce *target, int index, int *particle1, int *particle2, int *particle3, int *particle4, int *particle5, int *chiralCheckAtomIndex, int *gridIndex)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce *target, int index, int *particle1, int *particle2, int *particle3, int *particle4, int *particle5, int *chiralCheckAtomIndex, int *gridIndex)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce *target, int index, int *particle1, int *particle2, int *particle3, int *particle4, int *particle5, int *chiralCheckAtomIndex, int *gridIndex)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters$address() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(const OpenMM_AmoebaTorsionTorsionForce *target, int index, int *particle1, int *particle2, int *particle3, int *particle4, int *particle5, int *chiralCheckAtomIndex, int *gridIndex)
   *}
   */
  public static void OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(MemorySegment target, int index, MemorySegment particle1, MemorySegment particle2, MemorySegment particle3, MemorySegment particle4, MemorySegment particle5, MemorySegment chiralCheckAtomIndex, MemorySegment gridIndex) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters", target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
      }
      mh$.invokeExact(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce *target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce *target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce *target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters$address() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(OpenMM_AmoebaTorsionTorsionForce *target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex)
   *}
   */
  public static void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(MemorySegment target, int index, int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters", target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
      }
      mh$.invokeExact(target, index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern const OpenMM_3D_DoubleArray *OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce *target, int index)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern const OpenMM_3D_DoubleArray *OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce *target, int index)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern const OpenMM_3D_DoubleArray *OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid$address() {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern const OpenMM_3D_DoubleArray *OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(const OpenMM_AmoebaTorsionTorsionForce *target, int index)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(MemorySegment target, int index) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid", target, index);
      }
      return (MemorySegment) mh$.invokeExact(target, index);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce *target, int index, const OpenMM_3D_DoubleArray *grid)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce *target, int index, const OpenMM_3D_DoubleArray *grid)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce *target, int index, const OpenMM_3D_DoubleArray *grid)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid$address() {
    return OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(OpenMM_AmoebaTorsionTorsionForce *target, int index, const OpenMM_3D_DoubleArray *grid)
   *}
   */
  public static void OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(MemorySegment target, int index, MemorySegment grid) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid", target, index, grid);
      }
      mh$.invokeExact(target, index, grid);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(OpenMM_AmoebaTorsionTorsionForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static void OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(MemorySegment target, int periodic) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions", target, periodic);
      }
      mh$.invokeExact(target, periodic);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MethodHandle OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static MemorySegment OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(const OpenMM_AmoebaTorsionTorsionForce *target)
   *}
   */
  public static int OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  /**
   * Variadic invoker class for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeForce *OpenMM_DrudeForce_create()
   *}
   */
  public static class OpenMM_DrudeForce_create {
    private static final FunctionDescriptor BASE_DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER);
    private static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_create");

    private final MethodHandle handle;
    private final FunctionDescriptor descriptor;
    private final MethodHandle spreader;

    private OpenMM_DrudeForce_create(MethodHandle handle, FunctionDescriptor descriptor, MethodHandle spreader) {
      this.handle = handle;
      this.descriptor = descriptor;
      this.spreader = spreader;
    }

    /**
     * Variadic invoker factory for:
     * {@snippet lang = c:
     * extern OpenMM_DrudeForce *OpenMM_DrudeForce_create()
     *}
     */
    public static OpenMM_DrudeForce_create makeInvoker(MemoryLayout... layouts) {
      FunctionDescriptor desc$ = BASE_DESC.appendArgumentLayouts(layouts);
      Linker.Option fva$ = Linker.Option.firstVariadicArg(BASE_DESC.argumentLayouts().size());
      var mh$ = Linker.nativeLinker().downcallHandle(ADDR, desc$, fva$);
      var spreader$ = mh$.asSpreader(Object[].class, layouts.length);
      return new OpenMM_DrudeForce_create(mh$, desc$, spreader$);
    }

    /**
     * {@return the address}
     */
    public static MemorySegment address() {
      return ADDR;
    }

    /**
     * {@return the specialized method handle}
     */
    public MethodHandle handle() {
      return handle;
    }

    /**
     * {@return the specialized descriptor}
     */
    public FunctionDescriptor descriptor() {
      return descriptor;
    }

    public MemorySegment apply(Object... x0) {
      try {
        if (TRACE_DOWNCALLS) {
          traceDowncall("OpenMM_DrudeForce_create", x0);
        }
        return (MemorySegment) spreader.invokeExact(x0);
      } catch (IllegalArgumentException | ClassCastException ex$) {
        throw ex$; // rethrow IAE from passing wrong number/type of args
      } catch (Throwable ex$) {
        throw new AssertionError("should not reach here", ex$);
      }
    }
  }

  private static class OpenMM_DrudeForce_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_destroy(OpenMM_DrudeForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_destroy$descriptor() {
    return OpenMM_DrudeForce_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_destroy(OpenMM_DrudeForce *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_destroy$handle() {
    return OpenMM_DrudeForce_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_destroy(OpenMM_DrudeForce *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_destroy$address() {
    return OpenMM_DrudeForce_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_destroy(OpenMM_DrudeForce *target)
   *}
   */
  public static void OpenMM_DrudeForce_destroy(MemorySegment target) {
    var mh$ = OpenMM_DrudeForce_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_getNumParticles {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_getNumParticles");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumParticles(const OpenMM_DrudeForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_getNumParticles$descriptor() {
    return OpenMM_DrudeForce_getNumParticles.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumParticles(const OpenMM_DrudeForce *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_getNumParticles$handle() {
    return OpenMM_DrudeForce_getNumParticles.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumParticles(const OpenMM_DrudeForce *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_getNumParticles$address() {
    return OpenMM_DrudeForce_getNumParticles.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumParticles(const OpenMM_DrudeForce *target)
   *}
   */
  public static int OpenMM_DrudeForce_getNumParticles(MemorySegment target) {
    var mh$ = OpenMM_DrudeForce_getNumParticles.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_getNumParticles", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_getNumScreenedPairs {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_getNumScreenedPairs");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumScreenedPairs(const OpenMM_DrudeForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_getNumScreenedPairs$descriptor() {
    return OpenMM_DrudeForce_getNumScreenedPairs.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumScreenedPairs(const OpenMM_DrudeForce *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_getNumScreenedPairs$handle() {
    return OpenMM_DrudeForce_getNumScreenedPairs.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumScreenedPairs(const OpenMM_DrudeForce *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_getNumScreenedPairs$address() {
    return OpenMM_DrudeForce_getNumScreenedPairs.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_getNumScreenedPairs(const OpenMM_DrudeForce *target)
   *}
   */
  public static int OpenMM_DrudeForce_getNumScreenedPairs(MemorySegment target) {
    var mh$ = OpenMM_DrudeForce_getNumScreenedPairs.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_getNumScreenedPairs", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_addParticle {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_addParticle");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addParticle(OpenMM_DrudeForce *target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_addParticle$descriptor() {
    return OpenMM_DrudeForce_addParticle.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addParticle(OpenMM_DrudeForce *target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_addParticle$handle() {
    return OpenMM_DrudeForce_addParticle.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addParticle(OpenMM_DrudeForce *target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_addParticle$address() {
    return OpenMM_DrudeForce_addParticle.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addParticle(OpenMM_DrudeForce *target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static int OpenMM_DrudeForce_addParticle(MemorySegment target, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34) {
    var mh$ = OpenMM_DrudeForce_addParticle.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_addParticle", target, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
      }
      return (int) mh$.invokeExact(target, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_getParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_getParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getParticleParameters(const OpenMM_DrudeForce *target, int index, int *particle, int *particle1, int *particle2, int *particle3, int *particle4, double *charge, double *polarizability, double *aniso12, double *aniso34)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_getParticleParameters$descriptor() {
    return OpenMM_DrudeForce_getParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getParticleParameters(const OpenMM_DrudeForce *target, int index, int *particle, int *particle1, int *particle2, int *particle3, int *particle4, double *charge, double *polarizability, double *aniso12, double *aniso34)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_getParticleParameters$handle() {
    return OpenMM_DrudeForce_getParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getParticleParameters(const OpenMM_DrudeForce *target, int index, int *particle, int *particle1, int *particle2, int *particle3, int *particle4, double *charge, double *polarizability, double *aniso12, double *aniso34)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_getParticleParameters$address() {
    return OpenMM_DrudeForce_getParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getParticleParameters(const OpenMM_DrudeForce *target, int index, int *particle, int *particle1, int *particle2, int *particle3, int *particle4, double *charge, double *polarizability, double *aniso12, double *aniso34)
   *}
   */
  public static void OpenMM_DrudeForce_getParticleParameters(MemorySegment target, int index, MemorySegment particle, MemorySegment particle1, MemorySegment particle2, MemorySegment particle3, MemorySegment particle4, MemorySegment charge, MemorySegment polarizability, MemorySegment aniso12, MemorySegment aniso34) {
    var mh$ = OpenMM_DrudeForce_getParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_getParticleParameters", target, index, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
      }
      mh$.invokeExact(target, index, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_setParticleParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_setParticleParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setParticleParameters(OpenMM_DrudeForce *target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_setParticleParameters$descriptor() {
    return OpenMM_DrudeForce_setParticleParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setParticleParameters(OpenMM_DrudeForce *target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_setParticleParameters$handle() {
    return OpenMM_DrudeForce_setParticleParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setParticleParameters(OpenMM_DrudeForce *target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_setParticleParameters$address() {
    return OpenMM_DrudeForce_setParticleParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setParticleParameters(OpenMM_DrudeForce *target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34)
   *}
   */
  public static void OpenMM_DrudeForce_setParticleParameters(MemorySegment target, int index, int particle, int particle1, int particle2, int particle3, int particle4, double charge, double polarizability, double aniso12, double aniso34) {
    var mh$ = OpenMM_DrudeForce_setParticleParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_setParticleParameters", target, index, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
      }
      mh$.invokeExact(target, index, particle, particle1, particle2, particle3, particle4, charge, polarizability, aniso12, aniso34);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_addScreenedPair {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_addScreenedPair");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addScreenedPair(OpenMM_DrudeForce *target, int particle1, int particle2, double thole)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_addScreenedPair$descriptor() {
    return OpenMM_DrudeForce_addScreenedPair.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addScreenedPair(OpenMM_DrudeForce *target, int particle1, int particle2, double thole)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_addScreenedPair$handle() {
    return OpenMM_DrudeForce_addScreenedPair.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addScreenedPair(OpenMM_DrudeForce *target, int particle1, int particle2, double thole)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_addScreenedPair$address() {
    return OpenMM_DrudeForce_addScreenedPair.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_DrudeForce_addScreenedPair(OpenMM_DrudeForce *target, int particle1, int particle2, double thole)
   *}
   */
  public static int OpenMM_DrudeForce_addScreenedPair(MemorySegment target, int particle1, int particle2, double thole) {
    var mh$ = OpenMM_DrudeForce_addScreenedPair.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_addScreenedPair", target, particle1, particle2, thole);
      }
      return (int) mh$.invokeExact(target, particle1, particle2, thole);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_getScreenedPairParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_getScreenedPairParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getScreenedPairParameters(const OpenMM_DrudeForce *target, int index, int *particle1, int *particle2, double *thole)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_getScreenedPairParameters$descriptor() {
    return OpenMM_DrudeForce_getScreenedPairParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getScreenedPairParameters(const OpenMM_DrudeForce *target, int index, int *particle1, int *particle2, double *thole)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_getScreenedPairParameters$handle() {
    return OpenMM_DrudeForce_getScreenedPairParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getScreenedPairParameters(const OpenMM_DrudeForce *target, int index, int *particle1, int *particle2, double *thole)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_getScreenedPairParameters$address() {
    return OpenMM_DrudeForce_getScreenedPairParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_getScreenedPairParameters(const OpenMM_DrudeForce *target, int index, int *particle1, int *particle2, double *thole)
   *}
   */
  public static void OpenMM_DrudeForce_getScreenedPairParameters(MemorySegment target, int index, MemorySegment particle1, MemorySegment particle2, MemorySegment thole) {
    var mh$ = OpenMM_DrudeForce_getScreenedPairParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_getScreenedPairParameters", target, index, particle1, particle2, thole);
      }
      mh$.invokeExact(target, index, particle1, particle2, thole);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_setScreenedPairParameters {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_setScreenedPairParameters");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setScreenedPairParameters(OpenMM_DrudeForce *target, int index, int particle1, int particle2, double thole)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_setScreenedPairParameters$descriptor() {
    return OpenMM_DrudeForce_setScreenedPairParameters.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setScreenedPairParameters(OpenMM_DrudeForce *target, int index, int particle1, int particle2, double thole)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_setScreenedPairParameters$handle() {
    return OpenMM_DrudeForce_setScreenedPairParameters.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setScreenedPairParameters(OpenMM_DrudeForce *target, int index, int particle1, int particle2, double thole)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_setScreenedPairParameters$address() {
    return OpenMM_DrudeForce_setScreenedPairParameters.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setScreenedPairParameters(OpenMM_DrudeForce *target, int index, int particle1, int particle2, double thole)
   *}
   */
  public static void OpenMM_DrudeForce_setScreenedPairParameters(MemorySegment target, int index, int particle1, int particle2, double thole) {
    var mh$ = OpenMM_DrudeForce_setScreenedPairParameters.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_setScreenedPairParameters", target, index, particle1, particle2, thole);
      }
      mh$.invokeExact(target, index, particle1, particle2, thole);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_updateParametersInContext {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_updateParametersInContext");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_updateParametersInContext(OpenMM_DrudeForce *target, OpenMM_Context *context)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_updateParametersInContext$descriptor() {
    return OpenMM_DrudeForce_updateParametersInContext.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_updateParametersInContext(OpenMM_DrudeForce *target, OpenMM_Context *context)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_updateParametersInContext$handle() {
    return OpenMM_DrudeForce_updateParametersInContext.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_updateParametersInContext(OpenMM_DrudeForce *target, OpenMM_Context *context)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_updateParametersInContext$address() {
    return OpenMM_DrudeForce_updateParametersInContext.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_updateParametersInContext(OpenMM_DrudeForce *target, OpenMM_Context *context)
   *}
   */
  public static void OpenMM_DrudeForce_updateParametersInContext(MemorySegment target, MemorySegment context) {
    var mh$ = OpenMM_DrudeForce_updateParametersInContext.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_updateParametersInContext", target, context);
      }
      mh$.invokeExact(target, context);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(OpenMM_DrudeForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(OpenMM_DrudeForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions$handle() {
    return OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(OpenMM_DrudeForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions$address() {
    return OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(OpenMM_DrudeForce *target, OpenMM_Boolean periodic)
   *}
   */
  public static void OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions(MemorySegment target, int periodic) {
    var mh$ = OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_setUsesPeriodicBoundaryConditions", target, periodic);
      }
      mh$.invokeExact(target, periodic);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeForce_usesPeriodicBoundaryConditions {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeForce_usesPeriodicBoundaryConditions");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_DrudeForce_usesPeriodicBoundaryConditions(const OpenMM_DrudeForce *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeForce_usesPeriodicBoundaryConditions$descriptor() {
    return OpenMM_DrudeForce_usesPeriodicBoundaryConditions.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_DrudeForce_usesPeriodicBoundaryConditions(const OpenMM_DrudeForce *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeForce_usesPeriodicBoundaryConditions$handle() {
    return OpenMM_DrudeForce_usesPeriodicBoundaryConditions.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_DrudeForce_usesPeriodicBoundaryConditions(const OpenMM_DrudeForce *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeForce_usesPeriodicBoundaryConditions$address() {
    return OpenMM_DrudeForce_usesPeriodicBoundaryConditions.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_Boolean OpenMM_DrudeForce_usesPeriodicBoundaryConditions(const OpenMM_DrudeForce *target)
   *}
   */
  public static int OpenMM_DrudeForce_usesPeriodicBoundaryConditions(MemorySegment target) {
    var mh$ = OpenMM_DrudeForce_usesPeriodicBoundaryConditions.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeForce_usesPeriodicBoundaryConditions", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeIntegrator *OpenMM_DrudeIntegrator_create(double stepSize)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_create$descriptor() {
    return OpenMM_DrudeIntegrator_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeIntegrator *OpenMM_DrudeIntegrator_create(double stepSize)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_create$handle() {
    return OpenMM_DrudeIntegrator_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeIntegrator *OpenMM_DrudeIntegrator_create(double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_create$address() {
    return OpenMM_DrudeIntegrator_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_DrudeIntegrator *OpenMM_DrudeIntegrator_create(double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_create(double stepSize) {
    var mh$ = OpenMM_DrudeIntegrator_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_create", stepSize);
      }
      return (MemorySegment) mh$.invokeExact(stepSize);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_destroy(OpenMM_DrudeIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_destroy$descriptor() {
    return OpenMM_DrudeIntegrator_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_destroy(OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_destroy$handle() {
    return OpenMM_DrudeIntegrator_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_destroy(OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_destroy$address() {
    return OpenMM_DrudeIntegrator_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_destroy(OpenMM_DrudeIntegrator *target)
   *}
   */
  public static void OpenMM_DrudeIntegrator_destroy(MemorySegment target) {
    var mh$ = OpenMM_DrudeIntegrator_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_step {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_step");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_step(OpenMM_DrudeIntegrator *target, int steps)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_step$descriptor() {
    return OpenMM_DrudeIntegrator_step.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_step(OpenMM_DrudeIntegrator *target, int steps)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_step$handle() {
    return OpenMM_DrudeIntegrator_step.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_step(OpenMM_DrudeIntegrator *target, int steps)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_step$address() {
    return OpenMM_DrudeIntegrator_step.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_step(OpenMM_DrudeIntegrator *target, int steps)
   *}
   */
  public static void OpenMM_DrudeIntegrator_step(MemorySegment target, int steps) {
    var mh$ = OpenMM_DrudeIntegrator_step.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_step", target, steps);
      }
      mh$.invokeExact(target, steps);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_getDrudeTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_getDrudeTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getDrudeTemperature(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_getDrudeTemperature$descriptor() {
    return OpenMM_DrudeIntegrator_getDrudeTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getDrudeTemperature(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_getDrudeTemperature$handle() {
    return OpenMM_DrudeIntegrator_getDrudeTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getDrudeTemperature(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_getDrudeTemperature$address() {
    return OpenMM_DrudeIntegrator_getDrudeTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getDrudeTemperature(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeIntegrator_getDrudeTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeIntegrator_getDrudeTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_getDrudeTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_setDrudeTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_setDrudeTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setDrudeTemperature(OpenMM_DrudeIntegrator *target, double temp)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_setDrudeTemperature$descriptor() {
    return OpenMM_DrudeIntegrator_setDrudeTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setDrudeTemperature(OpenMM_DrudeIntegrator *target, double temp)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_setDrudeTemperature$handle() {
    return OpenMM_DrudeIntegrator_setDrudeTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setDrudeTemperature(OpenMM_DrudeIntegrator *target, double temp)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_setDrudeTemperature$address() {
    return OpenMM_DrudeIntegrator_setDrudeTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setDrudeTemperature(OpenMM_DrudeIntegrator *target, double temp)
   *}
   */
  public static void OpenMM_DrudeIntegrator_setDrudeTemperature(MemorySegment target, double temp) {
    var mh$ = OpenMM_DrudeIntegrator_setDrudeTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_setDrudeTemperature", target, temp);
      }
      mh$.invokeExact(target, temp);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_getMaxDrudeDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_getMaxDrudeDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getMaxDrudeDistance(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_getMaxDrudeDistance$descriptor() {
    return OpenMM_DrudeIntegrator_getMaxDrudeDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getMaxDrudeDistance(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_getMaxDrudeDistance$handle() {
    return OpenMM_DrudeIntegrator_getMaxDrudeDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getMaxDrudeDistance(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_getMaxDrudeDistance$address() {
    return OpenMM_DrudeIntegrator_getMaxDrudeDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeIntegrator_getMaxDrudeDistance(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeIntegrator_getMaxDrudeDistance(MemorySegment target) {
    var mh$ = OpenMM_DrudeIntegrator_getMaxDrudeDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_getMaxDrudeDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_setMaxDrudeDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_setMaxDrudeDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setMaxDrudeDistance(OpenMM_DrudeIntegrator *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_setMaxDrudeDistance$descriptor() {
    return OpenMM_DrudeIntegrator_setMaxDrudeDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setMaxDrudeDistance(OpenMM_DrudeIntegrator *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_setMaxDrudeDistance$handle() {
    return OpenMM_DrudeIntegrator_setMaxDrudeDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setMaxDrudeDistance(OpenMM_DrudeIntegrator *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_setMaxDrudeDistance$address() {
    return OpenMM_DrudeIntegrator_setMaxDrudeDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setMaxDrudeDistance(OpenMM_DrudeIntegrator *target, double distance)
   *}
   */
  public static void OpenMM_DrudeIntegrator_setMaxDrudeDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_DrudeIntegrator_setMaxDrudeDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_setMaxDrudeDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_setRandomNumberSeed {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_setRandomNumberSeed");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setRandomNumberSeed(OpenMM_DrudeIntegrator *target, int seed)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_setRandomNumberSeed$descriptor() {
    return OpenMM_DrudeIntegrator_setRandomNumberSeed.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setRandomNumberSeed(OpenMM_DrudeIntegrator *target, int seed)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_setRandomNumberSeed$handle() {
    return OpenMM_DrudeIntegrator_setRandomNumberSeed.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setRandomNumberSeed(OpenMM_DrudeIntegrator *target, int seed)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_setRandomNumberSeed$address() {
    return OpenMM_DrudeIntegrator_setRandomNumberSeed.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeIntegrator_setRandomNumberSeed(OpenMM_DrudeIntegrator *target, int seed)
   *}
   */
  public static void OpenMM_DrudeIntegrator_setRandomNumberSeed(MemorySegment target, int seed) {
    var mh$ = OpenMM_DrudeIntegrator_setRandomNumberSeed.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_setRandomNumberSeed", target, seed);
      }
      mh$.invokeExact(target, seed);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeIntegrator_getRandomNumberSeed {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_INT,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeIntegrator_getRandomNumberSeed");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeIntegrator_getRandomNumberSeed(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeIntegrator_getRandomNumberSeed$descriptor() {
    return OpenMM_DrudeIntegrator_getRandomNumberSeed.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeIntegrator_getRandomNumberSeed(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeIntegrator_getRandomNumberSeed$handle() {
    return OpenMM_DrudeIntegrator_getRandomNumberSeed.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern int OpenMM_DrudeIntegrator_getRandomNumberSeed(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeIntegrator_getRandomNumberSeed$address() {
    return OpenMM_DrudeIntegrator_getRandomNumberSeed.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern int OpenMM_DrudeIntegrator_getRandomNumberSeed(const OpenMM_DrudeIntegrator *target)
   *}
   */
  public static int OpenMM_DrudeIntegrator_getRandomNumberSeed(MemorySegment target) {
    var mh$ = OpenMM_DrudeIntegrator_getRandomNumberSeed.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeIntegrator_getRandomNumberSeed", target);
      }
      return (int) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeNoseHooverIntegrator *OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_create$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeNoseHooverIntegrator *OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_create$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeNoseHooverIntegrator *OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_create$address() {
    return OpenMM_DrudeNoseHooverIntegrator_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_DrudeNoseHooverIntegrator *OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_create(double temperature, double collisionFrequency, double drudeTemperature, double drudeCollisionFrequency, double stepSize, int chainLength, int numMTS, int numYoshidaSuzuki) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_create", temperature, collisionFrequency, drudeTemperature, drudeCollisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki);
      }
      return (MemorySegment) mh$.invokeExact(temperature, collisionFrequency, drudeTemperature, drudeCollisionFrequency, stepSize, chainLength, numMTS, numYoshidaSuzuki);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_destroy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_destroy$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_destroy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_destroy$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_destroy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_destroy$address() {
    return OpenMM_DrudeNoseHooverIntegrator_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_destroy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static void OpenMM_DrudeNoseHooverIntegrator_destroy(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(const OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(const OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(const OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance$address() {
    return OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(const OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_getMaxDrudeDistance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(OpenMM_DrudeNoseHooverIntegrator *target, double distance)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(OpenMM_DrudeNoseHooverIntegrator *target, double distance)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(OpenMM_DrudeNoseHooverIntegrator *target, double distance)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance$address() {
    return OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(OpenMM_DrudeNoseHooverIntegrator *target, double distance)
   *}
   */
  public static void OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance(MemorySegment target, double distance) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_setMaxDrudeDistance", target, distance);
      }
      mh$.invokeExact(target, distance);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy$address() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_computeDrudeKineticEnergy", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy$address() {
    return OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_computeTotalKineticEnergy", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature$address() {
    return OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_computeSystemTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature$descriptor() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature$handle() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature$address() {
    return OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(OpenMM_DrudeNoseHooverIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeNoseHooverIntegrator_computeDrudeTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeSCFIntegrator_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeSCFIntegrator_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeSCFIntegrator *OpenMM_DrudeSCFIntegrator_create(double stepSize)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeSCFIntegrator_create$descriptor() {
    return OpenMM_DrudeSCFIntegrator_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeSCFIntegrator *OpenMM_DrudeSCFIntegrator_create(double stepSize)
   *}
   */
  public static MethodHandle OpenMM_DrudeSCFIntegrator_create$handle() {
    return OpenMM_DrudeSCFIntegrator_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeSCFIntegrator *OpenMM_DrudeSCFIntegrator_create(double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_create$address() {
    return OpenMM_DrudeSCFIntegrator_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_DrudeSCFIntegrator *OpenMM_DrudeSCFIntegrator_create(double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_create(double stepSize) {
    var mh$ = OpenMM_DrudeSCFIntegrator_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeSCFIntegrator_create", stepSize);
      }
      return (MemorySegment) mh$.invokeExact(stepSize);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeSCFIntegrator_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeSCFIntegrator_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_destroy(OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeSCFIntegrator_destroy$descriptor() {
    return OpenMM_DrudeSCFIntegrator_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_destroy(OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeSCFIntegrator_destroy$handle() {
    return OpenMM_DrudeSCFIntegrator_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_destroy(OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_destroy$address() {
    return OpenMM_DrudeSCFIntegrator_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_destroy(OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static void OpenMM_DrudeSCFIntegrator_destroy(MemorySegment target) {
    var mh$ = OpenMM_DrudeSCFIntegrator_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeSCFIntegrator_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(const OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance$descriptor() {
    return OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(const OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance$handle() {
    return OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(const OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance$address() {
    return OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(const OpenMM_DrudeSCFIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(MemorySegment target) {
    var mh$ = OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(OpenMM_DrudeSCFIntegrator *target, double tol)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance$descriptor() {
    return OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(OpenMM_DrudeSCFIntegrator *target, double tol)
   *}
   */
  public static MethodHandle OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance$handle() {
    return OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(OpenMM_DrudeSCFIntegrator *target, double tol)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance$address() {
    return OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(OpenMM_DrudeSCFIntegrator *target, double tol)
   *}
   */
  public static void OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(MemorySegment target, double tol) {
    var mh$ = OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance", target, tol);
      }
      mh$.invokeExact(target, tol);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeSCFIntegrator_step {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeSCFIntegrator_step");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_step(OpenMM_DrudeSCFIntegrator *target, int steps)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeSCFIntegrator_step$descriptor() {
    return OpenMM_DrudeSCFIntegrator_step.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_step(OpenMM_DrudeSCFIntegrator *target, int steps)
   *}
   */
  public static MethodHandle OpenMM_DrudeSCFIntegrator_step$handle() {
    return OpenMM_DrudeSCFIntegrator_step.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_step(OpenMM_DrudeSCFIntegrator *target, int steps)
   *}
   */
  public static MemorySegment OpenMM_DrudeSCFIntegrator_step$address() {
    return OpenMM_DrudeSCFIntegrator_step.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeSCFIntegrator_step(OpenMM_DrudeSCFIntegrator *target, int steps)
   *}
   */
  public static void OpenMM_DrudeSCFIntegrator_step(MemorySegment target, int steps) {
    var mh$ = OpenMM_DrudeSCFIntegrator_step.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeSCFIntegrator_step", target, steps);
      }
      mh$.invokeExact(target, steps);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_create {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_create");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeLangevinIntegrator *OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_create$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_create.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeLangevinIntegrator *OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_create$handle() {
    return OpenMM_DrudeLangevinIntegrator_create.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern OpenMM_DrudeLangevinIntegrator *OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_create$address() {
    return OpenMM_DrudeLangevinIntegrator_create.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern OpenMM_DrudeLangevinIntegrator *OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_create(double temperature, double frictionCoeff, double drudeTemperature, double drudeFrictionCoeff, double stepSize) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_create.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_create", temperature, frictionCoeff, drudeTemperature, drudeFrictionCoeff, stepSize);
      }
      return (MemorySegment) mh$.invokeExact(temperature, frictionCoeff, drudeTemperature, drudeFrictionCoeff, stepSize);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_destroy {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_destroy");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_destroy(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_destroy$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_destroy.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_destroy(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_destroy$handle() {
    return OpenMM_DrudeLangevinIntegrator_destroy.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_destroy(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_destroy$address() {
    return OpenMM_DrudeLangevinIntegrator_destroy.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_destroy(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static void OpenMM_DrudeLangevinIntegrator_destroy(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_destroy.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_destroy", target);
      }
      mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_getTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_getTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getTemperature(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_getTemperature$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_getTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getTemperature(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_getTemperature$handle() {
    return OpenMM_DrudeLangevinIntegrator_getTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getTemperature(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_getTemperature$address() {
    return OpenMM_DrudeLangevinIntegrator_getTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getTemperature(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeLangevinIntegrator_getTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_getTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_getTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_setTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_setTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setTemperature(OpenMM_DrudeLangevinIntegrator *target, double temp)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_setTemperature$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_setTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setTemperature(OpenMM_DrudeLangevinIntegrator *target, double temp)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_setTemperature$handle() {
    return OpenMM_DrudeLangevinIntegrator_setTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setTemperature(OpenMM_DrudeLangevinIntegrator *target, double temp)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_setTemperature$address() {
    return OpenMM_DrudeLangevinIntegrator_setTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setTemperature(OpenMM_DrudeLangevinIntegrator *target, double temp)
   *}
   */
  public static void OpenMM_DrudeLangevinIntegrator_setTemperature(MemorySegment target, double temp) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_setTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_setTemperature", target, temp);
      }
      mh$.invokeExact(target, temp);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_getFriction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_getFriction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_getFriction$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_getFriction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_getFriction$handle() {
    return OpenMM_DrudeLangevinIntegrator_getFriction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_getFriction$address() {
    return OpenMM_DrudeLangevinIntegrator_getFriction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeLangevinIntegrator_getFriction(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_getFriction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_getFriction", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_setFriction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_setFriction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_setFriction$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_setFriction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_setFriction$handle() {
    return OpenMM_DrudeLangevinIntegrator_setFriction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_setFriction$address() {
    return OpenMM_DrudeLangevinIntegrator_setFriction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static void OpenMM_DrudeLangevinIntegrator_setFriction(MemorySegment target, double coeff) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_setFriction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_setFriction", target, coeff);
      }
      mh$.invokeExact(target, coeff);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_getDrudeFriction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_getDrudeFriction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_getDrudeFriction$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_getDrudeFriction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_getDrudeFriction$handle() {
    return OpenMM_DrudeLangevinIntegrator_getDrudeFriction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_getDrudeFriction$address() {
    return OpenMM_DrudeLangevinIntegrator_getDrudeFriction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(const OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeLangevinIntegrator_getDrudeFriction(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_getDrudeFriction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_getDrudeFriction", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_setDrudeFriction {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_DOUBLE
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_setDrudeFriction");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_setDrudeFriction$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_setDrudeFriction.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_setDrudeFriction$handle() {
    return OpenMM_DrudeLangevinIntegrator_setDrudeFriction.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_setDrudeFriction$address() {
    return OpenMM_DrudeLangevinIntegrator_setDrudeFriction.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(OpenMM_DrudeLangevinIntegrator *target, double coeff)
   *}
   */
  public static void OpenMM_DrudeLangevinIntegrator_setDrudeFriction(MemorySegment target, double coeff) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_setDrudeFriction.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_setDrudeFriction", target, coeff);
      }
      mh$.invokeExact(target, coeff);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_step {
    public static final FunctionDescriptor DESC = FunctionDescriptor.ofVoid(
        OpenMMNative.C_POINTER,
        OpenMMNative.C_INT
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_step");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_step(OpenMM_DrudeLangevinIntegrator *target, int steps)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_step$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_step.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_step(OpenMM_DrudeLangevinIntegrator *target, int steps)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_step$handle() {
    return OpenMM_DrudeLangevinIntegrator_step.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_step(OpenMM_DrudeLangevinIntegrator *target, int steps)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_step$address() {
    return OpenMM_DrudeLangevinIntegrator_step.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern void OpenMM_DrudeLangevinIntegrator_step(OpenMM_DrudeLangevinIntegrator *target, int steps)
   *}
   */
  public static void OpenMM_DrudeLangevinIntegrator_step(MemorySegment target, int steps) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_step.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_step", target, steps);
      }
      mh$.invokeExact(target, steps);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_computeSystemTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_computeSystemTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_computeSystemTemperature$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_computeSystemTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_computeSystemTemperature$handle() {
    return OpenMM_DrudeLangevinIntegrator_computeSystemTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_computeSystemTemperature$address() {
    return OpenMM_DrudeLangevinIntegrator_computeSystemTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_computeSystemTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_computeSystemTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }

  private static class OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature {
    public static final FunctionDescriptor DESC = FunctionDescriptor.of(
        OpenMMNative.C_DOUBLE,
        OpenMMNative.C_POINTER
    );

    public static final MemorySegment ADDR = SYMBOL_LOOKUP.findOrThrow("OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature");

    public static final MethodHandle HANDLE = Linker.nativeLinker().downcallHandle(ADDR, DESC);
  }

  /**
   * Function descriptor for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static FunctionDescriptor OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature$descriptor() {
    return OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature.DESC;
  }

  /**
   * Downcall method handle for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MethodHandle OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature$handle() {
    return OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature.HANDLE;
  }

  /**
   * Address for:
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static MemorySegment OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature$address() {
    return OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature.ADDR;
  }

  /**
   * {@snippet lang = c:
   * extern double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(OpenMM_DrudeLangevinIntegrator *target)
   *}
   */
  public static double OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(MemorySegment target) {
    var mh$ = OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature.HANDLE;
    try {
      if (TRACE_DOWNCALLS) {
        traceDowncall("OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature", target);
      }
      return (double) mh$.invokeExact(target);
    } catch (Error | RuntimeException ex) {
      throw ex;
    } catch (Throwable ex$) {
      throw new AssertionError("should not reach here", ex$);
    }
  }
}

