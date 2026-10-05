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
package ffx.openmm.ffm;

import java.lang.foreign.MemorySegment;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Base for wrappers around opaque native OpenMM handles.
 *
 * <p>This class tracks whether a Java wrapper has a usable native pointer, but does not choose the
 * native ownership policy: subclasses define whether {@link #destroy()} releases the handle or
 * merely invalidates the wrapper. The pointer segment is not copied or kept alive by this class;
 * callers that supply a scoped segment must keep its scope alive and ensure that the address
 * continues to refer to a live native object.</p>
 */
public abstract class OpenMMHandle implements AutoCloseable {

  private MemorySegment pointer;

  /**
   * Wrap a non-null native OpenMM handle.
   *
   * <p>This constructor validates the address but does not itself establish native ownership.</p>
   *
   * @param pointer non-null native handle with a nonzero address.
   * @throws IllegalArgumentException if the segment has address zero.
   * @throws NullPointerException if {@code pointer} is null.
   */
  protected OpenMMHandle(MemorySegment pointer) {
    this(pointer, false);
  }

  /**
   * Create a wrapper that may initially have no native handle when explicitly permitted.
   *
   * @param pointer   non-null native handle, or a null-address segment if {@code allowNull} is true.
   * @param allowNull whether to permit an address-zero segment.
   * @throws IllegalArgumentException if the address is zero and {@code allowNull} is false.
   * @throws NullPointerException if {@code pointer} is null.
   */
  protected OpenMMHandle(MemorySegment pointer, boolean allowNull) {
    this.pointer = Objects.requireNonNull(pointer, "OpenMM handle cannot be null.");
    if (!allowNull && pointer.address() == 0) {
      throw new IllegalArgumentException("OpenMM handle address cannot be null.");
    }
  }

  /**
   * Release or invalidate this wrapper using its subclass-defined destruction policy.
   *
   * @see #destroy()
   */
  @Override
  public final void close() {
    destroy();
  }

  /**
   * Apply this subclass's native-resource release or wrapper-invalidation policy.
   *
   * <p>Implementations should make repeated calls harmless.</p>
   */
  public abstract void destroy();

  /**
   * Get this wrapper's current native handle.
   *
   * @return live native handle.
   * @throws IllegalStateException if this wrapper has been destroyed or invalidated.
   */
  public final MemorySegment getPointer() {
    return requirePointer();
  }

  /**
   * Return whether this wrapper no longer has a usable native handle.
   *
   * @return {@code true} if the wrapper's handle is absent or has address zero.
   */
  public final boolean isDestroyed() {
    return pointer == null || pointer.address() == 0;
  }

  /**
   * Invoke the specified native destructor once, then invalidate this Java handle.
   *
   * @param destructor non-null native destructor to invoke for a present handle.
   */
  protected final void destroy(Consumer<MemorySegment> destructor) {
    MemorySegment current = pointer;
    pointer = null;
    if (current != null && current.address() != 0) {
      destructor.accept(current);
    }
  }

  /**
   * Install a native handle after deferred initialization.
   *
   * @param replacement non-null handle with a nonzero address; this wrapper must be invalidated.
   * @throws IllegalStateException if this wrapper currently has a live handle.
   * @throws IllegalArgumentException if {@code replacement} has address zero.
   * @throws NullPointerException if {@code replacement} is null.
   */
  protected final void replacePointer(MemorySegment replacement) {
    Objects.requireNonNull(replacement, "OpenMM handle cannot be null.");
    if (!isDestroyed()) {
      throw new IllegalStateException("Cannot replace a live OpenMM handle.");
    }
    if (replacement.address() == 0) {
      throw new IllegalArgumentException("Replacement OpenMM handle address cannot be null.");
    }
    pointer = replacement;
  }

  /**
   * Rebind this wrapper without releasing the previously referenced handle.
   *
   * <p>This is provided only for compatibility with legacy façade pointer setters. Callers remain
   * responsible for the displaced handle.</p>
   *
   * @param replacement non-null replacement handle with a nonzero address.
   * @throws IllegalArgumentException if {@code replacement} has address zero.
   * @throws NullPointerException if {@code replacement} is null.
   */
  protected final void rebindPointer(MemorySegment replacement) {
    Objects.requireNonNull(replacement, "OpenMM handle cannot be null.");
    if (replacement.address() == 0) {
      throw new IllegalArgumentException("Replacement OpenMM handle address cannot be null.");
    }
    pointer = replacement;
  }

  /**
   * Invalidate this wrapper after a different native owner releases its handle.
   *
   * <p>This does not invoke a native destructor. It is used when OpenMM transfers ownership of a
   * handle to another object, such as a {@link System}.</p>
   */
  protected final void invalidate() {
    pointer = null;
  }

  /**
   * Get the current handle, requiring a nonzero address.
   *
   * @return live native handle.
   * @throws IllegalStateException if this wrapper has been destroyed or invalidated.
   */
  protected final MemorySegment requirePointer() {
    if (pointer == null || pointer.address() == 0) {
      throw new IllegalStateException(getClass().getSimpleName() + " has been destroyed.");
    }
    return pointer;
  }
}
