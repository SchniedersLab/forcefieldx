# OpenMM JNA-to-FFM migration plan

## Current status

The Java façade implementation phase (Phase 1), potential energy term ports (Phase 2, Phases A–E),
and consumer wiring / dual-topology alchemical potentials (Phase F) are complete:
- FFM counterparts exist for core, AMOEBA, and Drude JNA façades.
- All potential energy terms (bonded interactions, restraint forces, fixed-charge Lennard-Jones/Coulomb/PME,
  continuum GB solvation, alchemical softcore terms, and the AMOEBA polarizable force field with 14-7 vdW,
  polarizable multipoles, GK, and WCA dispersion) have been implemented under `ffx.potential.ommffm`.
- Dynamics engines, integrators (Verlet, Langevin, Custom MTS, Custom MTS Langevin), and state polling
  have been implemented and verified.
- Dual-topology hybrid system setup (`OpenMMDualTopologySystem`) and energy evaluation (`OpenMMDualTopologyEnergy`)
  have been implemented with full force class dual-topology constructors and parameter updaters.
- Backend switch configuration (`ffx.openmm.backend=jna|ffm`) and consumer algorithm wiring (`MolecularDynamicsOpenMM`,
  `MinimizeOpenMM`) have been implemented and verified.
- Phase-specific unit test suites (`OpenMMFFMBondEnergyTest`, `OpenMMFFMFixedChargeNonbondedEnergyTest`,
  `OpenMMFFMAmoebaEnergyTest`, `OpenMMFFMIntegratorDynamicsTest`, `OpenMMFFMDualTopologyEnergyTest`, and
  `OpenMMFFMConsumerTest`) compile and pass cleanly on JDK 25 against the Reference platform.

The current active milestone is **Phase G (Lab Validation & Coexistence Period)**.

## Problem and approach

`ffx-openmm` currently compiles 78 Java façade classes against `jopenmm-fat` and JNA. The
fat artifact both supplies `OpenMMLibrary`/`OpenMMAmoebaLibrary`/`OpenMMDrudeLibrary` and
unpacks/loads platform native libraries. The project already compiles for Java 25, so its
stable Foreign Function & Memory API can replace JNA without a preview flag.

Generate Java FFM bindings from the installed OpenMM 8.3.4 C ABI headers, then retain the
existing `ffx.openmm` façade API while converting its internals from JNA pointers, structures,
and references to `MemorySegment`, generated layouts, and confined `Arena` allocations. Do
not attempt to run `jextract` over `OpenMMAmoeba.h` or `OpenMMDrude.h`: those headers expose
C++ classes. The installed `AmoebaOpenMMCWrapper.h` and `DrudeOpenMMCWrapper.h` are the
corresponding C ABI wrapper headers and are the correct generator inputs.

## Reproducible binding generation

The maintained `modules/openmm/jextract.sh` script resolves the repository root, requires
jextract 25, clears only the known generated binding package, and generates one binding from
the version-controlled aggregate header. Regenerate with:

```bash
JEXTRACT=/path/to/jextract-25/bin/jextract \
OPENMM_INSTALL=/path/to/openmm-install \
  modules/openmm/jextract.sh
```

`OpenMM_FFM.h` includes the core, AMOEBA, and Drude C wrapper headers so base ABI declarations
are available to extensions and shared declarations are generated once. Do not pass
`--library`: this keeps generated sources independent of the generation machine. Generated
Java is written under the ordinary `modules/openmm/src/main/java` Maven source root, so no
additional generated-source-root plugin is needed. At runtime,
`OpenMMRuntime` reads `FFX_OPENMM_LIB_DIR`, uses `System.mapLibraryName` to select the
platform-native extension, and loads core OpenMM, AMOEBA, and Drude in that order. It loads
plugins from `FFX_OPENMM_PLUGIN_DIR` or the `plugins` child of the library directory, staging
a temporary directory that excludes unsupported RPMD plugins. Confirm that generation is
byte-stable (or define the allowed generated-code tool-version metadata) and commit the
generated source rather than requiring normal Maven builds to have `jextract` installed.
Document the exact OpenMM and jextract versions used.

## Parallel JNA and FFM transition

Keep the stable JNA backend intact while the FFM backend is developed and validated. The two
backends must not exchange native handles: a JNA `Pointer` and an FFM `MemorySegment` can
refer to distinct OpenMM library instances, and crossing a `System`, `Context`, `Force`, or
array handle between them is unsafe. Compare backends in separate test processes and select
one backend for every production simulation object graph.

Use this source layout:

```text
ffx.openmm                    Existing JNA public façade; default during transition
ffx.openmm.amoeba             Existing JNA AMOEBA façade
ffx.openmm.drude              Existing JNA Drude façade
ffx.openmm.ffm                FFM core façade, runtime, and shared helpers
ffx.openmm.ffm.amoeba         FFM AMOEBA façade
ffx.openmm.ffm.drude          FFM Drude façade
ffx.openmm.ffm.bindings       jextract-generated code only
```

The parallel FFM packages use the same simple class names as their JNA counterparts because
their packages differentiate the backends: for example,
`ffx.openmm.ffm.Context` versus `ffx.openmm.Context`. Put reusable `MemorySegment` ownership,
C-string, `Vec3`,
boolean, and primitive out-parameter helpers directly in `ffx.openmm.ffm`, never in the
generated package.

Expose backend selection in the first downstream integration point using a property such as
`ffx.openmm.backend=jna|ffm`, defaulting to `jna`. Make an FFM selection fail early when
`FFX_OPENMM_LIB_DIR` is absent. Do not alter existing `ffx.potential.openmm` or algorithm
production paths until an all-FFM `System`/`Context` graph and parity tests are available.

### Porting sequence and parity gates

1. **FFM foundation and containers (complete):** `OpenMMRuntime`, `OpenMMHandle`, UTF-8
   string and boolean helpers, `Vec3`, and FFM `BondArray`, `DoubleArray`, `IntArray`,
   `IntSet`, `StringArray`, and `Vec3Array` now coexist with their JNA counterparts.
   `OpenMMFFMContainerTest` validates their native lifecycle, scalar/struct values, UTF-8
   strings, and output integers. `ffx.openmm.ffm.BondArray` is the reference pattern: it uses
   an explicit native handle, idempotent destroy, `AutoCloseable`, and confined `Arena`
   allocations for output integers.

   A follow-up public API audit against the matching JNA container classes found no missing
   operations. `getPointer()` is inherited from `OpenMMHandle`; pointer-taking constructors
   accept `MemorySegment`; and the JNA `BondArray.get(...)` out-parameter overloads are
   represented by the typed `BondArray.Bond` return value. The `Vec3Array` methods use the FFM
   `Vec3` value in place of JNA's `OpenMM_Vec3.ByValue`, with packed-array conversion retained.
   The helper classes (`OpenMMHandle`, `OpenMMStrings`, `OpenMMBooleans`, `OpenMMRuntime`, and
   `Vec3`) have no direct JNA façade counterparts. Their JavaDocs describe FFM ownership,
   encoding, initialization, and value semantics.
2. **Core lifecycle graph (complete):** FFM `Force`, `Integrator`, `VirtualSite`,
   `TabulatedFunction`, `System`, `Platform`, `Context`, and `State` now have explicit
   ownership contracts and method-level API documentation. `OpenMMFFMLifecycleTest` creates a
   Reference-platform context, verifies copied position/velocity state plus time and step
   values, and verifies idempotent context cleanup. A `Context` owns and destroys its
   integrator, matching the JNA façade lifecycle. `System.getForce()` and
   `System.getVirtualSite()` return only wrappers previously registered through the same FFM
   `System`; they intentionally return null for native handles introduced outside that wrapper.
   Registered OpenMM `Platform` instances are borrowed global handles: closing an FFM
   `Platform` wrapper invalidates only that wrapper and never destroys the native registry
   platform.

   A follow-up audit found no missing lifecycle operations. Inherited handle methods provide
   `getPointer()`, idempotent `destroy()`, and `close()`. The JNA constraint out-parameter
   methods are represented by `System.ConstraintParameters`; box-vector out parameters by
   `System.PeriodicBoxVectors`; and state parameter maps by borrowed FFM `MemorySegment`
   handles. Context string-pointer overloads accept `MemorySegment`, while ordinary Java callers
   use the safer UTF-8 `String` overloads. Native handle constructors and pointer rebinding use
   `MemorySegment` rather than exposing JNA types. JavaDocs were checked for these mappings,
   ownership/lifetime statements, and the JNA lifecycle behavior.
3. **Core forces and integrators (implementation complete; validation remains):** FFM
   `HarmonicBondForce` was the first
   completed ordinary force. It uses a typed `BondParameters` result record instead of JNA
   out-parameter objects, and its Reference-platform test verifies parameter round trips,
   harmonic energy, and `updateParametersInContext`. Adding a force to an FFM `System`
   transfers native ownership to that system; system destruction or force removal invalidates
   the associated Java wrapper to prevent a double native destroy. Port ordinary
   bonded/nonbonded forces, virtual sites, thermostats/barostats, tabulated functions,
   reporters, and integrators in similarly small batches. For each batch, compare JNA and FFM
   energy, forces, positions, velocities, and relevant parameters in separate JVMs using
   deterministic Reference-platform systems.

   Public method parity is audited for every JNA façade that currently has an FFM counterpart.
   JNA-specific pointer/out-parameter signatures are represented with `MemorySegment`, copied
   Java values, or typed records; legacy raw integer flags and pointer-rebinding entry points are
   retained where they are public API.    JavaDocs are checked against the corresponding JNA
   behavior and OpenMM units. Classes ported in Phase 3 batches receive the same audit.

   The conventional-bonded batch includes FFM `HarmonicAngleForce`,
   `PeriodicTorsionForce`, and `RBTorsionForce`. Their public APIs and JavaDocs have been ported
   and compared with the JNA façades and OpenMM headers. `OpenMMFFMConventionalBondedForceTest`
   verifies native construction and typed parameter round trips. Additional analytical
   energy/force parity and live-context update tests are deferred to the later testing pass and
   are not blockers for completing the implementation batch.

#### Phase 3 implementation record

All listed implementation batches and Batch 5 sub-batches are complete. Public API parity and
JavaDoc consistency have been reviewed against the JNA façades and OpenMM headers. Focused native
coverage exists for the implementation batches; broader analytical, live-context update, and
cross-backend behavioral comparisons remain part of validation, as noted below.

Each batch keeps its FFM façade classes in `ffx.openmm.ffm` and carries over reviewed
JavaDocs.

1. **Complete — conventional bonded forces:** `HarmonicAngleForce`, `PeriodicTorsionForce`,
   and `RBTorsionForce`. Typed parameter records replace multi-value native out parameters.
   Construction and parameter round trips are covered; analytical energies, forces, and
   `updateParametersInContext` remain for the deferred test pass.
2. **Complete — nonbonded and global force controls (implementation/API parity):**
   FFM `NonbondedForce`, `GBSAOBCForce`, `CMMotionRemover`, `AndersenThermostat`, and the Monte
   Carlo barostats (`MonteCarloBarostat`, `MonteCarloAnisotropicBarostat`,
   `MonteCarloFlexibleBarostat`, and `MonteCarloMembraneBarostat`) have their common parameter
   surfaces implemented with JavaDocs checked against their JNA façades and the headers. Particle,
   exception, PME, cutoff,
   switching, dielectric, and barostat-mode values use typed records/enums. Native tests cover
   common parameter round trips, all membrane Z modes, GBSA context updates, and a deterministic
   Reference-platform Lennard-Jones energy and live parameter update. Stochastic barostat
   trajectory behavior is intentionally not asserted.

   The remaining `NonbondedForce` C-wrapper operations are now exposed as well: LJ-PME and
   context-selected PME parameters, global parameters, particle/exception parameter offsets,
   reciprocal-space force-group selection, direct-space inclusion, and exception
   periodic-boundary controls. Typed Java results replace native out parameters. Further unit and
   analytical parity tests are deferred for a later testing pass.
3. **Complete — virtual sites and standard integrators:** `TwoParticleAverageSite`,
   `ThreeParticleAverageSite`, `OutOfPlaneSite`, `LocalCoordinatesSite`,
   `LangevinIntegrator`, `LangevinMiddleIntegrator`, `BrownianIntegrator`,
   `VariableVerletIntegrator`, and `VariableLangevinIntegrator`. FFM public methods and
   constructors are implemented with JavaDocs checked against `VirtualSite.h`,
   `LangevinMiddleIntegrator.h`, `BrownianIntegrator.h`, `VariableVerletIntegrator.h`, and
   `VariableLangevinIntegrator.h`. The native `LocalCoordinatesSite` constructors and getters
   are exposed with `IntArray`, `DoubleArray`, `Vec3`, and copied weight arrays. The legacy JNA
   array-based constructor does not match the C wrapper/header signature (it omits the required
   particle-index array and supplies a nonexistent z-weight array), so the FFM façade follows the
   actual native API instead of carrying that unusable signature forward. Native tests cover site
   weight/property round trips and integrator parameter round trips; stochastic trajectories and
   numerical integration parity remain deferred.
4. **Complete — advanced integration, tabulated functions, and reporting:** `CustomIntegrator`,
   `CompoundIntegrator`, `NoseHooverIntegrator`, `MinimizationReporter`, and
   `Continuous1D/2D/3DFunction` plus `Discrete1D/2D/3DFunction`. Value-returning FFM getters
   replace JNA output references; subsystem thermostats use the wrapper's actual `IntArray` and
   `BondArray` arguments. Added native tests cover continuous/discrete value and parameter
   round trips, CustomIntegrator global variables and expression metadata, CompoundIntegrator
   ownership transfer, Nose-Hoover property access, and reporter lifecycle.

   Two native-wrapper limitations are explicit in the façade: `CustomIntegrator.getComputationStep`
   is unavailable because the C wrapper declares `char**` outputs but writes C++ `std::string`
   objects to those addresses, making the declared ABI unsafe; and the C reporter factory creates
   a no-op native reporter, so Java callback overrides are not supported. Nose-Hoover's legacy
   JNA constructor order also differs from the native C wrapper's argument order; the FFM façade
   keeps the JNA Java parameter sequence and maps it to the documented native order. Analytical
   integration, callback, and minimization-progress behavior remain deferred pending wrapper-level
   support and focused behavioral tests.
5. **Custom and specialized forces:** `CustomAngleForce`, `CustomBondForce`,
   `CustomCentroidBondForce`, `CustomCompoundBondForce`, `CustomCVForce`,
   `CustomExternalForce`, `CustomGBForce`, `CustomHbondForce`, `CustomManyParticleForce`,
   `CustomNonbondedForce`, `CustomTorsionForce`, `CustomVolumeForce`, `CMAPTorsionForce`,
   `GayBerneForce`, `RMSDForce`, and `ATMForce`. This is the highest-risk batch: introduce
   typed result records for all parameter getters, retain dependent force/function wrappers
   while native parents use them, and add focused construction, parameter-update, expression,
   and nested-array lifetime tests per force. Work through it in these reviewable sub-batches:
   (a) **Complete — expression-driven scalar forces** (`CustomAngleForce`, `CustomBondForce`,
   `CustomTorsionForce`, and `CustomExternalForce`); (b) **Complete — multi-particle and
   interaction forces** (`CustomNonbondedForce`, `CustomGBForce`, `CustomHbondForce`, and
   `CustomManyParticleForce`);
   (c) **Complete — nested/composite forces** (`CustomCentroidBondForce`, `CustomCompoundBondForce`,
   `CustomCVForce`, and `CustomVolumeForce`); and (d) **Complete — specialized forces**
   (`CMAPTorsionForce`, `GayBerneForce`, `RMSDForce`, and `ATMForce`). Before each sub-batch,
   audit its JNA public
   constructors/methods, wrapper signatures, header semantics, ownership, and available generated
   symbols. Complete API/JavaDoc parity and focused native tests for a sub-batch before starting
   the next; defer full energy/force parity comparisons to the later validation pass.

   Sub-batch (a) is implemented in `ffx.openmm.ffm`. Typed result records replace the JNA
   out-reference and NIO buffer forms; parameter-vector inputs accept `double[]` and `DoubleArray`,
   with C strings handled through scoped UTF-8 strings (and `MemorySegment` equivalents for the
   legacy raw-pointer `CustomTorsionForce` overloads). Parameter counts/names/defaults,
   expressions, periodic-boundary controls available in the C wrapper, and context updates are
   exposed. Native tests verify construction, expression and parameter metadata, particle/index
   and parameter round trips, updates, and periodic-boundary properties for all four façades.

   Sub-batch (b) is implemented in `ffx.openmm.ffm`, with typed records/value-returning getters
   replacing JNA pointer out-parameters and NIO buffers. Array input accepts Java arrays, FFM
   array wrappers, or native `MemorySegment` handles; raw native string/handle overloads are
   retained where the JNA façade exposed them. The focused native tests cover construction,
   particle and group values, exclusions, interaction groups, type filters, force settings,
   and parameter updates. The C wrapper declares several custom-force string-output parameters
   as `char**` although their implementation writes C++ `std::string` objects; the affected
   function/computed-value getters fail explicitly with `UnsupportedOperationException` rather
   than calling that incompatible ABI. Fixing those wrappers or exposing safe C ABI accessors is
   required to restore those getters.

   Sub-batch (c) is implemented in `ffx.openmm.ffm`. Typed bond/group records and copied array
   values replace mutable JNA array-out usage, and raw native-handle overloads are available for
   pointer-based façade entry points. Collective-variable forces and tabulated functions invalidate
   their Java wrappers after ownership transfers to their native parent. Focused native tests cover
   group/bond and parameter round trips, expressions and settings, collective-variable access, and
   ownership transfer. `CustomCompoundBondForce.getFunctionParameters` remains explicitly
   unsupported because its generated `char**` output conflicts with the C++ `std::string` output
   used by the C wrapper; a safe wrapper ABI is needed to restore it.

   Sub-batch (d) is implemented in `ffx.openmm.ffm`. CMAP, Gay-Berne, RMSD, and ATM output
   parameters are exposed as typed records or copied Java arrays, with native array/handle
   overloads where appropriate. ATM nested-force ownership and RMSD borrowed-array handling are
   explicit. Focused native tests cover construction, specialized particle/map/exception values,
   setting round trips, and ownership. `GayBerneForce` maps its Java shape/strength parameter order
   to the actual C-wrapper order, which places strength values before shape values.

4. **AMOEBA and Drude:** Drude is implemented in `ffx.openmm.ffm.drude`: `DrudeForce`,
   `DrudeIntegrator`, `DrudeLangevinIntegrator`, `DrudeNoseHooverIntegrator`, and
   `DrudeSCFIntegrator` retain the JNA façade surface with typed FFM results. Focused native tests
   cover force parameter round trips, periodicity, and construction/property round trips for all
   integrators. The Drude Nose-Hoover C wrapper exposes one chain length, MTS count, and
   Yoshida-Suzuki count while the legacy Java signature labels the final integer arguments
   differently; the façade preserves its JNA constructor sequence and documents positional
   forwarding. Its kinetic-energy/temperature queries require an integrator attached to a live
   context. The JavaDoc audit against the JNA façades and upstream Drude headers is complete for
   all Drude FFM classes and the package description, including typed copied-value outputs and
   parameter documentation. The supported AMOEBA FFM façades are implemented in `ffx.openmm.ffm.amoeba`:
   `DoubleArray3D`, `VdwForce`, `WcaDispersionForce`, `GeneralizedKirkwoodForce`,
   `MultipoleForce`, `HippoNonbondedForce`, and `TorsionTorsionForce`. Typed records and
   `MemorySegment` arguments replace JNA references and buffers; arrays returned by native forces
   are documented as borrowed or copied into Java values. Native tests cover grid/torsion, VDW,
   WCA, generalized Kirkwood, multipole, and HIPPO parameter round trips. Two legacy surfaces are
   limited by the installed ABI: `GKCavitationForce` has no native C-wrapper symbols (the FFM
   façade reports this as unsupported instead of terminating the JVM), and
   `MultipoleForce.getCovalentMaps(int)` cannot safely return the legacy `IntArray` because the
   native function requires an opaque `OpenMM_2D_IntArray` for which the wrapper exports no
   creation/access functions; the FFM façade provides the raw-handle form and reports the
   unsupported convenience form explicitly. Complete context-level `updateParametersInContext`
   coverage and review native array ownership before integrating AMOEBA into potential-energy
   workflows. Class and method JavaDocs for all eight AMOEBA FFM classes were cross-checked
   against their JNA counterparts and the OpenMM headers in
   `openmm/openmm/plugins/amoeba/openmmapi/include/openmm/`; documented contracts include native
   units/defaults, enum meanings, context-update scope, and borrowed/copied handle lifetimes.
5. **Downstream opt-in and retirement:** add an all-FFM path in
   `ffx.potential.openmm`, keep JNA as the default through a laboratory validation period,
   and collect comparison results on supported macOS/Linux configurations. Remove JNA,
   `jopenmm-fat`, and the backend selector only after all supported workflows use FFM and the
   FFM result matrix meets agreed tolerances.

### FFM JavaDoc audit

JavaDoc review is complete for the FFM package, including the core lifecycle/container façades,
standard and specialized forces and integrators, custom forces, tabulated functions, and
reporters. Documentation was compared against matching JNA façades and the OpenMM core, AMOEBA,
and Drude headers. FFM-specific handle ownership, borrowed/copied values, confined-memory
lifetimes, units, defaults, enum values, and update-in-context restrictions are described.
Documented compatibility caveats include unsafe `char**`/C++ `std::string` outputs, unavailable
reporter callbacks, unexposed covalent-map array accessors, deprecated native APIs retained for
compatibility, and legacy Java signatures that differ from native wrapper meanings. No executable
API changes are part of this documentation pass.

### High-risk conversions

| Area | Classes/examples | FFM work and validation focus |
|---|---|---|
| Lifecycle and borrowed handles | `System`, `Context`, `State`, `Force`, `Integrator`, `Platform` | OpenMM returns a mix of owned and borrowed pointers. Model ownership explicitly so a returned system/force/state is not destroyed twice or outlives its owner. `System` assumes native ownership of added forces and virtual sites and invalidates their FFM wrappers when it releases them; registered `Platform` handles remain native-global and must never be destroyed. Preserve the current context/integrator destruction order. |
| Struct values and borrowed arrays | `Vec3Array`, `State`, `Context` | Replace JNA `OpenMM_Vec3` structures with generated layouts and `MemorySegment` accessors. Distinguish returned borrowed `Vec3Array`/parameter data from Java-owned arrays before copying values. |
| Strings, C legacy declarations, and plugins | `Platform`, force name/expression APIs | Allocate UTF-8 strings in a confined arena; copy returned C strings before the arena closes. `OpenMM_Platform_getPluginLoadFailures()` uses legacy C `()` syntax, so its generated binding is variadic. Preserve the runtime's explicit plugin diagnostics. |
| Many out parameters | `Custom*Force`, `DrudeForce`, `MultipoleForce`, `HippoNonbondedForce` | Replace JNA references and NIO output buffers with typed records/value objects backed by confined arenas. Use named result records rather than wide public `MemorySegment` output signatures. |
| Nested native arrays | `DoubleArray3D`, `TorsionTorsionForce`, AMOEBA covalent maps and multipoles | Specify whether Java or OpenMM owns every nested array, retain child handles through native calls, and test all array dimensions, resizing, and destruction. |
| Extension symbol overlap | `DoubleArray3D` across AMOEBA/Drude | The aggregate binding resolves shared symbols in library lookup order. Verify the C ABI is identical for every overlapping symbol and add extension-specific smoke tests before sharing a wrapper. |
| jextract code-generation quirks | `GayBerneForce`, generated plugin-failure API | Retain the scoped `ex` macro workaround: OpenMM's Gay-Berne parameters named `ex` collide with jextract's generated Java catch variable. Treat each new jextract version as an ABI/code-generation update requiring regeneration, compilation, and parity tests. |

## Implementation steps

1. **Runtime loading (bootstrap implemented; packaged distribution pending)**
   - Inventory the files and plugin libraries currently supplied by `jopenmm-fat`; create FFX
     platform-classified native artifacts/resources for macOS and Linux from the matching
     OpenMM installation rather than relying on JNA extraction.
   - Add one `ffx.openmm.ffm.OpenMMRuntime` bootstrapper that detects the supported OS/CPU,
     extracts libraries to a versioned temporary cache when necessary, loads `libOpenMM`
     before AMOEBA and Drude, and exposes the plugin directory. Make loading explicit and
     fail with a diagnostic that identifies the missing platform/library; do not silently
     fall back to JNA.
   - Initialize the runtime before the generated bindings are first referenced, call
     `OpenMM_Platform_loadPluginsFromDirectory` through the generated core binding, and
     surface plugin-load failures. Replace the current `OpenMMUtils.init()` and
     `OpenMMUtils.getPluginDirectory()` calls in tests and consumers.

2. **FFM adapter layer (implemented)**
   - Place hand-written helpers in `ffx.openmm.ffm`, separate from generated code:
     `OpenMMRuntime`, `OpenMMHandles`/ownership utilities, `OpenMMStrings`, `OpenMMStructs`,
     and narrowly scoped out-parameter helpers.
   - Represent all opaque OpenMM handles as non-NULL `MemorySegment` addresses and preserve
     each façade's current destroy/idempotency semantics. Use explicit `destroy()` as today;
     optionally add `AutoCloseable` only if all call sites can be migrated safely.
   - Use `Arena.ofConfined()` around each temporary C string, primitive out parameter, and
     short-lived vector structure. Decode returned `const char*` via UTF-8 before the arena
     closes. Keep native-owned strings/arrays as borrowed segments and never close or free
     them from Java.
   - Centralize `OpenMM_Vec3` layout construction/access and pointer-to-struct argument
     creation, replacing `OpenMM_Vec3`. Centralize boolean conversion (`int`/enum 0 or 1)
     and `long long`/Java `long` conversion. Replace JNA `Pointer`, `PointerByReference`,
     `IntByReference`, `DoubleByReference`, and NIO output-buffer overloads with strongly
     typed value-returning methods or `MemorySegment` adapters only where public
     compatibility requires them.

3. **Port façade classes (implementation complete; final parity validation remains)**
   - First migrate ownership roots and shared value containers: `Force`, `Integrator`,
     `VirtualSite`, `TabulatedFunction`, `System`, `Context`, `State`, `Platform`,
     `Vec3Array`, `DoubleArray`, `IntArray`, `StringArray`, `BondArray`, and `IntSet`.
     Replace the JNA-`Pointer` keyed identity maps in `System` with stable `MemorySegment`
     address keys.
   - `ffx.openmm.ffm.BondArray` is the initial standalone conversion prototype. It demonstrates opaque
     `MemorySegment` handle ownership and `Arena`-allocated integer out parameters before
     replacing the public JNA-backed `BondArray`.
   - Migrate the ordinary core forces, virtual sites, thermostats/barostats, integrators,
     tabulated functions, reporters, and serialization users by changing only their native
     call sites to the generated `OpenMMNative` methods plus FFM helpers.
   - Migrate `ffx.openmm.amoeba` and `ffx.openmm.drude` to their generated symbols in the
     aggregate `OpenMMNative` binding, including all pointer/out-parameter-heavy multipole and
     Drude APIs. Confirm every one of the 1,105 referenced `OpenMM_*` symbols is generated or
     deliberately removed as unused.
   - Update downstream OpenMM application classes in `modules/potential` and
     `modules/algorithms` only for intentional public type/API changes. Preserve simulation
     semantics and explicit native-resource lifecycle ordering.

4. **Dependencies and build integration (distribution work pending)**
   - Keep `jopenmm-fat` and JNA dependencies while the FFM backend is opt-in and JNA remains
     the default. Remove them from `modules/openmm/pom.xml` and then `modules/potential` only
     at the retirement gate defined above.
   - Generated bindings already reside in the standard Maven Java source tree. Complete
     native-resource/artifact packaging in the relevant Maven modules and ensure plugin
     libraries and dylib/so dependency paths work from IDE output and packaged distributions.
   - Update module documentation and configuration examples with required Java 25+, supported
     OS/architectures, native location override, plugin loading behavior, and the regeneration
     workflow.

5. **Validate correctness and decide on removal (in progress)**
   - Focused tests already cover FFM container lifecycle, strings, `Vec3`, booleans, out
     parameters, façade parameter round trips, and selected Reference-platform energies and
     context updates. Expand bootstrap/library-resolution and cross-backend tests, and fill
     the documented gaps in analytical energy/force and `updateParametersInContext` coverage.
   - Keep the original JNA `OpenMMHarmonicBondTest` unchanged during the parallel-backend
     period; retain its FFM counterpart and compare both in separate processes on the
     Reference platform. Add/verify core, AMOEBA, and Drude native smoke coverage through the
     aggregate generated `OpenMMNative` binding.
   - Run the OpenMM module target and affected potential/algorithm tests on macOS and Linux.
     Verify packaged plugin loading and native dependency resolution. Remove JNA and
     `jopenmm-fat` only after the downstream FFM workflows and parity matrix pass; then confirm
     no `com.sun.jna` or `edu.uiowa.jopenmm` imports remain in retired code.

## Phase 2: `ffx.potential.ommffm` port (planned)

### Findings from the existing `ffx.potential.openmm` package

- 35 classes, ~10.5k lines in `modules/potential`. Core: `OpenMMEnergy` (extends `ForceFieldEnergy`,
  implements `OpenMMPotential`), `OpenMMSystem` (extends `ffx.openmm.System`), `OpenMMContext`,
  `OpenMMState`, `OpenMMIntegrator`, `CustomMTS*Integrator`, and the dual-topology pair.
- ~23 force classes (`BondForce`, `AngleForce`, `Amoeba*Force`, `FixedCharge*`, `Restrain*`, ...) each
  **extend a JNA façade** (e.g. `BondForce extends CustomBondForce`) and expose
  `constructForce(OpenMMEnergy)`, `updateForce(OpenMMEnergy)` and dual-topology variants.
- Entry point: `ForceFieldEnergy.energyFactory` instantiates `OpenMMEnergy` for `PLATFORM` =
  `OMM`, `OMM_REF`, `OMM_CUDA`, `OMM_OPENCL`. Consumers in `modules/algorithms`
  (`MolecularDynamics`, `MolecularDynamicsOpenMM`, `Minimize`, `MinimizeOpenMM`, `TitrationManyBody`,
  `RotamerOptimization`, `EnergyExpansion`) and `potential` (`DualTopologyEnergy`, `LambdaGradient`)
  use `instanceof OpenMMEnergy` and the JNA-typed `OpenMMPotential` interface.
- Existing tests: only `ParentEnergyTest` touches OpenMM. There is no unit-level regression
  safety net, so parity tests must be written first (see below).

### Improvements to the proposed outline

1. **Subclass vs. compose.** The FFM façades are `final`, but the JNA forces subclass them. Decision:
   remove `final` only from the façades that ommffm extends (`System`, `CustomBondForce`, other
   `Custom*Force`, `AmoebaMultipoleForce`, ...), matching the `NoseHooverIntegrator` precedent. This
   keeps each port a mechanical change of imports and native calls. Do not introduce composition
   during the port; refactor afterward.
2. **Do not copy blindly: extract shared logic first.** Parameter/topology gathering (which atoms and
   bonds, units conversion, `use`/lambda flags, dual-topology indexing) is backend independent. Move
   it to backend-neutral helpers in `ffx.potential.ommffm` (or a shared package) as the force classes
   are ported, and leave the JNA classes untouched. Do not edit `ffx.potential.openmm` except to
   delete it at retirement.
3. **Backend-neutral seam before the first energy.** `instanceof OpenMMEnergy` and the JNA-typed
   `OpenMMPotential` block the algorithms layer. Add a small neutral interface (e.g.
   `ffx.potential.OpenMMBackedPotential`: `updateContext`, `setActiveAtoms`, `updateParameters`,
   `getTemperature`, ...) in a later phase; the first phases only need `energy` and
   `energyAndGradient`, so skip MD and Minimize until phase F.
4. **Backend selection.** Keep `PLATFORM` values unchanged; select the implementation with
   `ffx.openmm.backend=jna|ffm` (default `jna`) inside `energyFactory`. A JVM must use only one
   backend and never hand native handles across; compare in separate JVMs.
5. **Lifecycle.** `OpenMMEnergy.destroy()` and `OpenMMSystem.free()` map to `AutoCloseable`
   façades. `System.addForce` transfers ownership; after adding, the Java wrapper is invalid, so
   ommffm force classes must retain only Java-side parameters needed by `updateForce` and use the
   System's lookup (`getForce(index)`/typed accessors) for updates. Verify each updater this way
   (`updateParametersInContext`).
6. **Bulk transfer.** Positions/gradients move per atom through `Vec3Array`/State records. Add one
   helper pair that converts `double[] x` (Å) <-> `Vec3Array` (nm) and State forces -> gradient
   (kcal/mol/Å) in one pass, with the unit constants in a single class shared with the JNA version
   (the numeric constants must stay identical).
7. **Numerics gate.** Reference platform in double precision: |dE| < 1e-8 kcal/mol, per-atom gradient
   RMS diff < 1e-8; CUDA/OpenCL/CPU only a looser agreed tolerance. Run JNA and FFM in separate JVMs
   writing energies/gradients to a file, then compare. Include a displaced geometry, since
   equilibrium geometry hides gradient/sign errors.
8. **Known façade gaps to design around.** String-output getters, `OpenMM_2D_IntArray`
   (`MultipoleForce.getCovalentMaps(int)`), `GKCavitationForce` and reporter callbacks are
   unsupported in FFM. Ensure no ommffm path needs them (use the Java-side copy of the data); keep
   `AmoebaGKCavitationForce` last and behind the missing-symbol limitation (see "High-risk conversions").
9. **Platform plugins and loading.** `ffx.potential.ommffm.OpenMMContext.loadPlatform` must call
   `OpenMMRuntime` and honor `FFX_OPENMM_LIB_DIR` / `FFX_OPENMM_PLUGIN_DIR`, plus the existing
   CUDA/OpenCL device property logic (`getDefaultDevice`, `CUDA_DEVICES`).
10. **JPMS/module wiring.** Confirm `modules/potential` lists `ffx-openmm` and add
    `--enable-native-access=ALL-UNNAMED` (or the module name) to surefire/ffx launch scripts so FFM
    produces no restricted-method warnings.

### Phases (each ends with a green parity test and a checkbox here)

**Phase A – scaffold (bond-only energy).** [Completed]
- [x] Create `ffx.potential.ommffm` with `package-info`, `OpenMMContext`, `OpenMMSystem`
  (`addParticle`/masses only plus `addForces()` limited to bonds), `OpenMMState`, `BondForce`,
  `OpenMMEnergy` (energy and gradient only; the rest of the `ForceFieldEnergy` API is inherited).
- [x] Reduced FFM `OpenMMEnergy` is reachable only by direct construction, not yet by the factory.
- [x] Test 1: `OpenMMFFMBondEnergyTest` verifies equilibrium and displaced energies/gradients
  against analytic values and against pure-Java `ForceFieldEnergy`.
- [x] Defined `ffx.openmm.ffm.OpenMMUnits` to provide exact compile-time OpenMM unit conversion constants.
- [x] Exit: 1e-8 kcal/mol on Reference; close/destroy leaves no native leaks.

**Phase B – all bonded terms.** [Completed]
- [x] Implemented all bonded interactions and restraint forces in `ffx.potential.ommffm`:
  `AngleForce`, `InPlaneAngleForce`, `UreyBradleyForce`, `StretchBendForce`, `OutOfPlaneBendForce`,
  `TorsionForce`, `ImproperTorsionForce`, `PiOrbitalTorsionForce`, `StretchTorsionForce`,
  `AngleTorsionForce`, and `AmoebaTorsionTorsionForce` (with 3D bicubic splines via `DoubleArray3D`).
- [x] Added restraint terms: `RestrainPositionsForce`, `RestrainDistanceForce`, `RestrainTorsionsForce`,
  and `RestrainGroupsForce`.
- [x] Adopted exception-safe `try-with-resources` buffer allocation patterns across all force classes.
- [x] Exit: Validated against FFX pure-Java energy with full parameter update and gradient agreement.

**Phase C – fixed-charge nonbonded.** [Completed]
- [x] Implemented `FixedChargeNonbondedForce` (Lennard-Jones, Coulomb, 1-4 exceptions, cutoffs, PME).
- [x] Implemented `FixedChargeGBForce` (continuum OBC/ACE Generalized Born solvation).
- [x] Implemented `FixedChargeAlchemicalForces` (softcore sterics decoupling and scaled 1-4 interactions).
- [x] Wired `OpenMMSystem` controls: `AndersenThermostat`, `MonteCarloBarostat`, `CMMotionRemover`,
  periodic boundary conditions (`setPeriodicBoxVectors`), and dynamic $\lambda$ updates.
- [x] Exit: `OpenMMFFMFixedChargeNonbondedEnergyTest` verifies gas-phase and periodic energy/gradient parity.

**Phase D – AMOEBA.** [Completed]
- [x] Implemented `AmoebaVdwForce` (buffered 14-7 vdW, H-reduction, dispersion correction, alchemical $\lambda$-sterics).
- [x] Implemented `AmoebaMultipoleForce` (multipole moments, local frames, mutual/extrapolated polarization,
  1-2 through 1-5 covalent/polarization masks, PME electrostatics).
- [x] Implemented `AmoebaGeneralizedKirkwoodForce` and `AmoebaWcaDispersionForce` for polarizable continuum solvation.
- [x] `AmoebaGKCavitationForce` handled as an explicit unsupported stub matching C wrapper symbol limits.
- [x] Exit: `OpenMMFFMAmoebaEnergyTest` passes on water box and polypeptide systems on Reference platform.

**Phase E – integrators and State.** [Completed]
- [x] Create `OpenMMIntegrator` factory in `ffx.potential.ommffm`.
- [x] Implement `CustomMTSIntegrator` and `CustomMTSLangevinIntegrator` (Java-side step tracking is low-priority; keep simple).
- [x] Extend `OpenMMContext` and `OpenMMState` for particle velocity management, kinetic energy, and state polling (buffer reuse is optional/minimal since ~1000 native steps run between Java synchronizations).
- [x] Test: `OpenMMFFMIntegratorDynamicsTest` for NVE energy conservation and NVT temperature regulation.

**Phase F – consumers and dual topology.** [Completed]
- [x] Step 1: Dual Topology: `OpenMMDualTopologySystem`, `OpenMMDualTopologyEnergy`, and dual-topology overloads on all force classes (`BondForce`, `AngleForce`, `InPlaneAngleForce`, `StretchBendForce`, `UreyBradleyForce`, `OutOfPlaneBendForce`, `TorsionForce`, `ImproperTorsionForce`, `PiOrbitalTorsionForce`, `StretchTorsionForce`, `AngleTorsionForce`, `RestrainTorsionsForce`, `AmoebaTorsionTorsionForce`, `AmoebaVdwForce`, `AmoebaMultipoleForce`).
- [x] Test: `OpenMMFFMDualTopologyEnergyTest` verifying $E(\lambda)$ parity across intermediate $\lambda$ values against pure-Java `DualTopologyEnergy`.
- [x] Step 2: Backend Switch & Consumer Algorithm Wiring (wired `MolecularDynamicsOpenMM`, `MinimizeOpenMM`, and `ffx.openmm.backend=jna|ffm` in `ForceFieldEnergy` and `DualTopologyEnergy` factories).
- [x] Test: `OpenMMFFMConsumerTest` verifying minimization and molecular dynamics execution using `-Dffx.openmm.backend=ffm`.

**Phase G – Lab Validation & Coexistence Period.** [Active]
- [x] Backend switch (`ffx.openmm.backend=jna|ffm`, default `jna`) fully wired into `ForceFieldEnergy.energyFactory`, `DualTopologyEnergy.energyFactory`, `MinimizeOpenMM`, and `MolecularDynamicsOpenMM`.
- [x] High-level algorithm compatibility (`MolecularDynamicsOpenMM` bridge, `Minimize`, `RotamerOptimization`, `TitrationManyBody`, `EnergyExpansion`).
- [x] End-to-end integration tests (`OpenMMFFMConsumerTest`) passing for minimization and Verlet/Langevin molecular dynamics.
- [ ] Laboratory validation period: Enable lab members to test their production simulations, free energy calculations (e.g. BAR dual-topology, OST), and GPU platforms (CUDA, OpenCL) using `-Dffx.openmm.backend=ffm` or property `ffx.openmm.backend=ffm`.
- [ ] Review performance benchmarks and collect lab feedback across real-world workflows before deprecation.

**Phase H – Final JNA Backend Retirement.** [Deferred pending lab testing]
- [ ] Flip the default backend property to `ffm` (`ffx.openmm.backend=ffm` as default).
- [ ] Deprecate legacy JNA packages and classes (`ffx.potential.openmm`, `ffx.openmm` JNA façades).
- [ ] Remove `jopenmm-fat` and JNA dependencies from Maven build files (`modules/openmm/pom.xml`, `modules/potential/pom.xml`).
- [ ] Consolidate packages (optionally promote `ffx.potential.ommffm` to `ffx.potential.openmm` and `ffx.openmm.ffm` to `ffx.openmm`).

### Decisions (confirmed)

1. FFM façades that `ommffm` subclasses lose `final`; no composition during the port.
2. Backend selection is the simple `ffx.openmm.backend=jna|ffm` property read in `energyFactory` (defaulting to `jna`);
   no new `Platform` enum value. JNA will remain supported as the default throughout the lab validation period,
   and will only be removed once lab members have thoroughly tested and confirmed the FFM backend on their production workloads.
3. Testing targets macOS with the Reference platform only. Linux/CUDA runs come later.
4. Tests are **not** parameterized over backends. New tests are FFM-only and compare against
   analytic values, the pure-Java FFX energy (`ForceFieldEnergy`), or stored reference
   energies/gradients captured once from the JNA implementation. The separate-JVM JNA parity
   runs described above are optional one-time captures, not an ongoing test matrix.
   `ParentEnergyTest` is left unchanged.
5. In molecular dynamics, OpenMM runs natively for typically ~1000 steps per batch before returning
   control to Java for logging and trajectory I/O. Therefore, temporary buffer allocation overhead
   during state querying is negligible; buffer reuse in `OpenMMState` and `OpenMMContext` should only
   be implemented where trivial and maintainable.
6. For `CustomMTSIntegrator` and `CustomMTSLangevinIntegrator`, maintaining Java-side records of
   computation steps is low priority and should be kept minimal and simple.
7. `AmoebaGKCavitationForce` exists in local OpenMM source but has not yet been merged into upstream
   OpenMM repository distributions. The current Java stub that safely returns `null` from
   `constructForce()` is an intentional near-term workaround until upstream repository integration.
8. For AMOEBA dual-topology electrostatics, permanent multipoles in topology $t$ scale as $\sqrt{S_t}$
   (where $S_t$ is the topology scale factor) to guarantee exact linear scaling of the pairwise
   Coulomb/multipole interaction energy ($E \propto q_i q_j \propto \sqrt{S_t}\sqrt{S_t} = S_t$), with
   covalent exclusion masks (1-2 through 1-5 and 1-1 polarization groups) mapped and merged across both
   topologies using `atom.getTopologyAtomIndex()`.
9. Consumer integration & bridge architecture: High-level algorithms (`MolecularDynamicsOpenMM`, `MinimizeOpenMM`,
   `Minimize`, `MolecularDynamics`, `RotamerOptimization`, `TitrationManyBody`, `EnergyExpansion`) support both JNA
   and FFM backends polymorphically using bridge wrappers and dynamic dispatch to enable uninterrupted side-by-side
   comparison without disrupting legacy pipelines.

## Notes and considerations

- `OpenMMCWrapper.h`, `AmoebaOpenMMCWrapper.h`, and `DrudeOpenMMCWrapper.h` are the complete
  C binding inputs for this installation. `OpenMMAmoeba.h` and `OpenMMDrude.h` must not be
  generator inputs because `jextract` consumes C declarations, not OpenMM's C++ public API.
- The generated source is platform-independent and contains one aggregate `OpenMMNative`
  binding. A Linux release must verify it against version-matched headers and load
  `libOpenMM.so`, `libOpenMMAmoeba.so`, and `libOpenMMDrude.so`; native paths remain
  platform-specific runtime configuration.
- Keep generated bindings mechanically generated and hand-written policy/lifecycle code
  outside them so OpenMM upgrades are reviewable. Record the OpenMM ABI version and reject
  incompatible headers/libraries at generation or bootstrap time.
- The migration must not mix JNA and FFM in a completed façade: it risks distinct library
  instances, incompatible handle ownership, and invalid cross-API pointers. A temporary
  branch may port in batches, but the final dependency cleanup gates completion.
