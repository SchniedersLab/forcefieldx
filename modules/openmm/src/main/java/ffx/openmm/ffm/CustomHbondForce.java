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

import ffx.openmm.ffm.bindings.OpenMMNative;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * Custom interaction evaluated between each eligible donor and acceptor group.
 *
 * <p>The expression may use arbitrary distances, angles, and dihedrals among the six group
 * particle slots, written with OpenMM's {@code distance(p1,p2)}, {@code angle(p1,p2,p3)}, and
 * {@code dihedral(p1,p2,p3,p4)} functions, along with global, per-donor, and per-acceptor
 * parameters and tabulated functions. Distances are nm, angles/dihedrals radians, and energy
 * kJ/mol. Each donor/acceptor parameter array must contain exactly the values declared for its
 * corresponding role, in declaration order. Nonbonded method enum values are 0 {@code NoCutoff},
 * 1 {@code CutoffNonPeriodic}, and 2 {@code CutoffPeriodic}; cutoff distance is in nm.
 *
 * <p>The native header deprecates legacy continuous-function APIs in favor of tabulated
 * functions; matching Java overloads remain for compatibility.
 *
 * <p>Java arrays are copied to temporary native storage; native array inputs are borrowed
 * synchronously. Returned records contain copied parameter arrays. Adding a tabulated function
 * transfers ownership; the typed overload invalidates its wrapper, while a raw function handle
 * cannot invalidate an associated wrapper and must not be destroyed independently. Returned
 * handles are borrowed.
 * {@link #getFunctionParameters(int)} is deliberately unsupported because the wrapper's
 * {@code char**} output is implemented by writing C++ {@code std::string}, an unsafe incompatible
 * ABI. Updating an existing context affects supported donor/acceptor values and tabulated
 * functions, not the expression, role declarations, group topology, exclusions, or global
 * defaults; an uninitialized context is ignored.
 */
public class CustomHbondForce extends Force {

  /**
   * Snapshot of a donor or acceptor's three ordered particle slots and role-specific values.
   *
   * @param particle1 first slot (donor's hydrogen or acceptor atom as defined by the force setup)
   * @param particle2 second slot
   * @param particle3 third slot
   * @param parameters copied per-donor or per-acceptor values in declaration order
   */
  public record GroupParameters(int particle1, int particle2, int particle3, double[] parameters) {}
  /** @param donor donor-list index
   *  @param acceptor acceptor-list index */
  public record Exclusion(int donor, int acceptor) {}
  /** Legacy function snapshot; function retrieval is unavailable through the unsafe wrapper ABI.
   * @param name function name
   * @param function native function handle
   * @param min lower function-domain bound
   * @param max upper function-domain bound */
  public record FunctionParameters(String name, MemorySegment function, double min, double max) {}

  /** Create a force from an expression over donor/acceptor geometry and declared parameters.
   * @param energy custom interaction energy expression */
  public CustomHbondForce(String energy) { super(create(energy)); }

  /** Add an acceptor's three ordered particle slots and its role-specific values.
   * @param a1 first acceptor particle slot
   * @param a2 second acceptor particle slot
   * @param a3 third acceptor particle slot
   * @param parameters per-acceptor values in declaration order; Java array copied
   * @return acceptor-list index */
  public int addAcceptor(int a1,int a2,int a3,double[] parameters) {
    try(DoubleArray values=CustomForceParameters.toNative(parameters)){return addAcceptor(a1,a2,a3,values);}
  }
  /** Add an acceptor using a native parameter array borrowed during the call.
   * @param a1 first particle slot
   * @param a2 second particle slot
   * @param a3 third particle slot
   * @param parameters per-acceptor values in declaration order
   * @return acceptor-list index */
  public int addAcceptor(int a1,int a2,int a3,DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomHbondForce_addAcceptor(getPointer(),a1,a2,a3,parameters.getPointer());
  }
  /** Add an acceptor using a caller-owned native segment borrowed during the call.
   * @param a1 first particle slot
   * @param a2 second particle slot
   * @param a3 third particle slot
   * @param parameters native per-acceptor values
   * @return acceptor-list index */
  public int addAcceptor(int a1,int a2,int a3,MemorySegment parameters) {
    return OpenMMNative.OpenMM_CustomHbondForce_addAcceptor(getPointer(),a1,a2,a3,parameters);
  }
  /** Add a donor's three ordered particle slots and role-specific values.
   * @param d1 first donor particle slot
   * @param d2 second donor particle slot
   * @param d3 third donor particle slot
   * @param parameters per-donor values in declaration order; Java array copied
   * @return donor-list index */
  public int addDonor(int d1,int d2,int d3,double[] parameters) {
    try(DoubleArray values=CustomForceParameters.toNative(parameters)){return addDonor(d1,d2,d3,values);}
  }
  /** Add a donor using a native parameter array borrowed during the call.
   * @param d1 first particle slot
   * @param d2 second particle slot
   * @param d3 third particle slot
   * @param parameters per-donor values in declaration order
   * @return donor-list index */
  public int addDonor(int d1,int d2,int d3,DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomHbondForce_addDonor(getPointer(),d1,d2,d3,parameters.getPointer());
  }
  /** Add a donor using a caller-owned native segment borrowed during the call.
   * @param d1 first particle slot
   * @param d2 second particle slot
   * @param d3 third particle slot
   * @param parameters native per-donor values
   * @return donor-list index */
  public int addDonor(int d1,int d2,int d3,MemorySegment parameters) {
    return OpenMMNative.OpenMM_CustomHbondForce_addDonor(getPointer(),d1,d2,d3,parameters);
  }
  /** Exclude a donor/acceptor list-index pair from interaction evaluation.
   * @param donor donor-list index
   * @param acceptor acceptor-list index
   * @return exclusion index */
  public int addExclusion(int donor,int acceptor){return OpenMMNative.OpenMM_CustomHbondForce_addExclusion(getPointer(),donor,acceptor);}
  /** Declare one per-donor value used by the expression.
   * @param name parameter name
   * @return declaration index */
  public int addPerDonorParameter(String name){return withStringResult(name,n->OpenMMNative.OpenMM_CustomHbondForce_addPerDonorParameter(getPointer(),n));}
  /** Declare a per-donor value from a caller-owned NUL-terminated UTF-8 name.
   * @param name name segment borrowed during the call
   * @return declaration index */
  public int addPerDonorParameter(MemorySegment name){return OpenMMNative.OpenMM_CustomHbondForce_addPerDonorParameter(getPointer(),name);}
  /** Declare one per-acceptor value used by the expression.
   * @param name parameter name
   * @return declaration index */
  public int addPerAcceptorParameter(String name){return withStringResult(name,n->OpenMMNative.OpenMM_CustomHbondForce_addPerAcceptorParameter(getPointer(),n));}
  /** Declare a per-acceptor value from a caller-owned NUL-terminated UTF-8 name.
   * @param name name segment borrowed during the call
   * @return declaration index */
  public int addPerAcceptorParameter(MemorySegment name){return OpenMMNative.OpenMM_CustomHbondForce_addPerAcceptorParameter(getPointer(),name);}
  /** Declare a global parameter and the default used in new contexts.
   * @param name expression parameter name
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(String name,double defaultValue){return withStringResult(name,n->OpenMMNative.OpenMM_CustomHbondForce_addGlobalParameter(getPointer(),n,defaultValue));}
  /** Declare a global parameter from a caller-owned NUL-terminated UTF-8 name.
   * @param name name segment borrowed for the call
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(MemorySegment name,double defaultValue){return OpenMMNative.OpenMM_CustomHbondForce_addGlobalParameter(getPointer(),name,defaultValue);}
  /** Add a tabulated function and transfer its native ownership to this force.
   * @param name expression-visible name
   * @param function wrapper invalidated after transfer
   * @return function index */
  public int addTabulatedFunction(String name,TabulatedFunction function){
    int index=withStringResult(name,n->OpenMMNative.OpenMM_CustomHbondForce_addTabulatedFunction(getPointer(),n,function.getPointer()));
    function.invalidate();return index;
  }
  /** Add a function from a caller-owned name segment and native function handle. Native function
   * ownership transfers to this force, but an associated wrapper is not invalidated automatically.
   * @param name NUL-terminated UTF-8 name, borrowed during the call
   * @param function native function handle
   * @return function index */
  public int addTabulatedFunction(MemorySegment name,MemorySegment function){
    return OpenMMNative.OpenMM_CustomHbondForce_addTabulatedFunction(getPointer(),name,function);
  }
  /** Add a legacy continuous one-dimensional function (deprecated by the native API); prefer
   * {@link #addTabulatedFunction(String, TabulatedFunction)}.
   * @param name expression-visible name
   * @param values sampled values
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(String name,double[] values,double min,double max){
    try(DoubleArray nativeValues=CustomForceParameters.toNative(values)){
      return withStringResult(name,n->OpenMMNative.OpenMM_CustomHbondForce_addFunction(getPointer(),n,nativeValues.getPointer(),min,max));
    }
  }
  /** Add a continuous function from caller-owned name and value-array segments borrowed
   * synchronously.
   * @param name NUL-terminated UTF-8 name
   * @param values sampled values
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(MemorySegment name,MemorySegment values,double min,double max){
    return OpenMMNative.OpenMM_CustomHbondForce_addFunction(getPointer(),name,values,min,max);
  }
  @Override public void destroy(){destroy(OpenMMNative::OpenMM_CustomHbondForce_destroy);}

  /** @param index acceptor-list index
   *  @return three ordered particle slots and copied per-acceptor values */
  public GroupParameters getAcceptorParameters(int index){return getGroupParameters(index,false);}
  /** @param index donor-list index
   *  @return three ordered particle slots and copied per-donor values */
  public GroupParameters getDonorParameters(int index){return getGroupParameters(index,true);}
  /** @param index exclusion index
   *  @return donor and acceptor list indices */
  public Exclusion getExclusionParticles(int index){
    try(Arena arena=Arena.ofConfined()){
      MemorySegment donor=arena.allocate(ValueLayout.JAVA_INT),acceptor=arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomHbondForce_getExclusionParticles(getPointer(),index,donor,acceptor);
      return new Exclusion(donor.get(ValueLayout.JAVA_INT,0),acceptor.get(ValueLayout.JAVA_INT,0));
    }
  }
  /**
   * Unsupported: this C wrapper's {@code char**} output is implemented by writing a C++
   * {@code std::string}; that is an incompatible ABI and cannot be read safely.
   *
   * @param index function index
   * @return never returns
   * @throws UnsupportedOperationException always because of the output ABI mismatch
   */
  public FunctionParameters getFunctionParameters(int index){throw invalidStringOutput("getFunctionParameters");}
  /** @return copied current energy expression */
  public String getEnergyFunction(){return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomHbondForce_getEnergyFunction(getPointer()));}
  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index){return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomHbondForce_getGlobalParameterName(getPointer(),index));}
  /** @param index global-parameter index
   *  @return default for newly created contexts */
  public double getGlobalParameterDefaultValue(int index){return OpenMMNative.OpenMM_CustomHbondForce_getGlobalParameterDefaultValue(getPointer(),index);}
  /** @param index per-donor declaration index
   *  @return copied parameter name */
  public String getPerDonorParameterName(int index){return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomHbondForce_getPerDonorParameterName(getPointer(),index));}
  /** @param index per-acceptor declaration index
   *  @return copied parameter name */
  public String getPerAcceptorParameterName(int index){return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomHbondForce_getPerAcceptorParameterName(getPointer(),index));}
  /** @param index function index
   *  @return borrowed force-owned handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index){return OpenMMNative.OpenMM_CustomHbondForce_getTabulatedFunction(getPointer(),index);}
  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index){return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomHbondForce_getTabulatedFunctionName(getPointer(),index));}

  /** @return donor-group count */ public int getNumDonors(){return OpenMMNative.OpenMM_CustomHbondForce_getNumDonors(getPointer());}
  /** @return acceptor-group count */ public int getNumAcceptors(){return OpenMMNative.OpenMM_CustomHbondForce_getNumAcceptors(getPointer());}
  /** @return donor/acceptor exclusion count */ public int getNumExclusions(){return OpenMMNative.OpenMM_CustomHbondForce_getNumExclusions(getPointer());}
  /** @return per-donor parameter declarations */ public int getNumPerDonorParameters(){return OpenMMNative.OpenMM_CustomHbondForce_getNumPerDonorParameters(getPointer());}
  /** @return per-acceptor parameter declarations */ public int getNumPerAcceptorParameters(){return OpenMMNative.OpenMM_CustomHbondForce_getNumPerAcceptorParameters(getPointer());}
  /** @return global parameter declarations */ public int getNumGlobalParameters(){return OpenMMNative.OpenMM_CustomHbondForce_getNumGlobalParameters(getPointer());}
  /** @return legacy function count; prefer {@link #getNumTabulatedFunctions()} */
  public int getNumFunctions(){return OpenMMNative.OpenMM_CustomHbondForce_getNumFunctions(getPointer());}
  /** @return registered tabulated-function count */ public int getNumTabulatedFunctions(){return OpenMMNative.OpenMM_CustomHbondForce_getNumTabulatedFunctions(getPointer());}
  /** @return enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic} */
  public int getNonbondedMethod(){return OpenMMNative.OpenMM_CustomHbondForce_getNonbondedMethod(getPointer());}
  /** Set enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic}.
   * @param method native enum value */
  public void setNonbondedMethod(int method){OpenMMNative.OpenMM_CustomHbondForce_setNonbondedMethod(getPointer(),method);}
  /** @return cutoff distance in nm */ public double getCutoffDistance(){return OpenMMNative.OpenMM_CustomHbondForce_getCutoffDistance(getPointer());}
  /** Set cutoff distance in nanometers.
   * @param distance cutoff in nm */
  public void setCutoffDistance(double distance){OpenMMNative.OpenMM_CustomHbondForce_setCutoffDistance(getPointer(),distance);}

  /** Replace the expression; existing contexts must be recreated to use it.
   * @param energy custom donor/acceptor expression; geometry functions use nm/radians */
  public void setEnergyFunction(String energy){withString(energy,n->OpenMMNative.OpenMM_CustomHbondForce_setEnergyFunction(getPointer(),n));}
  /** Replace expression from a caller-owned NUL-terminated UTF-8 segment.
   * @param energy native expression segment borrowed during the call */
  public void setEnergyFunction(MemorySegment energy){OpenMMNative.OpenMM_CustomHbondForce_setEnergyFunction(getPointer(),energy);}
  /** Rename a global declaration; live contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new name */
  public void setGlobalParameterName(int index,String name){withString(name,n->OpenMMNative.OpenMM_CustomHbondForce_setGlobalParameterName(getPointer(),index,n));}
  /** Rename a global declaration using a caller-owned NUL-terminated UTF-8 segment.
   * @param index global-parameter index
   * @param name segment borrowed during the call */
  public void setGlobalParameterName(int index,MemorySegment name){OpenMMNative.OpenMM_CustomHbondForce_setGlobalParameterName(getPointer(),index,name);}
  /** Set the default for future contexts, not the live value in existing contexts.
   * @param index global-parameter index
   * @param defaultValue new default */
  public void setGlobalParameterDefaultValue(int index,double value){OpenMMNative.OpenMM_CustomHbondForce_setGlobalParameterDefaultValue(getPointer(),index,value);}
  /** Rename a per-donor declaration; existing contexts are not reconfigured.
   * @param index declaration index
   * @param name new name */
  public void setPerDonorParameterName(int index,String name){withString(name,n->OpenMMNative.OpenMM_CustomHbondForce_setPerDonorParameterName(getPointer(),index,n));}
  /** Rename from a caller-owned NUL-terminated UTF-8 segment.
   * @param index declaration index
   * @param name segment borrowed during the call */
  public void setPerDonorParameterName(int index,MemorySegment name){OpenMMNative.OpenMM_CustomHbondForce_setPerDonorParameterName(getPointer(),index,name);}
  /** Rename a per-acceptor declaration; existing contexts are not reconfigured.
   * @param index declaration index
   * @param name new name */
  public void setPerAcceptorParameterName(int index,String name){withString(name,n->OpenMMNative.OpenMM_CustomHbondForce_setPerAcceptorParameterName(getPointer(),index,n));}
  /** Rename from a caller-owned NUL-terminated UTF-8 segment.
   * @param index declaration index
   * @param name segment borrowed during the call */
  public void setPerAcceptorParameterName(int index,MemorySegment name){OpenMMNative.OpenMM_CustomHbondForce_setPerAcceptorParameterName(getPointer(),index,name);}
  /** Replace donor slots and role-specific values; Java array is copied.
   * @param index donor-list index
   * @param d1 first particle slot
   * @param d2 second particle slot
   * @param d3 third particle slot
   * @param parameters one per-donor value per declaration, in order */
  public void setDonorParameters(int index,int d1,int d2,int d3,double[] parameters){
    try(DoubleArray values=CustomForceParameters.toNative(parameters)){setDonorParameters(index,d1,d2,d3,values);}
  }
  /** Replace donor values using a native array borrowed for the call.
   * @param index donor-list index
   * @param d1 first particle slot
   * @param d2 second particle slot
   * @param d3 third particle slot
   * @param parameters values in per-donor declaration order */
  public void setDonorParameters(int index,int d1,int d2,int d3,DoubleArray parameters){
    OpenMMNative.OpenMM_CustomHbondForce_setDonorParameters(getPointer(),index,d1,d2,d3,parameters.getPointer());
  }
  /** Replace donor values using a caller-owned native segment borrowed for the call.
   * @param index donor-list index
   * @param d1 first particle slot
   * @param d2 second particle slot
   * @param d3 third particle slot
   * @param parameters native values in declaration order */
  public void setDonorParameters(int index,int d1,int d2,int d3,MemorySegment parameters){
    OpenMMNative.OpenMM_CustomHbondForce_setDonorParameters(getPointer(),index,d1,d2,d3,parameters);
  }
  /** Replace acceptor slots and role-specific values; Java array is copied.
   * @param index acceptor-list index
   * @param a1 first particle slot
   * @param a2 second particle slot
   * @param a3 third particle slot
   * @param parameters one per-acceptor value per declaration, in order */
  public void setAcceptorParameters(int index,int a1,int a2,int a3,double[] parameters){
    try(DoubleArray values=CustomForceParameters.toNative(parameters)){setAcceptorParameters(index,a1,a2,a3,values);}
  }
  /** Replace acceptor values using a native array borrowed during the call.
   * @param index acceptor-list index
   * @param a1 first particle slot
   * @param a2 second particle slot
   * @param a3 third particle slot
   * @param parameters values in declaration order */
  public void setAcceptorParameters(int index,int a1,int a2,int a3,DoubleArray parameters){
    OpenMMNative.OpenMM_CustomHbondForce_setAcceptorParameters(getPointer(),index,a1,a2,a3,parameters.getPointer());
  }
  /** Replace acceptor values using a caller-owned native segment borrowed during the call.
   * @param index acceptor-list index
   * @param a1 first particle slot
   * @param a2 second particle slot
   * @param a3 third particle slot
   * @param parameters native values in declaration order */
  public void setAcceptorParameters(int index,int a1,int a2,int a3,MemorySegment parameters){
    OpenMMNative.OpenMM_CustomHbondForce_setAcceptorParameters(getPointer(),index,a1,a2,a3,parameters);
  }
  /** Replace an exclusion pair.
   * @param index exclusion index
   * @param donor donor-list index
   * @param acceptor acceptor-list index */
  public void setExclusionParticles(int index,int donor,int acceptor){OpenMMNative.OpenMM_CustomHbondForce_setExclusionParticles(getPointer(),index,donor,acceptor);}
  /** Replace a legacy continuous 1-D function, deprecated by the native API; prefer the
   * tabulated-function parameter API.
   * @param index function index
   * @param name function name
   * @param values sampled values
   * @param min lower input bound
   * @param max upper input bound */
  public void setFunctionParameters(int index,String name,double[] values,double min,double max){
    try(DoubleArray nativeValues=CustomForceParameters.toNative(values)){
      withString(name,n->OpenMMNative.OpenMM_CustomHbondForce_setFunctionParameters(getPointer(),index,n,nativeValues.getPointer(),min,max));
    }
  }
  /** Replace a legacy function from caller-owned UTF-8 name and value-array segments borrowed
   * during the call.
   * @param index function index
   * @param name NUL-terminated UTF-8 function name
   * @param values native sampled values
   * @param min lower input bound
   * @param max upper input bound */
  public void setFunctionParameters(int index,MemorySegment name,MemorySegment values,double min,double max){
    OpenMMNative.OpenMM_CustomHbondForce_setFunctionParameters(getPointer(),index,name,values,min,max);
  }
  /** Apply only per-donor/per-acceptor values and tabulated-function values to a context. Donor or
   * acceptor particle slots cannot change and new groups cannot be added; function dimensions and
   * domain/range must remain unchanged.
   * @param context context containing this force; no-op without a native context pointer. Does
   *     not update expression, nonbonded method, cutoff, exclusions, declarations, or defaults. */
  public void updateParametersInContext(Context context){CustomForceParameters.updateContext(context,c->OpenMMNative.OpenMM_CustomHbondForce_updateParametersInContext(getPointer(),c));}
  /** @return whether periodic boundary conditions are enabled. */
  @Override public boolean usesPeriodicBoundaryConditions(){return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomHbondForce_usesPeriodicBoundaryConditions(getPointer()));}

  private GroupParameters getGroupParameters(int index,boolean donor){
    try(Arena arena=Arena.ofConfined();DoubleArray values=new DoubleArray(0)){
      MemorySegment a=arena.allocate(ValueLayout.JAVA_INT),b=arena.allocate(ValueLayout.JAVA_INT),c=arena.allocate(ValueLayout.JAVA_INT);
      if(donor)OpenMMNative.OpenMM_CustomHbondForce_getDonorParameters(getPointer(),index,a,b,c,values.getPointer());
      else OpenMMNative.OpenMM_CustomHbondForce_getAcceptorParameters(getPointer(),index,a,b,c,values.getPointer());
      return new GroupParameters(a.get(ValueLayout.JAVA_INT,0),b.get(ValueLayout.JAVA_INT,0),c.get(ValueLayout.JAVA_INT,0),CustomForceParameters.copy(values));
    }
  }
  private static MemorySegment create(String energy){return withStringResult(energy,n->{OpenMMRuntime.initialize();return OpenMMNative.OpenMM_CustomHbondForce_create(n);});}
  private static void withString(String value,java.util.function.Consumer<MemorySegment> action){OpenMMStrings.withUtf8String(value,action);}
  private static <T>T withStringResult(String value,java.util.function.Function<MemorySegment,T> action){return OpenMMStrings.withUtf8StringResult(value,action);}
  private static UnsupportedOperationException invalidStringOutput(String method){return new UnsupportedOperationException("OpenMM_CustomHbondForce_"+method+" has an incompatible char** output ABI in this OpenMM wrapper.");}
}
