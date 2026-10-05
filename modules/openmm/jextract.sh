#!/usr/bin/env bash
#
# Generate Java Foreign Function & Memory API bindings for the OpenMM C ABI.
#
# Requirements:
#   - jextract 25 (override with JEXTRACT=/path/to/jextract)
#   - OpenMM C headers matching the target native libraries (override with
#     OPENMM_INSTALL=/path/to/openmm-install)
#
# The generated files are written to modules/openmm/src/main/java. They
# are intended to be committed; ordinary Maven builds do not run this script.

set -euo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly REPOSITORY_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
readonly DEFAULT_JEXTRACT="/Library/Java/JavaVirtualMachines/jextract-25/bin/jextract"
readonly JEXTRACT="${JEXTRACT:-${DEFAULT_JEXTRACT}}"
readonly OPENMM_INSTALL="${OPENMM_INSTALL:-/Users/mjschnie/Data/ffx-project/forcefieldx/openmm/openmm-install}"
readonly JAVA_SOURCE_DIRECTORY="${SCRIPT_DIR}/src/main/java"
readonly GENERATED_DIRECTORY="${JAVA_SOURCE_DIRECTORY}/ffx/openmm/ffm/bindings"
readonly FFM_HEADER="${SCRIPT_DIR}/src/main/include/OpenMM_FFM.h"
readonly PACKAGE="ffx.openmm.ffm.bindings"

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

[[ -x "${JEXTRACT}" ]] || die "jextract 25 is not executable: ${JEXTRACT}"
[[ -d "${OPENMM_INSTALL}" ]] || die "OpenMM installation does not exist: ${OPENMM_INSTALL}"

JEXTRACT_VERSION="$("${JEXTRACT}" --version 2>&1)"
[[ "${JEXTRACT_VERSION}" == jextract\ 25* ]] ||
  die "jextract 25 is required; found: ${JEXTRACT_VERSION}"

readonly INCLUDE_DIRECTORY="${OPENMM_INSTALL}/include"
readonly CORE_HEADER="${INCLUDE_DIRECTORY}/OpenMMCWrapper.h"
readonly AMOEBA_HEADER="${INCLUDE_DIRECTORY}/AmoebaOpenMMCWrapper.h"
readonly DRUDE_HEADER="${INCLUDE_DIRECTORY}/DrudeOpenMMCWrapper.h"

for required_path in \
  "${FFM_HEADER}" \
  "${CORE_HEADER}" \
  "${AMOEBA_HEADER}" \
  "${DRUDE_HEADER}"; do
  [[ -f "${required_path}" ]] || die "required OpenMM file does not exist: ${required_path}"
done

# Clear only the generated binding package, never handwritten Java sources.
rm -rf "${GENERATED_DIRECTORY}"
mkdir -p "${JAVA_SOURCE_DIRECTORY}"

"${JEXTRACT}" \
  --output "${JAVA_SOURCE_DIRECTORY}" \
  --target-package "${PACKAGE}" \
  --header-class-name OpenMMNative \
  --include-dir "${INCLUDE_DIRECTORY}" \
  "${FFM_HEADER}"
