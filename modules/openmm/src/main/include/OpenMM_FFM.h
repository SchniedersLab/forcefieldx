#ifndef FFX_OPENMM_FFM_H
#define FFX_OPENMM_FFM_H

/* Avoid a jextract-generated catch-variable collision in Gay-Berne APIs. */
#define ex epsilonX
#include "OpenMMCWrapper.h"
#undef ex
#include "AmoebaOpenMMCWrapper.h"
#include "DrudeOpenMMCWrapper.h"

#endif
