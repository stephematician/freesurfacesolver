#ifndef INTEGRAL_TYPES_H
#define INTEGRAL_TYPES_H

#include "real_vector.h"


/** Parameters required by linearised approximations to far-stream
    contributions. */
typedef struct {
    real lambda; /**< Decay rate of velocity. */
    real gamma; /**< Optional transform parameter. */
    real D; /**< Coefficient of solution. */
    real phi_match; /**< Location of farstream matching. */
} ff_params;


/** \internal Expression for the integrand of the Cauchy P.V. integral. */ 
typedef void (*convolution_func)(const int,           const real_vector &,
                                 const real_vector &, const real_vector &,
                                 const real_vector &,       real_vector &);


/** \internal Expression for a removed singularity from a Cauchy P.V. integral.
*/ 
typedef void (*removed_singularity_func)(const real_vector&,
                                         const real_vector&,
                                         const real_vector&,
                                         real_vector&);

#endif

