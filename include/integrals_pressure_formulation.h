#include "vector.h"

#ifndef pi
#define pi M_PI
#endif

/* Integrands and singularities for the 'general' wave formulation. */
void compute_convolution_s(const int,      const vector &,
                           const vector &, const vector &,
                           const vector &,       vector &);

void compute_far_downstream_integral_s(const vector &, const real,
                                       const real, const real, vector &);

void compute_far_upstream_integral_s(const vector &, const real,
                                     const real, const real, vector &);

void removed_singularity_general_s(const vector&, const vector&, const vector&, vector&);

/* Integrands and singularities for the 'symmetric' wave formulation. */
void compute_convolution_symmetric_s(const int,      const vector &,
                                     const vector &, const vector &,
                                     const vector &,       vector &);

void removed_singularity_symmetric_s(const vector&, const vector&, const vector&, vector&);
