/** \file integrals_sbox_formulation.h

    C implementation of integrands (such as convolutions and the farstream 
    linearised solution) and singularities for pressure disturbances
    on the free surface.

    The pressure term does not appear in these calculations, and the topography
    is assumed flat, thus there are no functions for calculating integration
    over the topography.
    
    See Chapter 3 of Thesis "Very Steep Forced Solitary Waves in Two-Dimensional
    Free Surface Flow".

*/
#ifndef INTEGRALS_PRESSURE_FORMULATION_H
#define INTEGRALS_PRESSURE_FORMULATION

#include "real_vector.h"

#ifndef pi
#define pi M_PI
#endif


/** Computes the terms for convolution sum/integral of \f$v_s\f$ for
    evaluating \f$u_s\f$.

    \todo Detailed description.

    \param i index of current convolution
    \param v_s (n)-vector  the vertical component of velocity of the free
               surface at \f$\beta_{i}\f$.
    \param phi (n)-vector the value of \f$\phi_{i}\f$.
    \param dphi (n)-vector the value of \f$\phi_{i}'\f$.
    \param phi_mid (n-1)-vector the value of \f$\phi_{i+1/2}\f$.
    \param o (n)-vector the values for the convolution sum.

*/
void compute_convolution_s(const int                 i,
                           const real_vector &     v_s,
                           const real_vector &     phi,
                           const real_vector &    dphi,
                           const real_vector & phi_mid,
                           real_vector &             o);


/** Computes the terms for convolution sum/integral of \f$v_s\f$ for
    evaluating \f$u_s\f$ for symmetric problems.

    \todo Detailed description.

    \param i index of current convolution
    \param v_s (n)-vector  the vertical component of velocity of the free
               surface at \f$\beta_{i}\f$.
    \param phi (n)-vector the value of \f$\phi_{i}\f$.
    \param dphi (n)-vector the value of \f$\phi_{i}'\f$.
    \param phi_mid (n-1)-vector the value of \f$\phi_{i+1/2}\f$.
    \param o (n)-vector the values for the convolution sum.

*/
void compute_convolution_symmetric_s(const int                 i,
                                     const real_vector &     v_s,
                                     const real_vector &     phi,
                                     const real_vector &    dphi,
                                     const real_vector & phi_mid,
                                     real_vector &             o);


/** Computes the integrand for the far downstream contribution to \f$\u_s\f$
    using linearised solution.

    The integrand
    \f[ \int_{e^{\pi\phi_m}}^{\infty} \frac{1}{\pi}
            \frac{v_s(\hat{\alpha})}{\hat{\alpha} + \alpha} \,
            \mathrm{d} \hat{\alpha}. \f]
    is computed using the approximation that for very large
    \f$\hat{\alpha}=-e^{\pi\hat{\phi}}\f$ the velocity obeys
    \f[ v_s(\hat{\alpha}) = D_v\sin(\lambda) \hat{\alpha}^{-\lambda/\pi}. \f]

    \param phi_mid (n-1) vector of \f$\phi(\beta_{i+1/2})$.
    \param phi_match The value of \f$\phi_m\f$.
    \param Dv The real value of \f$D_v\f$ on the downstream linear solution.
    \param lambda The real value of \f$\lambda\f$ on the downstream linear
                  solution.
    \param o (n-1) real vector of the integral evaluated at \f$\beta_{i+1/2}\f$.

*/
void far_downstream_integral_s(const real_vector & phi_mid,
                               const real        phi_match,
                               const real               Dv,
                               const real           lambda,
                               real_vector &             o);

/** Computes the integrand for the far upstream contribution to \f$\u_s\f$
    using linearised solution.

    The integrand
    \f[ \int_{0}^{-e^{-\pi\phi_m}} \frac{1}{\pi}
            \frac{v_s(\hat{\alpha})}{\hat{\alpha} + \alpha}
                \, \mathrm{d} \hat{\alpha}. \f]
    is computed using the approximation that for very small
    \f$\hat{\alpha}=-e^{\pi\hat{\phi}}\f$ the velocity obeys
    \f[ v_s(\hat{\alpha}) = -D_v\sin(\lambda) \hat{\alpha}^{\lambda/\pi}. \f]

    \param phi_mid (n-1) vector of \f$\phi(\beta_{i+1/2})$.
    \param phi_match The value of \f$\phi_m\f$.
    \param Dv The real value of \f$D_v\f$ on the upstream linear solution.
    \param lambda The real value of \f$\lambda\f$ on the upstream linear
                  solution.
    \param o (n-1) real vector of the integral evaluated at \f$\beta_{i+1/2}\f$.

*/
void far_upstream_integral_s(const real_vector & phi_mid,
                             const real        phi_match,
                             const real               Dv,
                             const real           lambda,
                             real_vector &             o);

/** Computes the removed singularity term for numerical convolution of
    \f$v_s\f$.

    See Chapter 3 of thesis.

    \todo document properly

    \param v_s (n) real vector of \f$v(\phi(\beta_{i}))\f$.
    \param phi (n) real vector of \f$\phi(\beta_{i})\f$.
    \param phi_mid (n-1) real vector of \f$\phi(\beta_{i+1/2})\f$.
    \param o (n-1) vector of removed singularity terms at \f$\beta_{i+1/2}\f$.

*/
void removed_singularity_s(const real_vector &     v_s,
                           const real_vector &     phi,
                           const real_vector & phi_mid,
                           real_vector &             o);

/** Computes the removed singularity term for numerical convolution of \f$v_s\f$
    for symmetric problems.

    See Chapter 3 of thesis.

    \todo document properly

    \param v_s (n) real vector of \f$v(\phi(\beta_{i}))\f$.
    \param phi (n) real vector of \f$\phi(\beta_{i})\f$.
    \param phi_mid (n-1) real vector of \f$\phi(\beta_{i+1/2})\f$.
    \param o (n-1) vector of removed singularity terms at \f$\beta_{i+1/2}\f$.

*/
void removed_singularity_symmetric_s(const real_vector &     v_s,
                                     const real_vector &     phi,
                                     const real_vector & phi_mid,
                                     real_vector &             o);

#endif
