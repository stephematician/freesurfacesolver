/** \file integrals_sbox_formulation.h

  C implementation of integrands (such as convolutions and the farstream linearised solution)
  and singularities for smooth topographical disturbances described by a function \f$f(\phi)\f$,
  such that
  \f[ y_b(x) \approx f(\phi), \quad \theta_b(x) \approx \arctan \frac{\mathrm{d} y_b}{\mathrm{d} \phi} \f]

  Two cases have been implemented, a smoothed box where \f$ f(\phi) \sim \mathrm{arctanh} \f$, and a
  Gaussian shape where \f$f(\phi) \sim \exp^{-\phi^2} \f$.

*/
#include "vector.h"

#ifndef pi
#define pi M_PI
#endif

/** Evaluates \f$\theta_b(\phi)\f$ for the gaussian topography.

  \f[ \theta_b(\phi) = \arctan \frac{\mathrm{d}}{\mathrm{d} \phi} \frac{AB}{\sqrt{\pi}} e^{-(B\phi)^2}. \f]
  
 */
real gaussian_theta(const real phi, /**< \f$\phi\f$ */
		    const vector& topography /**< vector containing \f$A\f$ and \f$B\f$. */);

/** Evaluates \f$\theta_b(\phi)\f$ for the smoothed box topography.

  \f[ \theta_b(\phi) = \arctan \frac{\mathrm{d}}{\mathrm{d} \phi} \frac{A}{2\phi_c}
                       \left( \mathrm{tanh} B(\phi+\phi_c/2) - \mathrm{tanh} B(\phi-\phi_c/2) \right). \f]
  
 */
real smoothbox_theta(const real phi, /**< \f$\phi\f$ */
		     const vector& topography /**< vector containing \f$A, B\f$ and \f$\phi_c\f$. */);

/* Integrands and singularities for the 'general' wave formulation. */



/** Computes the terms for convolution sum/integral of \f$\theta_s\f$ for evaluating \f$\tau_s\f$.

  \todo Detailed description.
 */
void compute_convolution_s(const int           i, /**< index of the current convolution */
                           const vector& theta_s, /**< (n)-vector  the angle of the free surface at \f$\beta_{i}\f$. */
                           const vector&     phi, /**< (n)-vector  the value of \f$\phi_{i}\f$. */
                           const vector&    dphi, /**< (n)-vector  the value of \f$\phi_{i}'\f$. */
                           const vector& phi_mid, /**< (n-1)-vector  the value of \f$\phi_{i+1/2}\f$. */
                           vector&             o  /**< (n)-vector  the values for the convolution sum. */);

/** Computes the terms for convolution sum/integral of \f$\theta_s\f$ for evaluating \f$\tau_s\f$ for symmetric problems.

  \todo Detailed description.
 */
void compute_convolution_symmetric_s(const int           i, /**< index of the current convolution */
                                     const vector& theta_s, /**< (n)-vector  the angle of the free surface at \f$\beta_{i}\f$. */
                                     const vector&     phi, /**< (n)-vector  the value of \f$\phi_{i}\f$. */
                                     const vector&    dphi, /**< (n)-vector  the value of \f$\phi_{i}'\f$. */
                                     const vector& phi_mid, /**< (n-1)-vector  the value of \f$\phi_{i+1/2}\f$. */
                                     vector&             o  /**< (n)-vector  the values for the convolution sum. */);

/** Computes the integrand for the far upstream contribution to \f$\tau_s\f$ using linearised solution.

  The integrand
  \f[ \int_{0}^{-e^{-\pi\phi_m}} \frac{1}{\pi}
      \frac{\arctan(v_s(\alpha)/u_s(\alpha))}{\hat{\alpha} + \alpha} \, \mathrm{d} \hat{\alpha}. \f]
  is computed using the approximation that for very small \f$\hat{\alpha}=-e^{\pi\hat{\phi}}\f$ the velocity obeys
  \f[ u_s(\hat{\alpha}) = 1 + D\cos(\lambda) \hat{\alpha}^{\lambda/\pi}, \f] and
  \f[ v_s(\hat{\alpha}) = -D\sin(\lambda) \hat{\alpha}^{\lambda/\pi}. \f]

 */
real far_upstream_integrand_s(real x, /**< The value of \f$x \in [0,e^{-\pi\phi_M}]\f$.*/
                              void *params /**<  A real array containing \f$D,\lambda\f$ and \f$\phi\f$.*/);

/** Computes the integrand for the far downstream contribution to \f$\tau_s\f$ using linearised solution.

  The integrand
  \f[ \int_{e^{\pi\phi_m}}^{\infty} \frac{1}{\pi}
      \frac{\arctan(v_s(\alpha)/u_s(\alpha))}{\hat{\alpha} + \alpha} \, \mathrm{d} \hat{\alpha}. \f]
  is computed using the approximation that for very large \f$\hat{\alpha}=-e^{\pi\hat{\phi}}\f$ the velocity obeys
 \f[ u_s(\hat{\alpha}) = 1 + D\cos(\lambda) \hat{\alpha}^{-\lambda/\pi}, \f] and
  \f[ v_s(\hat{\alpha}) = D\sin(\lambda) \hat{\alpha}^{-\lambda/\pi}. \f]

  The transform
  \f$ \hat{\alpha} = e^{\pi\phi_m} + \frac{1 - x^{\gamma}}{x^{\gamma}} \f$ is used in order to map
  the infinite integral to a finite interval.
 */
real far_downstream_integrand_s(real x,      /**< The value of \f$x \in [0,1]\f$. */
				void *params /**< A real array containing \f$D,\lambda,\phi,\phi_m\f$ and \f$\gamma\f$ */);


/** Computes the integrand for the far upstream integral using linearised solution for symmetric problem.

  Actually computes the same integrand as \ref far_upstream_integrand_s.
 */
real far_upstream_symmetric_integrand_s(real, void *);

/** Computes the integrand for the far downstream integral using linearised solution for symmetric problem.

  Actually computes the same integrand as \ref far_downstream_integrand_s.
 */
real far_downstream_symmetric_integrand_s(real, void *);


/** Removal of singularity for the integral over the free surface for \f$\tau_s\f$.

  Calculates all the removed singularities for each mid-point calculation performed;
  \f[ \theta_s(\phi_{i+1/2}) \ln \left| 
     \frac{e^{\pi\phi_b} - e^{\pi\phi_{i+1/2}}}{e^{\pi\phi_a} - e^{\pi\phi_{i+1/2}}} \right|, \f]
  where \f$\phi_a\f$ and \f$\phi_b\f$ are the left and right hand end points of the truncated
  domain.

 */
void removed_singularity_general_s(const vector& theta_s, /**< (n)-vector  the angle of the free surface \f$\theta_s\f$ at \f$\phi_{i+1/2}\f$. */
                                   const vector&     phi, /**< (n)-vector  the values of \f$\phi_i\f$. */
                                   const vector& phi_mid, /**< (n-1)-vector  the values of \f$\phi_{i+1/2}\f$. */
                                   vector&             o  /**< (n-1)-vector  the values of the singularity. */);

/** Removal of singularity for the integral over the free surface for \f$\tau_s\f$ for symmetric problems.

  Calculates all the removed singularities for each mid-point calculation performed;
  \f[ \theta_s(\phi_{i+1/2}) \ln \left| 
     \frac{\mathrm{sinh} (\pi(\phi_b - \phi_{i+1/2})/2)}{\mathrm{sinh} (\pi(\phi_b + \phi_{i+1/2})/2)} \right|, \f]
  where \f$\phi_b\f$ is right hand end point of the truncated (symmetric) domain.

 */
void removed_singularity_symmetric_s(const vector& theta_s, /**< (n)-vector  the angle of the free surface \f$\theta_s\f$ at \f$\phi_{i+1/2}\f$. */
                                     const vector&     phi, /**< (n)-vector  the values of \f$\phi_i\f$. */
                                     const vector& phi_mid, /**< (n-1)-vector  the values of \f$\phi_{i+1/2}\f$. */
                                     vector&             o  /**< (n-1)-vector  the values of the singularity. */);

/** Integrand for convolution over topography for computing \f$\tau_s\f$.

  Integrand of
  \f[ \int_{-\infty}^{\infty} \frac{\theta_b(\hat{\phi})}{e^{\pi\hat{\phi}} + e^{\pi\phi}} \, \mathrm{d} \hat{\phi}. \f]
  where \f$\theta_b\f$ is given by the Gaussian function with parameters \f$A\f$ and \f$B\f$.

  Uses a transform from a semi-infinite domain to the finite domain \f$(-1,1)\f$,
  \f[ \hat\phi = B \sqrt{2} \frac{\tanh^{-1}(\beta)}{\tanh^{-1}(1/2)}. \f]
 */
real gaussian_b(real x,      /**< Value of \f$x\f$ in \f$[0,1]\f$. */
                void *params /**< Real array containing \f$A, B\f$ and \f$\phi\f$.*/);

/** Integrand for convolution over topography for computing \f$\tau_s\f$
 
  Calculates the integrand
  \f[ \int_{-\infty}^{\infty} \frac{\theta_b(\phi)}{\exp^{\hat{\phi}} + \alpha} e^{\pi\hat{\phi}} \, \mathrm{d} \hat{\phi}. \f]

  Uses a transform from a infinite domain \f$\mathbb{R}\f$ to the finite domain \f$(-1,1)\f$.
  \f[ \hat\phi = \phi_c \frac{\tanh^{-1}(\beta)}{\tanh^{-1}(1/2)}. \f]
 */
real smoothbox_b(real x, /**< Value of \f$x\f$ in \f$(-1,1)\f$ */
                 void *params /**<  Real array containing \f$A, B, \phi_c\f$ and \f$\phi\f$.*/);

/* NEW STUFF!

 What I need
 1. OVER FREE SURFACE:
    * Far upstream - done
    * Far downstream - done
    - Integrand in middle part - done
 3. OVER BOTTOM TOPOGRAPHY:
    * functions as per matlab tests that computes the integral over the bottom. - done
 */

/** Computes the terms for convolution sum/integral of \f$\theta_s\f$ for evaluating \f$\tau_b\f$.

  \f[ \int_{\phi_a}^{\phi_b} \frac{\theta_s(\hat{\phi})}{e^{\pi\hat{\phi}} + e^{\pi\phi}} e^{\pi\hat{\phi}} 
      \, \mathrm{d} \hat{\phi}. \f]

  \todo Detailed description.
 */
void topography_compute_convolution_s(const real      phi_b, /**<value of \f$\phi\f$ (on topography) in convolution. */
                                      const vector& theta_s, /**< (n)-vector  the angle of the free surface at \f$\beta_{i}\f$. */
                                      const vector&     phi, /**< (n)-vector  the value of \f$\hat{\phi}_{i}\f$. */
                                      const vector&    dphi, /**< (n)-vector  the value of \f$\hat{\phi}_{i}'\f$. */
                                      vector&             o  /**< (n)-vector  the values for the convolution sum.*/);

/** Computes the terms for convolution sum/integral of \f$\theta_s\f$ for evaluating \f$\tau_b\f$.

  \f[ \int_{\phi_a}^{\phi_b} \frac{\theta_s(\hat{\phi})}{e^{\pi\hat{\phi}} + e^{\pi\phi}} e^{\pi\hat{\phi}} 
      \, \mathrm{d} \hat{\phi}. \f]

  \todo Detailed description.
 */
void topography_compute_symmetric_convolution_s(const real      phi_b, /**<value of \f$\phi\f$ (on topography) in convolution. */
                                                const vector& theta_s, /**< (n)-vector  the angle of the free surface at \f$\beta_{i}\f$. */
                                                const vector&     phi, /**< (n)-vector  the value of \f$\hat{\phi}_{i}\f$. */
                                                const vector&    dphi, /**< (n)-vector  the value of \f$\hat{\phi}_{i}'\f$. */
                                                vector&             o  /**< (n)-vector  the values for the convolution sum.*/);
/** Computes the integrand for the far upstream contribution to \f$\tau_b\f$ using a linearised solution.

  The integrand
  \f[ \int_{0}^{-e^{\pi\phi_a}} \frac{1}{\pi}
      \frac{\arctan(v_s(\alpha)/u_s(\alpha))}{\hat{\alpha} + \alpha} \, \mathrm{d} \hat{\alpha}. \f]
  is computed using the approximation that for very small \f$\hat{\alpha}=e^{\pi\hat{\phi}}\f$ the velocity obeys
  \f[ u_s(\hat{\alpha}) = 1 + D\cos(\lambda) \hat{\alpha}^{\lambda/\pi}, \f] and
  \f[ v_s(\hat{\alpha}) = -D\sin(\lambda) \hat{\alpha}^{\lambda/\pi}. \f]

 */
real topography_far_upstream_s(real x,      /**< The value of \f$x \in [0,e^{\pi\phi_a}]\f$.*/
                               void *params /**<  A real array containing \f$D,\lambda\f$ and \f$\phi\f$.*/);

/** Computes the integrand for the far downstream contribution to \f$\tau_b\f$ using a linearised solution.

  The integrand
  \f[ \int_{-e^{\pi\phi_b}}^{-\infty} \frac{1}{\pi}
      \frac{\arctan(v_s(\alpha)/u_s(\alpha))}{\hat{\alpha} + \alpha} \, \mathrm{d} \hat{\alpha}. \f]
  is computed using the approximation that for very large \f$\hat{\alpha}=e^{\pi\hat{\phi}}\f$ the velocity obeys
  \f[ u_s(\hat{\alpha}) = 1 + D\cos(\lambda) \hat{\alpha}^{-\lambda/\pi}, \f] and
  \f[ v_s(\hat{\alpha}) = -D\sin(\lambda) \hat{\alpha}^{-\lambda/\pi}. \f]

  The transform
  \f$ \hat{\alpha} = e^{\pi\phi_b} + \frac{1 - x^{\gamma}}{x^{\gamma}} \f$ is used in order to map
  the infinite integral to a finite interval.
 */
real topography_far_downstream_s(real x,      /**< The value of \f$x \in [0,1]\f$.*/
                                 void *params /**< A real array containing \f$D,\lambda,\phi,\phi_b\f$ and \f$\gamma\f$.*/);


/** Removal of singularity for the integral over the bottom for \f$\tau_b\f$.

  \todo detailed description?

 */
void topography_smoothbox_removed_singularity_b(const real           phi, /**< value of \f$\phi\f$. */
                                                const real         phi_m, /**< limit of domain of integration, \f$\phi_m\f$. */
                                                const vector& topography, /**< vector containing \f$A, B\f$ and \f$\phi_c\f$. */
                                                real&                  o  /**< value of removed singularity. */);

/** Integrand of convolution of \f$\theta_b\f$ for evaluating \f$\tau_b\f$.

 */
real topography_smoothbox_integrand_b(real, void *);


/** Removal of singularity for the integral over the bottom for \f$\tau_b\f$.

  \todo detailed description?

 */
void topography_gaussian_removed_singularity_b(const real           phi, /**< value of \f$\phi\f$. */
                                               const real         phi_m, /**< limit of domain of integration, \f$\phi_m\f$. */
                                               const vector& topography, /**< vector containing \f$A, B\f$ and \f$\phi_c\f$. */
                                               real&                  o  /**< value of removed singularity. */);

/** Integrand of convolution of \f$\theta_b\f$ for evaluating \f$\tau_b\f$.

 */
real topography_gaussian_integrand_b(real, void *);

/* Just use the non-symmetric functions for these */
/*
real topography_gaussian_upstream_symmetric_s(real, void *);
real topography_gaussian_downstream_symmetric_s(real, void *);

real topography_smoothbox_upstream_symmetric_s(real, void *);
real topography_smoothbox_downstream_symmetric_s(real, void *);
*/ 


