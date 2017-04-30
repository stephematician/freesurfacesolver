/** \file integrals_sbox.h 
  
    Common functions used by all methods of calculating flows with box shape
    disturbances. 
  
*/
#ifndef INTEGRALS_SBOX_H
#define INTEGRALS_SBOX_H

#include <stdlib.h> // For exit() and EXIT_FAILURE

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "pwpoly.h" /* transform_func */
#include "integrals_types.h"
#include <gsl/gsl_integration.h>

#ifndef pi
#define pi M_PI /**< Constant value of \f$\pi\f$. */
#endif

typedef enum {
    USE_COSINE,
    USE_SINE
} cosine_or_sine;


/** \internal Approximation to a contribution from farstream to a Cauchy P.V.
              integral. */
typedef real (*farstream_integrand_func)(real, void*);


/** \internal Expression for a removed singularity from a Cauchy P.V. integral.
*/ 
typedef void (*removed_singularity_func)(const real_vector &,
                                         const real_vector &,
                                         const real_vector &,
                                         real_vector &);


/** \internal Function to calculate \f$\theta_b(\phi)\f$ given topography
              parameters. */
typedef real (*tpgraphy_theta_func)(const real, const real_vector&);


/** \internal Removal of singularity on topography bottom. */
typedef void (*removed_singularity_b_func)(const real, const real,
                                           const real_vector &, real &);


/** \internal Expression for the integrand for the topographical term used to
              calculate tau on the free surface. */
typedef real (*tpgraphy_integrand_func)(real, void*);


/** \internal Expression for the integrand of the Cauchy P.V. integral. */ 
typedef void (*tpgraphy_convolution_func)(const real,
                                          const real_vector &,
                                          const real_vector &,
                                          const real_vector &,
                                          real_vector &);

/** Functions required by the Cauchy P.V. boundary integrals for free surface
    calculations with topographical disturbances. 
    
    \param transform_phi Transform (arbitrary) relating \f$\beta\f$ to
                         \f$\phi\f$.
    \param transform_dphi Transform (arbitrary) relating \f$\mathrm{d} \beta\f$
                          to \f$\mathrm{d} \phi\f$.
    \param convolution_terms Expression for the integrand of the Cauchy P.V.
                             integral.
    \param removed_singularity Optional term for the removed singularity from 
                               the Cauchy P.V. integral.
    \param linear_upstream_int Optional approximation to the contribution from
                               far upstream.
    \param linear_downstream_int Optional approximation to the contribution from
                                 far downstream.
    \param topography_integrand Expression for the integrand of the channel
                                topography.
*/
typedef struct {
    transform_func           transform_phi;
    transform_func           transform_dphi;
    convolution_func         convolution_terms;
    removed_singularity_func removed_singularity;
    farstream_integrand_func linear_upstream_int;
    farstream_integrand_func linear_downstream_int;
    tpgraphy_integrand_func  topography_integrand;
} method_s_funcs;


/** Functions required by the Cauchy P.V. boundary integrals for geometry of
    topography calculations.

    \param transform_phi Transform (arbitrary) relating \f$\beta\f$ to 
                           \f$\phi\f$.
    \param transform_dphi Transform (arbitrary) relating \f$\mathrm{d} \beta\f$
                          to \f$\mathrm{d} \phi\f$.
    \param theta_func The type of topography being used expressed using
                      \f$\theta\f$.
    \param convolution_s Expression for the integrand of the integral over the
                         free-surface (non P.V.).
    \param topography_integrand Expression for the integrand of the integral 
                                over the topography. (P.V).
    \param removed_singularity_b Optional term for the removed singularity from
                                 the Cauchy P.V. integral over topography.
    \param linear_upstream_s Optional approximation to the contribution from
                             far upstream.
    \param linear_downstream_s Optional approximation to the contribution from
                               far downstream.
*/
typedef struct {
    transform_func             transform_phi;
    transform_func             transform_dphi;
    tpgraphy_theta_func        theta_func;
    tpgraphy_convolution_func  convolution_s;
    tpgraphy_integrand_func    topography_integrand;
    removed_singularity_b_func removed_singularity_b;
    farstream_integrand_func   linear_upstream_s;
    farstream_integrand_func   linear_downstream_s;
} method_b_funcs;


/// \cond EXCLUDE_SYMBOLS
typedef struct {
    method_b_funcs method_b;
    real_vector theta_s;
    real_vector beta_s;
    real_vector beta_sub;
    real_vector phi_sub;
    real_vector dphi_sub;
    ff_params ffus_params;
    ff_params ffds_params;
    real homotopy_s;
    real_vector topography;
    cosine_or_sine cos_sin;
    gsl_integration_workspace** ws;
} z_b_integrand_params;
/// \endcond


/** Generic entry point for Lua loader

    Registers the library with Lua, creates a GSL integration workspace for use
    by the boundary integral method functions, and sets up the garbage
    collection for the workspace.
 */
int luaopen_integrals_sbox(lua_State *L);


/** \todo document error handler */
void integrals_sbox_handler(const char *reason, 
                            const char   *file, 
                            int           line, 
                            int      gsl_errno);


/** Generalised trapezoidal method for calculating \f$\tau_s(\beta_{i+1/2})\f$
    based on Cauchy P.V. boundary integral methods.

    \param theta_s (n)-real_vector the angle of the free surface at
                   \f$\beta_{i}\f$.
    \param phi_sub (m)-real_vector of 'clustering' locations in computational
                   grid; \f$\phi^{[j]}\f$.
    \param dphi_sub (m)-real_vector of prescribed value of \f$\phi'\f$ at
                    \f$\phi^{[j]}\f$.
    \param beta (n)-real_vector of \f$\beta_{i}\f$.
    \param beta_sub (m)-real_vector of \f$\beta^{[j]}\f$.
    \param topography parameters with mathematically describe topography.
    \param ffus_params parameters for linearised region upstream.
    \param ffds_params parameters for linearised region downstream.
    \param s optional homotopy parameter.
    \param ws integration workspace for linearised regions.
    \param method current method functions.
    \param tau_s (n-1)-vector to store returned \f$\tau_s(\beta_{i+1/2})\f$.
 */
void tau_s_worker(const real_vector &    theta_s, 
                  const real_vector &    phi_sub,
                  const real_vector &   dphi_sub,
                  const real_vector &       beta,
                  const real_vector &   beta_sub,
                  const real_vector & topography,
                  const ff_params &  ffus_params,
                  const ff_params &  ffds_params,
                  const real                   s,
                  gsl_integration_workspace** ws,
                  const method_s_funcs &  method,
                  real_vector &            tau_s);


/** Generalised mid-point method for calculating \f$x_s(\beta_{i+1/2})\f$ r
    \f$y_s(\beta_{i+1/2})\f$ based on \f$\theta\f$ and \f$\tau\f$.

  Uses the mid-point integration rule to evaluate \f$x(\beta_{j+1/2})\f$, by
  averaging the value at \f$j\f$ and \f$j+1\f$. Whether to calculate \f$x_s\f$
  or \f$y_s\f$ is determine by using the cosine or sine function. 
 
  For \f$x_s\f$ let
  \f[ W^{j} = \exp(-\tau_s(\beta_j)) \cos(\theta_s(\beta_j)) \f]
  otherwise for \f$y_s\f$
  \f[ W^{j} = \exp(-\tau_s(\beta_j)) \cos(\theta_s(\beta_j)) \f]
  then the integral is given by
  \f[ I_s(\beta_{j-1/2}) = I_s(\beta_{j+1/2}) - \frac{1}{2}
                           \left( (\beta_{j+1}-\beta_{j}) W^{j+1/2} +
                                  (\beta_{j}-\beta_{j-1}) W^{j-1/2} \right). \f]

    \param theta_mid \f$(n-1)\f$-real_vector; \f$\theta\f$ on free surface at
                     \f$\beta_{j+1/2}\f$.
    \param tau_mid \f$(n-1)\f$-real_vector; \f$\tau\f$ on free surface at
                   \f$\beta_{j+1/2}\f$.
    \param beta \f$(n)\f$-real_vector; \f$\beta_{j}\f$.
    \param beta_sub \f$(m)\f$-real_vector; value of \f$\beta\f$ at
                    \f$\phi^{[j]}\f$.
    \param phi_sub \f$(m)\f$-real_vector; 'clustering' locations in
                   computational grid; \f$\phi^{[j]}\f$.
    \param dphi_sub \f$(m)\f$-real_vector; value of \f$\phi'\f$ at
                    \f$\phi^{[j]}\f$.
    \param s optional homotopy parameter.
    \param transform_dphi function to calculate \f$\phi'\f$ from
                          \f$\beta\f$.
    \param cos_sin use cosine (calculate \f$x_s\$f) or sine (calculate
                   \f$y_s\f$).
    \param x0 return value of \f$x\f$ or \f$y\f$ at end-point of integral
              \f$\beta^{[0]}\f$.
    \param z_s (n-1)-vector to store returned \f$x_s(\beta_{j+1/2})\f$ or
               \f$y_s(\beta_{j+1/2})\f$.

*/
void z_s_worker(const real_vector&  theta_mid,
                const real_vector&    tau_mid, 
                const real_vector&       beta, 
                const real_vector&   beta_sub, 
                const real_vector&    phi_sub, 
                const real_vector&   dphi_sub,
                const real                  s,
                transform_func transform_dphi, 
                const cosine_or_sine  cos_sin,
                real&                      z0,
                real_vector&              z_s);


/** \todo Document b_worker */
void z_b_worker(const real_vector&           phi,
                const real_vector&       theta_s,
                const real_vector&        beta_s,
                const real_vector&      beta_sub,
                const real_vector&       phi_sub,
                const real_vector&      dphi_sub,
                const real                     s,
                gsl_integration_workspace**   ws,
                method_b_funcs            method,
                const real_vector&    topography,
                const ff_params&     ffus_params,
                const ff_params&     ffds_params,
                const cosine_or_sine     cos_sin,
                real_vector&                 z_b);

/** \todo document z_b_integrand */
real z_b_integrand(real phi_b, void* _params);

#endif

