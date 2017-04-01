/** \file integrals_sbox.h 
  
  Common functions used by all methods of calculating flows with box shape disturbances. 
  
  \todo At some point I will have to consider how to calculate the topography along the bottom.

*/

#include <stdlib.h> // For exit() and EXIT_FAILURE

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "vector.h"
#include <gsl/gsl_integration.h>

#ifndef pi
#define pi M_PI /**< Constant value of \f$\pi\f$. */
#endif

/** Parameters required by linearised approximations to far-stream contributions. */
typedef struct {
  real lambda; /**< Decay rate of velocity. */
  real gamma; /**< Optional transform parameter. */
  real D; /**< Coefficient of solution. */
  real phi_match; /**< Location of farstream matching. */
} linear_params;

/** \internal Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
typedef void (*transform_func)(const vector &, const vector &,
                               const vector &, const vector &,
                                     vector &, const real);

/* For free surface integration */
/** \internal Expression for the integrand of the Cauchy P.V. integral. */ 
typedef void (*convolution_func)(const int,      const vector &,
                                 const vector &, const vector &,
                                 const vector &,       vector &);
/** \internal Approximation to a contribution from farstream to a Cauchy P.V. integral. */
typedef real (*linear_farstream_func)(real, void*);
/** \internal Expression for a removed singularity from a Cauchy P.V. integral. */ 
typedef void (*removed_singularity_func)(const vector&, const vector&, const vector&, vector&);

/* For topography integration */
/** \internal Integrand for x_b*/
// typedef real (*x_b_integrand_func)(real, void*);

/* NEW STUFF!

 What I need
 1. GENERAL:
    - a function that calculates the angle of the topography, theta
 2. OVER FREE SURFACE:
    * Far upstream
    * Far downstream
    - Integrand in middle part
 3. OVER BOTTOM TOPOGRAPHY:
    * function as per matlab tests that computes the integral over the bottom.
 */
/*smoothbox_theta();

void removed_singularity_b(const real,  
                                     const real,
                          const tpgraphy_param&,
				          real&);


real far_upstream_integrand_b(real, void *);
real far_downstream_integrand_b(real, void *);

real topography_smoothbox_integrand_b(real, void *);
*/

/** \internal Function to calculate \f$\theta_b(\phi)\f$ given topography parameters. */
typedef real (*tpgraphy_theta_func)(const real, const vector&);
/** \internal Removal of singularity on topography bottom. */
typedef void (*removed_singularity_b_func)(const real, const real,
                                           const vector &, real &);
/** \internal Expression for the integrand for the topographical term used to calculate tau on the free surface. */
typedef real (*tpgraphy_integrand_func)(real, void*);
/* \internal Contribution from integral over topography far up/down stream */
// typedef real (*tpgraphy_farstream_func)(real, void*);
/** \internal Expression for the integrand of the Cauchy P.V. integral. */ 
typedef void (*tpgraphy_convolution_func)(const real,     const vector &,
                                          const vector &, const vector &,
                                          vector &);

/** Functions required by the Cauchy P.V. boundary integrals for free surface calculations 
    with topographical disturbances. */
typedef struct {
  transform_func transform_phi; /**< Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
  transform_func transform_dphi; /**< Transform (arbitrary) relating \f$\mathrm{d} \beta\f$ to \f$\mathrm{d} \phi\f$. */
  convolution_func convolution_terms; /**< Expression for the integrand of the Cauchy P.V. integral. */
  removed_singularity_func removed_singularity; /**< Optional term for the removed singularity from the Cauchy P.V. integral. */ 
  linear_farstream_func linear_upstream_int; /**< Optional approximation to the contribution from far upstream. */
  linear_farstream_func linear_downstream_int; /**< Optional approximation to the contribution from far downstream. */
  tpgraphy_integrand_func topography_integrand; /**<  Expression for the integrand of the channel topography. */
} method_s_funcs;

/** Functions required by the Cauchy P.V. boundary integrals for geometry of topography calculations. */
typedef struct {
  transform_func transform_phi; /**< Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
  transform_func transform_dphi; /**< Transform (arbitrary) relating \f$\mathrm{d} \beta\f$ to \f$\mathrm{d} \phi\f$. */
  tpgraphy_theta_func theta_func; /**< The type of topography being used expressed using \f$\theta\f$. */
  tpgraphy_convolution_func convolution_s; /**< Expression for the integrand of the integral over the free-surface (non P.V.). */
  tpgraphy_integrand_func topography_integrand; /**< Expression for the integrand of the integral over the topography. (P.V.) */
  removed_singularity_b_func removed_singularity_b; /**< Optional term for the removed singularity from the Cauchy P.V. integral over topography. */ 
  linear_farstream_func linear_upstream_s; /**< Optional approximation to the contribution from far upstream. */
  linear_farstream_func linear_downstream_s; /**< Optional approximation to the contribution from far downstream. */
} method_b_funcs;

/// \cond EXCLUDE_SYMBOLS
typedef struct {
  method_b_funcs method_b; /**<  */
  vector theta_s; /**< */
  vector beta_s; /**< */
  vector beta_sub; /**< */
  vector phi_sub; /**< */
  vector dphi_sub; /**< */
  linear_params lin_ustream;
  linear_params lin_dstream;
  real homotopy_s;
  vector topography; /**< */
  gsl_integration_workspace** ws;
  real temp_phi;
} x_b_integrand_params;
/// \endcond

/** Generic entry point for Lua loader

  Registers the library with Lua, creates a GSL integration workspace for use by the
  boundary integral method functions, and sets up the garbage collection for the workspace.
 */
int luaopen_integrals_sbox(lua_State *L);

void integrals_sbox_handler(const char *reason, 
                            const char *file, 
                            int line, 
                            int gsl_errno);


/** Generalised trapezoidal method for calculating \f$\tau_s(\beta_{i+1/2})\f$ based on Cauchy P.V. boundary integral methods.

 */
vector tau_s_worker(const vector&             theta_s, /**< (n)-vector  the angle of the free surface at \f$\beta_{i}\f$. */
                    const vector&             phi_sub, /**< (m)-vector of 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                    const vector&            dphi_sub, /**< (m)-vector of prescribed value of \f$\phi'\f$ at \f$\phi^{[j]}\f$. */
                    const vector&                beta, /**< (n)-vector of \f$\beta_{i}\f$. */
                    const vector&            beta_sub, /**< (m)-vector of \f$\beta^{[j]}\f$ */
                    const vector&          topography, /**< parameters with mathematically describe topography */
                    const linear_params&  lin_ustream, /**< parameters for linearised region upstream. */
                    const linear_params&  lin_dstream, /**< parameters for linearised region downstream. */
                    const real                      s, /**< optional homotopy parameter. */
                    gsl_integration_workspace**    ws, /**< integration workspace for linearised regions. */
                    method_s_funcs&            method  /**< current method functions. */);


/** Generalised mid-point method for calculating \f$y_s(\beta_{i+1/2})\f$ based on \f$\theta\f$ and \f$\tau\f$.

  \f[ W^{j} = \exp(-\tau_s(\beta_j)) \sin(\theta_s(\beta_j)) \f]
  \f[ y_s(\beta_{j-1/2}) = y_s(\beta_{j+1/2}) - \frac{1}{2}
                           \left( (\beta_{j+1}-\beta_{j}) W^{j+1/2} +
                                  (\beta_{j}-\beta_{j-1}) W^{j-1/2} \right). \f]

 */
vector y_s_worker(const vector&       theta_mid, /**< \f$(n-1)\f$-vector; \f$\theta\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector&         tau_mid, /**< \f$(n-1)\f$-vector; \f$\tau\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector&            beta, /**< \f$(n)\f$-vector; \f$\beta_{j}\f$. */
                  const vector&        beta_sub, /**< \f$(m)\f$-vector; value of \f$\beta\f$ at \f$\phi^{[j]}\f$. */
                  const vector&         phi_sub, /**< \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                  const vector&        dphi_sub, /**< \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$. */
                  const real                  s, /**< optional homotopy parameter. */
                  transform_func transform_dphi, /**< function to calculate \f$\phi'\f$ from \f$\beta\f$. */
                  real&                    eta0  /**< return value of \f$y\f$ at end-point of integral \f$\beta^{[0]}\f$. */);

/** Generalised mid-point method for calculating \f$x_s(\beta_{i+1/2})\f$ based on \f$\theta\f$ and \f$\tau\f$.

  Uses the mid-point integration rule to evaluate \f$x(\beta_{j+1/2})\f$, by
  averaging the value at \f$j\f$ and \f$j+1\f$.
  
  \f[ W^{j} = \exp(-\tau_s(\beta_j)) \cos(\theta_s(\beta_j)) \f]
  \f[ x_s(\beta_{j-1/2}) = x_s(\beta_{j+1/2}) - \frac{1}{2}
                           \left( (\beta_{j+1}-\beta_{j}) W^{j+1/2} +
                                  (\beta_{j}-\beta_{j-1}) W^{j-1/2} \right). \f]

 */
vector x_s_worker(const vector&       theta_mid, /**< \f$(n-1)\f$-vector; \f$\theta\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector&         tau_mid, /**< \f$(n-1)\f$-vector; \f$\tau\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector&            beta, /**< \f$(n)\f$-vector; \f$\beta_{j}\f$. */
                  const vector&        beta_sub, /**< \f$(m)\f$-vector; value of \f$\beta\f$ at \f$\phi^{[j]}\f$. */
                  const vector&         phi_sub, /**< \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                  const vector&        dphi_sub, /**< \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$. */
                  const real                  s, /**< optional homotopy parameter. */
                  transform_func transform_dphi, /**< function to calculate \f$\phi'\f$ from \f$\beta\f$. */
                  real&                      x0  /**< return value of \f$x\f$ at end-point of integral \f$\beta^{[0]}\f$. */);

/** \todo Document x_b_integrand */
real x_b_integrand(real, void*);

/** \todo Document x_b_worker */
vector x_b_worker(const vector&                 phi,
                  const vector&             theta_s,
                  const vector&              beta_s,
                  const vector&            beta_sub,
                  const vector&             phi_sub,
                  const vector&            dphi_sub,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  method_b_funcs             method,
                  const vector&          topography,
                  const linear_params&  lin_ustream,
                  const linear_params&  lin_dstream);

/** \todo Document y_b_integrand */
real y_b_integrand(real, void*);

/** \todo Document x_b_worker */
vector y_b_worker(const vector&                 phi,
                  const vector&             theta_s,
                  const vector&              beta_s,
                  const vector&            beta_sub,
                  const vector&             phi_sub,
                  const vector&            dphi_sub,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  method_b_funcs             method,
                  const vector&          topography,
                  const linear_params&  lin_ustream,
                  const linear_params&  lin_dstream);
