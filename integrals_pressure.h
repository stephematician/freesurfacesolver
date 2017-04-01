/** \file integrals_pressure.h 
  
  Common functions used by all methods of calculating flows with box shape disturbances. 
  
*/

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

/** \internal Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
typedef void (*transform_func)(const vector &, const vector &,
                               const vector &, const vector &,
                                     vector &, const real);
/** \internal Expression for the integrand of the Cauchy P.V. integral. */ 
typedef void (*convolution_func)(const int,      const vector &,
                                 const vector &, const vector &,
                                 const vector &,       vector &);
/** \internal Approximation to a contribution from farstream to a Cauchy P.V. integral. */
typedef void (*linear_farstream_func)(const vector &, const real,
                                      const real, const real, vector &);
/** \internal Expression for a removed singularity from a Cauchy P.V. integral. */ 
typedef void (*removed_singularity_func)(const vector&, const vector&, const vector&, vector&);

/** Typical functions required by the Cauchy P.V. boundary integrals for box calculations.


 */
typedef struct {
  transform_func transform_phi; /**< Transform (arbitrary) relating \f$\beta\f$ to \f$\phi\f$. */
  transform_func transform_dphi; /**< Transform (arbitrary) relating \f$\mathrm{d} \beta\f$ to \f$\mathrm{d} \phi\f$. */
  convolution_func convolution_terms; /**< Expresison for the integrand of the Cauchy P.V. integral. */
  removed_singularity_func removed_singularity; /**< Optional term for the removed singularity from the Cauchy P.V. integral. */ 
  linear_farstream_func linear_upstream_int; /**< Optional approximation to the contribution from far upstream. */
  linear_farstream_func linear_downstream_int; /**< Optional approximation to the contribution from far downstream. */
} method_funcs;

/** Parameters required by linearised approximations to far-stream contributions.

 
 */
typedef struct {
  real lambda; /**< Decay rate of velocity. */
  real D; /**< Vertical velocity constant. */
  real phi_match; /**< Location of farstream matching. */
} linear_params;

/** Generic entry point for Lua loader

  Registers the library with Lua, creates a GSL integration workspace for use by the
  boundary integral method functions, and sets up the garbage collection for the workspace.
 */
int luaopen_integrals_pressure(lua_State *L);

/** Generalised trapezoidal method for calculating \f$\tau_s\f$ based on Cauchy P.V. boundary integral methods.

 */
vector u_s_worker(const vector&                u_s, /**< (n)-vector of theta on free surface beta \f$\beta_{i}\f$. */
                  const vector&            phi_sub, /**< (m)-vector of 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                  const vector&           dphi_sub, /**< (m)-vector of prescribed value of $\phi'$ at \f$\phi^{[j]}\f$. */
                  const vector&               beta, /**< (n)-vector of */
                  const vector&           beta_sub, /**< (m)-vector */
                  const linear_params& lin_ustream, /**< parameters for linearised region upstream. */
                  const linear_params& lin_dstream, /**< parameters for linearised region downstream. */
                  const real                     s, /**< optional homotopy parameter. */
                  gsl_integration_workspace**   ws, /**< integration workspace for linearised regions. */
                  method_funcs&             method  /**< current method functions. */);

/** Generalised mid-point method for calculating \f$y_s(\beta_{i+1/2})\f$ based on \f$u\f$ and \f$v\f$.

  \f[ W^{j} = \frac{v_s(\beta_j)}{(u_s(\beta_j))^2 + (v_s(\beta_j))^2} \f]
  \f[ y_s(\beta_{j-1/2}) = y_s(\beta_{j+1/2}) - \frac{1}{2}
                           \left( (\beta_{j+1}-\beta_{j}) W^{j+1/2} +
                                  (\beta_{j}-\beta_{j-1}) W^{j-1/2} \right). \f]
 */
vector y_s_worker(const vector           &u_mid, /**< \f$(n-1)\f$-vector; \f$u\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector           &v_mid, /**< \f$(n-1)\f$-vector; \f$v\f$ on free surface at \f$\beta_{j+1/2}\f$. */
                  const vector            &beta, /**< \f$(n)\f$-vector; \f$\beta_{j}\f$. */
                  const vector        &beta_sub, /**< \f$(m)\f$-vector; value of \f$\beta\f$ at \f$\phi^{[j]}\f$. */
                  const vector         &phi_sub, /**< \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                  const vector        &dphi_sub, /**< \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$. */
                  const real                  s, /**< optional homotopy parameter. */
                  transform_func transform_dphi, /**< function to calculate \f$\phi'\f$ from \f$\beta\f$. */
                  real                    &eta0  /**< return value of \f$y\f$ at end-point of integral \f$\beta^{[0]}\f$. */);

/** Generalised mid-point method for calculating \f$x_s(\beta_{i+1/2})\f$ using \f$u\f$ and \f$v\f$.

  Uses the mid-point integration rule to evaluate
  \f$x(\beta_{j+1/2})\f$, by averaging the value at \f$j\f$ and \f$j+1\f$.
  
  \f[ W^{j} = \frac{u_s(\beta_j)}{(u_s(\beta_j))^2 + (v_s(\beta_j))^2} \f]
  \f[ x_s(\beta_{j-1/2}) = x_s(\beta_{j+1/2}) - \frac{1}{2}
                           \left( (\beta_{j+1}-\beta_{j}) W^{j+1/2} +
                                  (\beta_{j}-\beta_{j-1}) W^{j-1/2} \right). \f]

 */
vector x_s_worker(const vector           &u_mid, /**< \f$(n-1)\f$-vector; \f$u\f$ on free surface at \f$\beta_{i+1/2}\f$. */
                  const vector           &v_mid, /**< \f$(n-1)\f$-vector; \f$v\f$ on free surface at \f$\beta_{i+1/2}\f$. */
                  const vector            &beta, /**< \f$(n)\f$-vector; \f$\beta_{j}\f$. */
                  const vector        &beta_sub, /**< \f$(m)\f$-vector; value of \f$\beta\f$ at \f$\phi^{[j]}\f$. */
                  const vector         &phi_sub, /**< \f$(m)\f$-vector; 'clustering' locations in computational grid; \f$\phi^{[j]}\f$. */
                  const vector        &dphi_sub, /**< \f$(m)\f$-vector; value of \f$\phi'\f$ at \f$\phi^{[j]}\f$. */
                  const real                  s, /**< optional homotopy parameter. */
                  transform_func transform_dphi, /**< function to calculate \f$\phi'\f$ from \f$\beta\f$. */
                  real                      &x0  /**< return value of \f$x\f$ at end-point of integral \f$\beta^{[0]}\f$. */);
