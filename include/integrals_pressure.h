/** \file integrals_pressure.h 
  
    Common functions used by all methods of calculating flows with pressure
    disturbances.

    See Chapter 3 of thesis "Very Steep Forced Solitary Waves in Two-Dimensional
    Free Surface Flow".
  
*/

#ifndef INTEGRALS_PRESSURE_H
#define INTEGRALS_PRESSURE_H

#include <stdlib.h> // For exit() and EXIT_FAILURE

extern "C" {
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

#include "pwpoly.h"
#include "integrals_types.h"
#include <gsl/gsl_integration.h>

#ifndef pi
#define pi M_PI /**< Constant value of \f$\pi\f$. */
#endif


typedef enum {
    USE_U,
    USE_V
} use_u_or_v;


/** \internal Approximation to a contribution from farstream to a Cauchy P.V.
              integral. */
typedef void (*farstream_integral_func)(const real_vector &, const real,
                                                 const real, const real,
                                              real_vector &);

/** Functions required by the Cauchy P.V. boundary integrals for free surface
    calculations with pressure disturbances. 
    
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
*/
typedef struct {
    transform_func transform_phi;
    transform_func transform_dphi;
    convolution_func convolution_terms;
    removed_singularity_func removed_singularity;
    farstream_integral_func linear_upstream_int;
    farstream_integral_func linear_downstream_int;
} method_funcs;


/** Generic entry point for Lua loader

    Registers the library with Lua, creates a GSL integration workspace for use
    by the boundary integral method functions, and sets up the garbage
    collection for the workspace.
 */
int luaopen_integrals_pressure(lua_State *L);


/** Generalised trapezoidal method for calculating \f$\u_s(\beta_{i+1/2})\f$
    based on Cauchy P.V. boundary integral methods.

    \param v_s (n)-real_vector the horizontal component of the free-surface
               velocity at \f$\beta_{i}\f$.
    \param phi_sub (m)-real_vector of 'clustering' locations in computational
                   grid; \f$\phi^{[j]}\f$.
    \param dphi_sub (m)-real_vector of prescribed value of \f$\phi'\f$ at
                    \f$\phi^{[j]}\f$.
    \param beta (n)-real_vector of \f$\beta_{i}\f$.
    \param beta_sub (m)-real_vector of \f$\beta^{[j]}\f$.
    \param topography parameters with mathematically describe topography.
    \param lin_ustream parameters for linearised region upstream.
    \param lin_dstream parameters for linearised region downstream.
    \param s optional homotopy parameter.
    \param ws integration workspace for linearised regions.
    \param method current method functions.
    \param u_s (n-1)-vector to store returned \f$u_s(\beta_{i+1/2})\f$.
*/
void u_s_worker(const real_vector &             v_s,
                  const real_vector &       phi_sub, 
                  const real_vector &      dphi_sub,
                  const real_vector &          beta,
                  const real_vector &      beta_sub,
                  const ff_params &     lin_ustream,
                  const ff_params &     lin_dstream,
                  const real                      s,
                  gsl_integration_workspace**    ws,
                  method_funcs&              method,
                  real_vector &                 u_s);


/** Generalised mid-point method for calculating \f$x_s(\beta_{i+1/2})\f$ r
    \f$y_s(\beta_{i+1/2})\f$ based on \f$\theta\f$ and \f$\tau\f$.

    Uses the mid-point integration rule to evaluate \f$x(\beta_{j+1/2})\f$, by
    averaging the value at \f$j\f$ and \f$j+1\f$. Whether to calculate \f$x_s\f$
    or \f$y_s\f$ is determine by using the cosine or sine function. 
 
    For \f$x_s\f$ let
    \f[ W^{j} = \frac{u_s(\beta_j)}{(u_s(\beta_j))^2 + (v_s(\beta_j))^2} \f]
    otherwise for \f$y_s\f$
    \f[ W^{j} = \frac{v_s(\beta_j)}{(u_s(\beta_j))^2 + (v_s(\beta_j))^2} \f]
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
void z_s_worker(const real_vector &      u_mid,
                const real_vector &      v_mid, 
                const real_vector &       beta, 
                const real_vector &   beta_sub, 
                const real_vector &    phi_sub, 
                const real_vector &   dphi_sub,
                const real                   s,
                transform_func  transform_dphi, 
                const use_u_or_v        u_or_v,
                real &                      z0,
                real_vector &              z_s);


#endif

